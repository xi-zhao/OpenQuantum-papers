"""Audited local two-mode projection; keeps earlier exploratory data intact.

Left columns store the transpose dual: W.T R = I, without conjugation.
All dense contractions use einsum to avoid the local Accelerate matmul
floating-point warnings; sparse factorization is SciPy SuperLU.
"""
import argparse
import json
from pathlib import Path
import time
import warnings
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigs, splu, LinearOperator
from biorthogonal_generator import generators, populations

ROOT = Path(__file__).resolve().parents[1]


def mm(a, b):
  return np.einsum('ik,kj->ij', a, b, optimize=False)


def hermitize(v, s):
  diag = np.flatnonzero(s['trace'])
  pivot = diag[np.argmax(abs(v[diag]))]
  v = v * np.exp(-1j*np.angle(v[pivot]))
  defect = np.linalg.norm(v-v[s['swap']].conj())/np.linalg.norm(v)
  return (v+v[s['swap']].conj())/2, float(defect)


def project(L0, L1, s, detuning, seed, target):
  L = (L0+detuning*L1).tocsc()
  rng = np.random.default_rng(seed)
  lu = splu(L-target*sp.eye(L.shape[0], format='csc'))
  op = LinearOperator(L.shape, matvec=lu.solve, dtype=complex)
  opt = LinearOperator(L.shape, matvec=lambda x: lu.solve(x, trans='T'), dtype=complex)
  w, Rall = eigs(L, k=6, sigma=target, OPinv=op, tol=2e-14,
                 v0=rng.normal(size=L.shape[0])+1j*rng.normal(size=L.shape[0]))
  wl, Wall = eigs(L.T, k=6, sigma=target, OPinv=opt, tol=2e-14,
                  v0=rng.normal(size=L.shape[0])+1j*rng.normal(size=L.shape[0]))
  select = [i for i in np.argsort(abs(w-target)) if abs(w[i].imag) < 1e-7][:2]
  if len(select) != 2:
    raise ValueError('Did not find two real modes')
  select = sorted(select, key=lambda i: w[i].real)
  lam = w[select]
  match = [int(np.argmin(abs(wl-value))) for value in lam]
  if len(set(match)) != 2:
    raise ValueError('Left/right matching is not one-to-one')
  R, W = Rall[:, select].copy(), Wall[:, match].copy()
  res_r = [np.linalg.norm(L@R[:,i]-lam[i]*R[:,i])/np.linalg.norm(R[:,i]) for i in range(2)]
  res_l = [np.linalg.norm(L.T@W[:,i]-wl[match[i]]*W[:,i])/np.linalg.norm(W[:,i]) for i in range(2)]
  herm = []
  for i in range(2):
    R[:,i], dr = hermitize(R[:,i], s)
    W[:,i], dl = hermitize(W[:,i], s)
    herm.append([dr,dl])
  gram = mm(W.T,R)
  W = mm(W,np.linalg.inv(gram).T)
  A = mm(W.T,L@R)
  B = mm(W.T,L1@R)
  if np.max(abs(B.imag)) > 1e-8:
    raise ValueError('Hermitian-basis projected derivative is not real')
  slopes,T = np.linalg.eig(B.real)
  if np.max(abs(slopes.imag)) > 1e-8:
    raise ValueError('No real derivative eigenbasis')
  order = np.argsort(abs(slopes))
  T = T[:,order].real
  Rd = mm(R,T)
  Wd = mm(W,np.linalg.inv(T).T)
  scales = []
  for i in range(2):
    pop = populations(Rd[:,i],s).real
    scales.append(pop[np.argmax(abs(pop))])
  Rd /= np.array(scales)
  Wd *= np.array(scales)
  H = mm(Wd.T,L@Rd)
  H1 = mm(Wd.T,L1@Rd)
  H0 = mm(Wd.T,L0@Rd)
  product = H[0,1]*H[1,0]
  discriminant = ((H[0,0]-H[1,1])/2)**2+product
  eps = max(float(np.max(abs(H.imag))),float(np.max(abs(H1.imag))))
  diagnostics = dict(
    residual_right=res_r,residual_left=res_l,
    left_right_eigenvalue_mismatch=float(max(abs(lam-wl[match]))),
    hermiticity_defects=herm,gram_condition=float(np.linalg.cond(gram)),
    right_pair_condition=float(np.linalg.cond(R)),
    derivative_basis_condition=float(np.linalg.cond(T)),
    biorthogonality_error=float(np.linalg.norm(mm(Wd.T,Rd)-np.eye(2))),
    invariant_right_error=float(np.linalg.norm(L@Rd-mm(Rd,H))/np.linalg.norm(Rd)),
    invariant_left_error=float(np.linalg.norm(L.T@Wd-mm(Wd,H.T))/np.linalg.norm(Wd)),
    trace_preservation_defect=float(max(abs(s['trace']@L))),
    right_mode_trace=[float(abs(np.sum(s['trace']*Rd[:,i]))) for i in range(2)],
    projected_imaginary_max=eps,
    derivative_offdiagonal=float(abs(H1[0,1])+abs(H1[1,0])),
    halfgap_identity_error=float(abs(np.sqrt(discriminant)-(lam[1]-lam[0])/2)),
    diagonal_correction_fraction=float(abs(((H[0,0]-H[1,1])/2)**2/product)),
    all_finite=bool(all(np.isfinite(v).all() for v in (Rd,Wd,H,H1))))
  out = dict(detuning=detuning,lam=lam.real.tolist(),g=float((lam[1]-lam[0]).real/2),
             H=H.real.tolist(),H1=H1.real.tolist(),H0=H0.real.tolist(),
             product=[float(product.real),float(product.imag)],diagnostics=diagnostics)
  return out,Rd,Wd,L


def bins(values,labels,size):
  return (np.bincount(labels,weights=values.real,minlength=size)
          +1j*np.bincount(labels,weights=values.imag,minlength=size))


def encode(values):
  return [[float(z.real),float(z.imag)] for z in values]


def decomposition(L,L0,L1,R,W,s):
  N=s['N']; ne=s['ne']; row=s['row']; col=s['col']
  # Entry labels are a coordinate diagnostic, not a positive probability.
  coordinates=dict(excitation_sum=ne[row]+ne[col],
                   coherence_order=abs(ne[row]-ne[col]),twoj=s['twoj'][row])
  result={}
  for name,a,b in [('ab',0,1),('ba',1,0)]:
    term=W[:,a]*(L@R[:,b])
    net=np.sum(term)
    entry=L.tocoo()
    raw=W[entry.row,a]*entry.data*R[entry.col,b]
    shifted=W[:,a]*((L-sp.eye(L.shape[0])*0.5*
                      np.trace(mm(W.T,L@R)))@R[:,b])
    values=dict(net=[float(net.real),float(net.imag)],
                row_abs_sum=float(np.sum(abs(term))),
                edge_abs_sum=float(np.sum(abs(raw))),
                cancellation_ratio=float(np.sum(abs(raw))/abs(net)),
                row_cancellation_ratio=float(np.sum(abs(term))/abs(net)),
                shifted_row_abs_sum=float(np.sum(abs(shifted))),
                L0_element=encode([np.sum(W[:,a]*(L0@R[:,b]))])[0],
                L1_element=encode([np.sum(W[:,a]*(L1@R[:,b]))])[0])
    for label,ids in coordinates.items():
      signed=bins(term,ids,2*N+1)
      absolute=np.bincount(ids,weights=abs(term),minlength=2*N+1)
      values[label]=dict(signed=encode(signed),absolute=absolute.tolist())
    result[name]=values
  prof={}
  for name,vec in [('ra',R[:,0]),('rb',R[:,1]),('la',W[:,0]),('lb',W[:,1])]:
    prof[name]=dict(population=populations(vec,s).real.tolist(),
                   excitation_l1=np.bincount(coordinates['excitation_sum'],
                     weights=abs(vec),minlength=2*N+1).tolist(),
                   max_abs=float(max(abs(vec))),l2=float(np.linalg.norm(vec)))
  return result,prof


def run(args):
  start=time.monotonic()
  inputs=json.loads((ROOT/'data/biorthogonal_audit/inputs.json').read_text())
  old={int(r[0]):r for r in inputs['previous_crossings']}
  detuning=old[args.N][1]
  L0,L1,s=generators(args.N,backend=args.backend)
  iteration=[]
  with warnings.catch_warnings(record=True) as captured:
    warnings.simplefilter('always')
    for _ in range(4 if args.refine else 1):
      out,R,W,L=project(L0,L1,s,detuning,args.seed,args.target)
      H,H1=np.array(out['H']),np.array(out['H1'])
      correction=-(H[0,0]-H[1,1])/(H1[0,0]-H1[1,1])
      iteration.append(dict(detuning=detuning,correction=float(correction),g=out['g']))
      if abs(correction)<1e-11:
        break
      if args.refine:
        detuning+=correction
    if args.refine and abs(correction) >= 1e-11:
      raise RuntimeError('Crossing refinement did not converge')
    diagnostics=out['diagnostics']
    if (not diagnostics['all_finite'] or
        max(diagnostics['residual_right']+diagnostics['residual_left']) > 1e-9 or
        diagnostics['biorthogonality_error'] > 1e-8 or
        diagnostics['projected_imaginary_max'] > 1e-8):
      raise RuntimeError(f'Projection checks failed: {diagnostics}')
    out['decomposition'],out['profiles']=decomposition(L,L0,L1,R,W,s)
  out.update(N=args.N,omega=2.0,interaction=10.0,backend=args.backend,
             seed=args.seed,target=args.target,refine=args.refine,
             iterations=iteration,previous_g=old[args.N][2],
             warnings=sorted(set(str(w.message) for w in captured)),
             seconds=time.monotonic()-start)
  outdir=ROOT/'data/biorthogonal_audit'
  outdir.mkdir(exist_ok=True)
  tag=f'N{args.N}_{args.backend}_s{args.seed}_t{args.target:g}'
  (outdir/f'{tag}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
  if args.vectors:
    np.savez_compressed(outdir/f'{tag}.npz',R=R,W=W,detuning=out['detuning'])
  print(json.dumps({k:out[k] for k in ('N','backend','detuning','g','H','diagnostics','seconds')},indent=2),flush=True)
  print('cancellation',[(k,out['decomposition'][k]['cancellation_ratio']) for k in ('ab','ba')],flush=True)


if __name__=='__main__':
  p=argparse.ArgumentParser(description=__doc__)
  p.add_argument('N',type=int)
  p.add_argument('--backend',choices=['float64','piqs'],default='float64')
  p.add_argument('--seed',type=int,default=17)
  p.add_argument('--target',type=float,default=-0.76)
  p.add_argument('--refine',action='store_true')
  p.add_argument('--vectors',action='store_true')
  run(p.parse_args())
