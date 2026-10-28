"""Exact Schur continued fraction in signed coherence order.

For real spectral parameter z, negative-q blocks are conjugate to positive-q
blocks. The diagonal-sector nonlinear operator is K(z)=A0+2 Re Sigma0(z).
Roots are computed without the exact full-L eigenvectors or gaps. The exact
crossing detuning remains an input. This is a numerical spectral reduction,
not an independent large-N action prediction.
"""
from pathlib import Path
import argparse
import json
import time
import warnings
import numpy as np
from scipy.linalg import solve,eigvals
from scipy.optimize import minimize_scalar,brentq
from biorthogonal_generator import generators
from biorthogonal_audit import mm

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/coherence_mechanism'


def checked_root(residual,left,right,tolerance=1e-10):
  """Reject false Brent convergence at a nearest-branch discontinuity."""
  root=brentq(residual,left,right,xtol=2e-13)
  error=float(residual(root))
  if not np.isfinite(error) or abs(error)>tolerance:
    raise RuntimeError(f'Root failed residual check: z={root}, f={error}')
  return root


def nearest_real_index(values,z):
  """Use only a real LAPACK eigenvalue; an unresolved complex pair is rejected."""
  index=int(np.argmin(abs(values-z)))
  if values[index].imag!=0:
    raise RuntimeError(f'Nearest effective eigenvalue is complex: {values[index]}')
  return index


class Fraction:
  def __init__(self,N,detuning,qmax=None):
    if not np.isreal(detuning):
      raise ValueError('Conjugate-sector reduction requires real detuning')
    self.qmax=N if qmax is None else qmax
    if not isinstance(self.qmax,(int,np.integer)) or not 0<=self.qmax<=N:
      raise ValueError('qmax must be an integer between zero and N')
    L0,L1,s=generators(N);L=L0+float(np.real(detuning))*L1
    q=s['ne'][s['row']]-s['ne'][s['col']]
    ids=[np.flatnonzero(q==k) for k in range(self.qmax+1)]
    self.A=[L[v][:,v].toarray() for v in ids]
    self.up=[L[ids[k]][:,ids[k+1]].toarray() for k in range(self.qmax)]
    self.down=[L[ids[k+1]][:,ids[k]].toarray() for k in range(self.qmax)]
    self.N=N
    self.cache={}
    self.samples=[]

  def operator(self,z,depth=None):
    if not np.isreal(z):
      raise ValueError('The 2 Re Sigma formula is restricted to real z')
    top=self.qmax if depth is None else depth
    if not isinstance(top,(int,np.integer)) or not 0<=top<=self.qmax:
      raise ValueError('depth must be an available nonnegative coherence order')
    sigma=np.zeros_like(self.A[top])
    for k in range(top,0,-1):
      B=z*np.eye(len(self.A[k]))-self.A[k]-sigma
      sigma=mm(self.up[k-1],solve(B,self.down[k-1],assume_a='gen'))
    return self.A[0].real+2*sigma.real

  def residual(self,z):
    key=float(z)
    if key not in self.cache:
      w=eigvals(self.operator(z))
      ix=nearest_real_index(w,z); value=w[ix]
      self.cache[key]=float(value.real-z)
      self.samples.append([key,self.cache[key]])
    return self.cache[key]


def main(args):
  start=time.monotonic()
  ref=json.loads((ROOT/f'data/biorthogonal_audit/N{args.N}_float64_s17_t-0.76.json').read_text())
  d=ref['detuning']
  model=Fraction(args.N,d,args.qmax)
  with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    result=minimize_scalar(model.residual,bounds=(-.81,-.71),method='bounded',
                           options={'xatol':2e-11,'maxiter':45})
    center=float(result.x);value=float(result.fun)
    if not result.success or value>=0:
      raise RuntimeError(f'No resolved root pair: {result}')
    brackets=[]
    for sign in (-1,1):
      step=1e-5
      for _ in range(16):
        endpoint=center+sign*step
        if model.residual(endpoint)>0:
          brackets.append(tuple(sorted([center,endpoint])))
          break
        step*=1.8
      else:
        raise RuntimeError('Could not bracket a real root near the minimum')
    left=checked_root(model.residual,*brackets[0])
    right=checked_root(model.residual,*brackets[1])
    frozen=model.operator(center)
    w=eigvals(frozen); selected=w[np.argsort(abs(w-center))[:4]]
    fd={}
    for dz in (.01,.003,.001):
      f0=model.residual(center)
      fp=model.residual(center+dz);fm=model.residual(center-dz)
      curvature=(fp+fm-2*f0)/dz**2
      fd[str(dz)]=dict(curvature=curvature,quadratic_g=float(np.sqrt(-2*f0/curvature)) if curvature>0 else None)
  DATA.mkdir(exist_ok=True)
  out=dict(N=args.N,detuning=d,qmax=model.qmax,center=center,minimum_residual=value,
    roots=[left,right],g=(right-left)/2,reference_g=ref['g'],
    reference_relative_error=float((right-left)/(2*ref['g'])-1),
    frozen_eigenvalues=[[float(z.real),float(z.imag)] for z in selected],
    local_quadratic=fd,samples=model.samples,
    warnings=sorted(set(str(w.message) for w in caught)),seconds=time.monotonic()-start)
  (DATA/f'fraction_N{args.N}_q{model.qmax}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
  print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':
  parser=argparse.ArgumentParser(description=__doc__)
  parser.add_argument('N',type=int)
  parser.add_argument('--qmax',type=int)
  main(parser.parse_args())
