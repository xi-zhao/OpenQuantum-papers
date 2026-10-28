"""Physical initial-state projections and mode continuation for the PRA draft.

The stored Dicke blocks sum representation multiplicities. Product states are
constructed as multiplicity_j * det(rho)^(N/2) * exp(2 eta rhat.J_j).
Residues in an observable are signed and are not modal probabilities.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import time
import numpy as np
from scipy.linalg import eigh, eig
from scipy.sparse.linalg import eigs, expm_multiply
from biorthogonal_generator import generators, populations
from biorthogonal_audit import project, mm, hermitize
from mf import roots

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/pra_spectral'
OMEGA, V = 2., 10.
DF, NF = 3.599636890457168, .15692704039732738


def fixed(detuning, n):
  c = (-1j*OMEGA/2*(1-2*n))/(.5-1j*(detuning-V*n))
  return np.array([n, c.real, c.imag])


def jacobian(detuning, state):
  n, u, v = state
  d = detuning-V*n
  return np.array([[-1., 0., -OMEGA], [V*v, -.5, -d],
                   [OMEGA-V*u, d, -.5]])


def product_state(state, space):
  n, u, v = state
  bloch = np.array([2*u, -2*v, 2*n-1])
  length = np.linalg.norm(bloch)
  if length >= 1:
    raise ValueError('Product-state seed must be strictly positive')
  logdet = np.log((1-length**2)/4)
  eta_over_r = np.arctanh(length)/length if length else 1.
  N, d = space['N'], space['d']
  full = np.zeros((d, d), complex)
  for tj in range(N, -1, -2):
    j = tj/2
    m = np.arange(j, -j-1, -1)
    jp = np.diag(np.sqrt((j-m[1:])*(j+m[1:]+1)), 1)
    jx = (jp+jp.T)/2
    jy = (jp-jp.T)/(2j)
    h = bloch[0]*jx+bloch[1]*jy+bloch[2]*np.diag(m)
    w, vectors = eigh(h)
    k = (N-tj)//2
    multiplicity = math.comb(N, k)-(math.comb(N, k-1) if k else 0)
    weights = np.exp(np.log(multiplicity)+N/2*logdet+2*eta_over_r*w)
    block = mm(vectors*weights[None, :], vectors.conj().T)
    ix = np.flatnonzero(space['twoj'] == tj)
    full[np.ix_(ix, ix)] = block
  vector = full[space['row'], space['col']]
  trace = np.dot(space['trace'], vector)
  if abs(trace-1) > 1e-10:
    raise RuntimeError(f'Product state trace mismatch: {trace}')
  return vector/trace


def product_validation():
  records = []
  for N in (1, 2, 4):
    L0, L1, s = generators(N)
    state = np.array([.29, .13, -.18])
    p = product_state(state, s)
    pop = populations(p, s)
    expected = np.array([math.comb(N, k)*state[0]**k*(1-state[0])**(N-k)
                         for k in range(N+1)])
    full = np.zeros((s['d'], s['d']), complex)
    full[s['row'], s['col']] = p
    jz = s['ne']-N/2
    jp = np.zeros_like(full)
    for tj in range(N, -1, -2):
      ids = np.flatnonzero(s['twoj'] == tj)
      j = tj/2; m = np.arange(j, -j-1, -1)
      jp[np.ix_(ids, ids)] = np.diag(np.sqrt((j-m[1:])*(j+m[1:]+1)), 1)
    jy = (jp-jp.T)/(2j)
    population_error = float(np.max(abs(pop-expected)))
    u = float(np.trace(mm(full, s['jx'].toarray())).real/N)
    v = float(-np.trace(mm(full, jy)).real/N)
    coherence_error = max(abs(u-state[1]), abs(v-state[2]))
    if max(population_error, coherence_error) > 1e-12:
      raise RuntimeError('Product-state convention check failed')
    records.append(dict(N=N, population_error=population_error,
                        coherence_error=coherence_error))
  return records


def oscillatory_pair(L, s, seed):
  rng = np.random.default_rng(seed)
  shift = -.72+2.85j
  values, vectors = eigs(L, k=4, sigma=shift, tol=2e-14,
    v0=rng.normal(size=L.shape[0])+1j*rng.normal(size=L.shape[0]))
  left_values, left_vectors = eigs(L.T, k=4, sigma=shift, tol=2e-14,
    v0=rng.normal(size=L.shape[0])+1j*rng.normal(size=L.shape[0]))
  index = int(np.argmin(abs(values-shift)))
  value = values[index]
  if value.imag < 1.:
    raise RuntimeError('Could not identify the oscillatory pair')
  left_index = int(np.argmin(abs(left_values-value)))
  right = vectors[:, index]
  left = left_vectors[:, left_index]
  R = np.column_stack([right, right[s['swap']].conj()])
  W = np.column_stack([left, left[s['swap']].conj()])
  W = mm(W, np.linalg.inv(mm(W.T, R)).T)
  lam = np.array([value, value.conjugate()])
  residual = max(float(np.linalg.norm(L@R[:, i]-lam[i]*R[:, i])/np.linalg.norm(R[:, i]))
                 for i in range(2))
  left_residual = max(float(np.linalg.norm(L.T@W[:, i]-lam[i]*W[:, i])/np.linalg.norm(W[:, i]))
                      for i in range(2))
  if max(residual, left_residual) > 1e-9:
    raise RuntimeError('Oscillatory mode residual failed')
  return lam, R, W, dict(right_residual=residual, left_residual=left_residual,
    candidates=[[float(z.real), float(z.imag)] for z in values])


def run(N, dynamics=False, quick=False):
  start = time.monotonic()
  DATA.mkdir(exist_ok=True)
  reference = json.loads((ROOT/f'data/biorthogonal_audit/N{N}_float64_s17_t-0.76.json').read_text())
  star = reference['detuning']
  L0, L1, s = generators(N)
  diagonal = s['trace'].astype(bool)
  obs = s['trace']*s['ne'][s['row']]/N
  qfold = fixed(DF, NF)
  fold_seed = product_state(qfold, s)
  offsets = [-.06, 0., .06] if quick else np.unique(np.r_[np.linspace(-.06, .06, 13), -.003, .003])
  records = []
  for number, offset in enumerate(offsets):
    d = star+offset
    audit, Rd, Wd, L = project(L0, L1, s, d, 701+N+number, -.76)
    h = np.array(audit['H'])
    lam, transform = eig(h)
    order = np.argsort(lam.real)
    lam = lam[order]
    transform = transform[:, order]
    if max(abs(lam.imag)) > 1e-10:
      raise RuntimeError('Target branch is not real')
    R = mm(Rd, transform)
    W = mm(Wd, np.linalg.inv(transform).T)
    n_hi = roots(d, OMEGA, V)[-1]
    qhi = fixed(d, n_hi)
    j = jacobian(d, qhi)
    mfvalues, mfvectors = eig(j)
    realidx = np.flatnonzero(abs(mfvalues.imag) < 1e-10)
    slow = realidx[np.argmax(mfvalues[realidx].real)]
    direction = mfvectors[:, slow].real
    direction *= np.sign(direction[0])/np.linalg.norm(direction)
    step = .005
    high_plus = product_state(qhi+step*direction, s)
    high_minus = product_state(qhi-step*direction, s)
    high_delta = high_plus-high_minus
    w0, ss0 = eigs(L, k=1, sigma=1e-9, tol=2e-14,
                   v0=np.random.default_rng(990+N).normal(size=L.shape[0]))
    ss, _ = hermitize(ss0[:, 0], s)
    ss /= np.dot(s['trace'], ss)
    ss_residual = float(np.linalg.norm(L@ss)/np.linalg.norm(ss))
    folded_delta = fold_seed-ss
    inputs = np.column_stack([high_delta, folded_delta])
    coefficients = mm(W.T, inputs)
    visibility = np.dot(obs, R)
    residues = visibility[:, None]*coefficients
    if max(abs(residues.imag).ravel()) > 1e-7:
      raise RuntimeError('Physical residues are not real')
    normpop = []
    location = []
    for i in range(2):
      profile = populations(R[:, i], s)
      profile *= np.exp(-1j*np.angle(profile[np.argmax(abs(profile))]))
      if max(abs(profile.imag)) > 1e-8:
        raise RuntimeError('Hermitian profile convention failed')
      profile = profile.real/max(abs(profile))
      normpop.append(profile.tolist())
      weights = abs(profile)
      location.append(dict(absolute_mean=float(np.dot(np.arange(N+1)/N, weights)/sum(weights)),
        peak=float(np.argmax(weights)/N),
        low_region_absolute_fraction=float(sum(weights[np.arange(N+1)/N < .22])/sum(weights))))
    record = dict(N=N, detuning=float(d), offset=float(offset), eigenvalues=lam.real.tolist(),
      high_meanfield_state=qhi.tolist(), high_meanfield_rate=float(mfvalues[slow].real),
      meanfield_all_eigenvalues=[[float(z.real), float(z.imag)] for z in mfvalues],
      eigenmode_population_profiles=normpop, eigenmode_population_location=location,
      diabatic_slopes=np.diag(audit['H1']).tolist(), H=audit['H'],
      observable='excitation density <Ne>/N',
      high_perturbation_step=step, high_perturbation_direction=direction.tolist(),
      residues=residues.real.tolist(),
      initial_signal=np.dot(obs, inputs).real.tolist(), stationary_density=float(np.dot(obs, ss).real),
      steady_state_residual=ss_residual, diagnostics=audit['diagnostics'])
    if ss_residual > 1e-9:
      raise RuntimeError('Stationary reference failed')
    if dynamics and (abs(offset)<1e-12 or abs(abs(offset)-.06)<1e-12):
      times = np.linspace(0., 12., 121)
      evolved = expm_multiply(L, inputs, start=0., stop=12., num=len(times),
                             traceA=L.diagonal().sum())
      signal = np.einsum('i,tij->tj', obs, evolved, optimize=False).real
      pair = np.einsum('ti,ij->tj', np.exp(times[:, None]*lam.real[None, :]), residues.real,
                       optimize=False)
      oscillatory_lam, Ro, Wo, oscillatory_checks = oscillatory_pair(L, s, 201+N+number)
      oscillatory_residues = np.dot(obs, Ro)[:, None]*mm(Wo.T, inputs)
      oscillatory_signal = np.einsum('ti,ij->tj', np.exp(times[:, None]*oscillatory_lam[None, :]),
        oscillatory_residues, optimize=False)
      if np.max(abs(oscillatory_signal.imag)) > 1e-9:
        raise RuntimeError('Oscillatory pair does not give a real signal')
      four = pair+oscillatory_signal.real
      windows = {}
      for lower in (0., 2., 4., 6.):
        mask = times >= lower
        windows[str(lower)] = dict(maximum_absolute_error=np.max(abs(signal[mask]-pair[mask]),axis=0).tolist(),
          relative_l2_error=(np.linalg.norm(signal[mask]-pair[mask],axis=0)/np.linalg.norm(signal[mask],axis=0)).tolist(),
          four_mode_maximum_absolute_error=np.max(abs(signal[mask]-four[mask]),axis=0).tolist(),
          four_mode_relative_l2_error=(np.linalg.norm(signal[mask]-four[mask],axis=0)/np.linalg.norm(signal[mask],axis=0)).tolist())
      record['dynamics'] = dict(times=times.tolist(), exact=signal.tolist(), two_modes=pair.tolist(),
        four_modes=four.tolist(), oscillatory_eigenvalues=[[float(z.real),float(z.imag)] for z in oscillatory_lam],
        oscillatory_residues=[[[float(z.real),float(z.imag)] for z in row] for row in oscillatory_residues],
        oscillatory_checks=oscillatory_checks,
        windows=windows, trace_error=float(np.max(abs(np.einsum('i,tij->tj',s['trace'],evolved,optimize=False)))))
    records.append(record)
    print(json.dumps({k:record[k] for k in ('N','offset','eigenvalues','high_meanfield_rate','residues')}),flush=True)
    output = dict(N=N, crossing=star, fold=dict(detuning=DF, state=qfold.tolist()),
      product_state_validation=product_validation() if len(records)==1 else output['product_state_validation'],
      records=records, seconds=time.monotonic()-start,
      interpretation='Signed spectral residues; no probabilistic mode weights or universal ghost identity are assumed.')
    (DATA/f'modes_N{N}{"_quick" if quick else ""}.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
  return output


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('N',type=int)
  parser.add_argument('--dynamics',action='store_true')
  parser.add_argument('--quick',action='store_true')
  args = parser.parse_args()
  run(args.N,args.dynamics,args.quick)
