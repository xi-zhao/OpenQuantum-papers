"""Predict fixed-preparation detuning response without fitting scan curves.

Parameter tangents propagate under an augmented full PI generator. The small-g
bound concerns only replacing two exponentials by their matching repeated-pole
shape; it is not a physical alternative generator or a metrological bound.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.sparse import bmat
from scipy.sparse.linalg import expm_multiply
from biorthogonal_generator import generators
from coherence_response import real_checked
from physical_effects import ROOT, DATA, TIMES, preparations


def tangent_prediction(scan):
  N, detuning, width = scan['N'], scan['crossing'], scan['mixing_width']
  L0, derivative, space = generators(N)
  L = (L0+detuning*derivative).tocsc()
  source, _ = preparations(scan['preparations'], space)
  # y0'=Ly0; y1'=Ly1+L1 y0; y2'=Ly2+2 L1 y1.
  augmented = bmat([[L, None, None], [derivative, L, None],
    [None, 2*derivative, L]], format='csc')
  initial = np.vstack([source, np.zeros_like(source), np.zeros_like(source)])
  values = expm_multiply(augmented, initial, start=0., stop=TIMES[-1], num=len(TIMES),
    traceA=3*L.diagonal().sum())
  d = len(source)
  observable = space['trace']*space['ne'][space['row']]/N
  signal = [real_checked(np.einsum('i,tij->tj', observable, values[:, i*d:(i+1)*d],
    optimize=False), 'parameter tangent') for i in range(3)]
  central = next(r for r in scan['records'] if r['offset'] == 0.)
  center_error = float(np.max(abs(signal[0]-np.asarray(central['full']))))
  if center_error > 1e-10:
    raise RuntimeError(f'Augmented propagation failed at center: {center_error}')
  points = []
  for offset in (-width, width):
    record = min(scan['records'], key=lambda r: abs(r['offset']-offset))
    actual = np.asarray(record['full'])
    linear = signal[0]+offset*signal[1]
    quadratic = linear+offset**2/2*signal[2]
    denom = np.linalg.norm(actual, axis=0)
    points.append(dict(offset=offset,
      linear_maximum_absolute_error=np.max(abs(linear-actual), axis=0).tolist(),
      quadratic_maximum_absolute_error=np.max(abs(quadratic-actual), axis=0).tolist(),
      linear_relative_l2=(np.linalg.norm(linear-actual, axis=0)/denom).tolist(),
      quadratic_relative_l2=(np.linalg.norm(quadratic-actual, axis=0)/denom).tolist()))
  modes = central['modes']
  kappa, g = -modes['mean_rate'], modes['g']
  A, B = np.asarray(modes['A']), np.asarray(modes['B'])
  if kappa <= g or g <= 0:
    raise RuntimeError('Uniform splitting bound requires kappa > g > 0')
  bound = g**2*(2*abs(A)/(np.e**2*(kappa-g)**2)+
    4.5*abs(B)/(np.e**3*(kappa-g)**3))
  measured = np.asarray(central['finite_splitting_effect']['maximum_absolute_change'])
  if np.any(measured > bound):
    raise RuntimeError('Sampled two-pole difference exceeds analytic bound')
  n, u, v = scan['preparations']['fold']
  exact_second = 1.5*2*v-2*(detuning-10*(N-1)/N*n)*u
  actual_second = float(real_checked(np.dot(observable, L@(L@source[:, 0])), 'second derivative'))
  if abs(exact_second-actual_second) > 1e-10:
    raise RuntimeError('Exact initial second derivative failed')
  return dict(N=N, mixing_width=width, central_reproduction_error=center_error,
    first_derivative=signal[1].tolist(), second_derivative=signal[2].tolist(),
    maximum_first_derivative=np.max(abs(signal[1]), axis=0).tolist(),
    maximum_second_derivative=np.max(abs(signal[2]), axis=0).tolist(), predictions=points,
    splitting_all_time_bound=bound.tolist(), splitting_sampled_difference=measured.tolist(),
    initial_second_derivative=dict(analytic=exact_second, full_generator=actual_second),
    bound_scope='All t>=0 for this pair with fixed A,B and mean rate; no claim about all states, '
      'measurement noise, physical realizability of the coalesced shape, or the remaining spectrum.')


def main(args):
  path, output = DATA/args.input, DATA/args.output
  if output.exists():
    raise FileExistsError(f'Preserve evidence; choose a new output: {output}')
  data = json.loads(path.read_text())
  for p, sha in data['source_sha256'].items():
    if hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != sha:
      raise RuntimeError(f'Source changed: {p}')
  records = []
  for scan in data['scans']:
    record = tangent_prediction(scan)
    records.append(record)
    print(json.dumps({k: record[k] for k in ('N', 'maximum_first_derivative',
      'maximum_second_derivative', 'predictions', 'splitting_all_time_bound')}), flush=True)
  result = dict(created='2026-11-01', input=str(path.relative_to(ROOT)),
    input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), records=records,
    interpretation='No-fit full-generator parameter derivatives at the crossing; finite sizes only.')
  output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--input', default='effects_2026-11-01.json')
  parser.add_argument('--output', default='sensitivity_2026-11-01.json')
  main(parser.parse_args())
