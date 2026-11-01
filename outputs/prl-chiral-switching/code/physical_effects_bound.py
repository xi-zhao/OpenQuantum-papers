"""Check full-signal CPTP parameter bounds using frozen density records.

For fixed Hermitian x=rho_a-rho_b, D0=||x||_1/2 and 0<=O<=I,
|d_Delta s|<=N t D0 and |d_Delta^2 s|<=N^2 t^2 D0.
The commutator is bounded after centering Ne by N/2. These are ordinary
Duhamel/trace-contraction bounds, not a claim of a new theorem.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/physical_effects'


def main(args):
  path, output = DATA/args.input, DATA/args.output
  if output.exists():
    raise FileExistsError(f'Preserve evidence; choose another output: {output}')
  data = json.loads(path.read_text())
  records = []
  for scan in data['scans']:
    N, width, distance = scan['N'], scan['mixing_width'], scan['initial_trace_distance']
    center = next(r for r in scan['records'] if r['offset'] == 0.)
    lo = min(scan['records'], key=lambda r: abs(r['offset']+width))
    hi = min(scan['records'], key=lambda r: abs(r['offset']-width))
    endpoint, curvature = [0., 0.], [0., 0.]
    for i, t in enumerate(data['times']):
      for k in range(2):
        change = abs(hi['full'][i][k]-lo['full'][i][k])
        bend = abs(hi['full'][i][k]+lo['full'][i][k]-2*center['full'][i][k])
        if t == 0:
          if max(change, bend) > 1e-12:
            raise RuntimeError('A detuning scan changed the initial signal')
          continue
        endpoint[k] = max(endpoint[k], change/(2*N*width*t*distance[k]))
        curvature[k] = max(curvature[k], bend/((N*width*t)**2*distance[k]))
    if max(endpoint+curvature) > 1+1e-8:
      raise RuntimeError('Full-signal contraction bound failed')
    records.append(dict(N=N, maximum_endpoint_bound_fraction=endpoint,
      maximum_curvature_bound_fraction=curvature))
  result = dict(created='2026-11-01', input=str(path.relative_to(ROOT)),
    input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    records=records, scope='Fixed physical preparations and bounded density observable at real detuning. '
      'Finite-time bounds only; no assertion of asymptotic mixing-width scaling or optimal detectability.')
  output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
  print(json.dumps(result, indent=2))


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--input', default='effects_2026-11-01.json')
  parser.add_argument('--output', default='bounds_2026-11-01.json')
  main(parser.parse_args())
