"""Package the PRA manuscript and its explicitly selected numerical evidence.

No downloads, installation, or research calculations are performed. The ZIP is
built from an allowlist; exploratory manuscripts and unrelated outputs stay out.
"""
from pathlib import Path
import ast
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parents[1]
CODE = (
  'biorthogonal_generator', 'biorthogonal_audit', 'biorthogonal_validate',
  'biorthogonal_analysis', 'biorthogonal_flux', 'coherence_fraction',
  'coherence_fraction_derivatives', 'coherence_validate', 'coherence_recurrence',
  'mf', 'pra_modes', 'pra_meanfield', 'pra_ep_verify', 'pra_figures', 'pra_validate',
  'coherence_response', 'coherence_response_validate', 'pra_package',
)
FIGURES = ('spectrum', 'residues', 'couplings', 'coherence', 'dynamics', 'ep_loop')


def payload():
  paths = [ROOT/'PRA_README.md', ROOT/'PRA_REPORT_2026-10-28.md',
    ROOT/'PRA_REVISION_2026-10-29.md', ROOT/'PRA_RESPONSE_2026-10-30.md', ROOT/'requirements-pra.txt',
    ROOT/'paper/references_suggested_2026-10-29.bib']
  paths += [ROOT/f'code/{name}.py' for name in CODE]
  paths += [ROOT/f'paper/{name}.{ext}' for name in ('main', 'supplement') for ext in ('tex', 'pdf')]
  paths += [ROOT/f'figs/pra_{name}.{ext}' for name in FIGURES for ext in ('pdf', 'png')]
  audit = ROOT/'data/biorthogonal_audit'
  paths += [audit/name for name in ('inputs.json', 'analysis.json', 'boundary_flux.json',
    'generator_validation.json', 'N56_piqs_s17_t-0.76.json')]
  paths += [audit/f'N{N}_float64_s17_t-0.76.{ext}' for N in range(24, 57, 4) for ext in ('json', 'npz')]
  coherence = ROOT/'data/coherence_mechanism'
  paths += [coherence/f'fraction_N{N}_q{N}.json' for N in (24, 32, 40, 48)]
  paths += [coherence/f'fraction_derivatives_N{N}.json' for N in (24, 32, 40, 48)]
  paths += [coherence/name for name in ('validation.json', 'failure_path_validation.json',
    'derivative_validation.json', 'coherence_N40.json', 'coherence_N56.json')]
  pra = ROOT/'data/pra_spectral'
  paths += [pra/name for name in ('meanfield.json', 'modes_N32.json', 'modes_N40_quick.json',
    'mode_stability.json', 'validation_repeat.json')]
  paths += sorted(pra.glob('ep_*.json'))+sorted(pra.glob('ep_*_vectors.npz'))
  paths += [ROOT/'data/coherence_response'/name for name in
    ('response_2026-10-30.json', 'validation_2026-10-30.json')]
  return sorted(set(paths))


def main():
  paths = payload()
  if any(not p.is_file() for p in paths):
    raise FileNotFoundError([str(p) for p in paths if not p.is_file()])
  for p in paths:
    if p.suffix == '.py':
      ast.parse(p.read_text(), filename=str(p))
    if p.suffix == '.json':
      json.loads(p.read_text())
  total = sum(p.stat().st_size for p in paths)
  if total > 100*1024**2:
    raise RuntimeError('Unexpected package growth; review output allocation before packaging')
  files = {str(p.relative_to(ROOT)): dict(bytes=p.stat().st_size,
    sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths}
  manifest = dict(title='Finite-size decay-mode hybridization near a mean-field fold in an open Rydberg ensemble',
    created='2026-10-30', files=files, total_bytes=total,
    evidence={
      'real gaps and gauge dependence': 'data/biorthogonal_audit/analysis.json',
      'signed boundary flux': 'data/biorthogonal_audit/boundary_flux.json',
      'physical residues and transients': 'data/pra_spectral/modes_N32.json',
      'independent size control': 'data/pra_spectral/modes_N40_quick.json',
      'small-N physical evolution and repeatability': 'data/pra_spectral/validation_repeat.json',
      'complex EPs and monodromy': 'data/pra_spectral/ep_summary.json',
      'Schur back-substitution and failure paths': 'data/coherence_mechanism/validation.json',
      'literature attribution and manuscript build checks': 'PRA_REVISION_2026-10-29.md',
      'initial-coherence source and observable residues': 'data/coherence_response/response_2026-10-30.json',
      'dephasing preparation and independent integrator': 'data/coherence_response/validation_2026-10-30.json'},
    interpretation='Finite-size numerical evidence; no established sigma=Phi identity or thermodynamic-limit mechanism.',
    review_status='Local mathematical/data/layout checks and independent review of the existing finite-size mathematical claims completed; scope and response review are recorded in PRA_RESPONSE_2026-10-30.md.')
  manifest_path = ROOT/'data/pra_spectral/artifact_manifest.json'
  manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
  outdir = ROOT/'deliverables'
  outdir.mkdir(exist_ok=True)
  archive = outdir/'pra-decay-mode-source-2026-10-30.zip'
  with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
    for p in paths+[manifest_path]:
      bundle.write(p, p.relative_to(ROOT))
  with zipfile.ZipFile(archive) as bundle:
    assert bundle.testzip() is None
    for name, record in files.items():
      assert hashlib.sha256(bundle.read(name)).hexdigest() == record['sha256']
  print(json.dumps(dict(archive=str(archive), files=len(files)+1,
    uncompressed_bytes=total, archive_bytes=archive.stat().st_size,
    sha256=hashlib.sha256(archive.read_bytes()).hexdigest())))


if __name__ == '__main__':
  main()
