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
  'physical_effects', 'physical_effects_validate', 'physical_effects_sensitivity',
  'physical_effects_bound', 'physical_effects_review', 'physical_effects_figures',
)
RELEASE = '2026-11-02-r2'
FIGURES = ('spectrum', 'residues', 'couplings', 'coherence', 'dynamics', 'ep_loop')


def payload():
  paths = [ROOT/'PRA_README.md', ROOT/'PRA_REPORT_2026-10-28.md',
    ROOT/'PRA_REVISION_2026-10-29.md', ROOT/'PRA_RESPONSE_2026-10-30.md',
    ROOT/'PRA_REVISION_2026-11-02.md', ROOT/'PRA_CONCISE_2026-11-02.md',
    ROOT/'PHYSICAL_EFFECT_TEST_2026-11-01.md',
    ROOT/'PHYSICAL_EFFECT_REPORT_2026-11-01.md', ROOT/'requirements-pra.txt',
    ROOT/'paper/references_suggested_2026-10-29.bib',
    ROOT/'data/pra_spectral/artifact_manifest.json',
    ROOT/'data/pra_spectral/artifact_manifest_2026-11-02.json']
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
  paths += [ROOT/f'figs/physical_effects_2026-11-01_v2.{ext}' for ext in ('pdf', 'png', 'csv')]
  paths += [ROOT/'data/physical_effects'/name for name in (
    'artifact_manifest_2026-11-01.json', 'effects_2026-11-01.json',
    'pilot_2026-11-01.json', 'validation_2026-11-01.json',
    'validation_failure_2026-11-01.json', 'sensitivity_2026-11-01.json',
    'bounds_2026-11-01.json', 'review_2026-11-01.json')]
  return sorted(set(paths))


def main():
  manifest_path = ROOT/f'data/pra_spectral/artifact_manifest_{RELEASE}.json'
  archive = ROOT/f'deliverables/pra-decay-mode-source-{RELEASE}.zip'
  if manifest_path.exists() or archive.exists():
    raise FileExistsError('Preserve the existing release; choose a new RELEASE before packaging')
  frozen = json.loads((ROOT/'data/physical_effects/artifact_manifest_2026-11-01.json').read_text())
  for name, record in frozen['files'].items():
    if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != record['sha256']:
      raise RuntimeError(f'Frozen physical-effect evidence changed: {name}')
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
    created=RELEASE[:10], release=RELEASE, files=files, total_bytes=total,
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
      'dephasing preparation and independent integrator': 'data/coherence_response/validation_2026-10-30.json',
      'fixed-preparation signals and independent-atom controls': 'data/physical_effects/effects_2026-11-01.json',
      'full-generator detuning derivatives and two-pole bound': 'data/physical_effects/sensitivity_2026-11-01.json',
      'full-spin, ODE, and phase-quadrature checks': 'data/physical_effects/validation_2026-11-01.json',
      'retained phase-quadrature failure': 'data/physical_effects/validation_failure_2026-11-01.json',
      'full-signal trace-contraction bounds': 'data/physical_effects/bounds_2026-11-01.json',
      'spectral identity and independent parameter differences': 'data/physical_effects/review_2026-11-01.json',
      'PRA interpretation and manuscript checks': 'PRA_REVISION_2026-11-02.md',
      'main-text concision and evidence placement': 'PRA_CONCISE_2026-11-02.md'},
    interpretation='Finite-size spectral hybridization with no strong narrow density response in the tested '
      'fixed preparations and time window; no established sigma=Phi identity or thermodynamic-limit mechanism.',
    historical_manifests={'data/pra_spectral/artifact_manifest.json':
      'Frozen 2026-10-30 snapshot; its manuscript hashes describe the earlier archive, not this release.',
      'data/pra_spectral/artifact_manifest_2026-11-02.json':
      'Frozen first 2026-11-02 revision; its manuscript hashes precede the main-text concision.'},
    review_status='Calculation reviews and manuscript checks have distinct scopes; see '
      'PHYSICAL_EFFECT_REPORT_2026-11-01.md, PRA_REVISION_2026-11-02.md, '
      'and PRA_CONCISE_2026-11-02.md.')
  manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
  archive.parent.mkdir(exist_ok=True)
  with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
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
