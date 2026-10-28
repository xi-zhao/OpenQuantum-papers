"""Build the explicit, verified source bundle; exclude exploratory files."""
from pathlib import Path
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parents[1]
FILES = [
  "README.md", "PRECHECK_B.md",
  "paper/manuscript.tex", "paper/references.bib", "paper/manuscript.pdf",
  "code/spectral_analysis.py", "code/make_figures.py", "code/package_source.py",
  "review/block_parity_review.py", "review/prior_art.md", "review/semiclassical_review.md",
  "data/verified_algebra.json", "data/verified_perron.json",
  "data/verified_precision.json", "data/verified_blocks.json",
  "data/verification_summary.json", "data/scaling_table.csv", "data/splitting_table.csv",
  "figs/scaling.pdf", "figs/splitting.pdf", "figs/scaling.png", "figs/splitting.png",
]


def main():
  missing = [name for name in FILES if not (ROOT / name).is_file()]
  if missing:
    raise FileNotFoundError(f"Deliverable incomplete: {missing}")
  manifest = {
    "format_version": 1,
    "files": {name: {"sha256": hashlib.sha256((ROOT / name).read_bytes()).hexdigest(),
                     "bytes": (ROOT / name).stat().st_size} for name in FILES},
  }
  (ROOT / "data/artifact_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
  destination = ROOT / "pra-ep-source.zip"
  with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for name in FILES + ["data/artifact_manifest.json"]:
      archive.write(ROOT / name, arcname=f"pra-ep-scaling/{name}")
  with zipfile.ZipFile(destination) as archive:
    assert archive.testzip() is None
    for name, record in manifest["files"].items():
      assert hashlib.sha256(archive.read(f"pra-ep-scaling/{name}")).hexdigest() == record["sha256"]
  print(f"Verified bundle: {destination} ({destination.stat().st_size} bytes)")


if __name__ == "__main__":
  main()
