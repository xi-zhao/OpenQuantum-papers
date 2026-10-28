"""Package this research note without downloaded literature or build products."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parents[1]

def deliverables():
  paths = [ROOT/'.gitignore', ROOT/'README.md', ROOT/'sources/README.md']
  for folder in ['code', 'data', 'figures', 'paper', 'review']:
    paths += [p for p in (ROOT/folder).iterdir() if p.is_file()
              and p.suffix in {'.py','.json','.tsv','.pdf','.png','.tex','.bib','.md'}
              and p.name != 'artifact_manifest.json']
  return sorted(paths)


def main():
  paths = deliverables()
  manifest = {str(p.relative_to(ROOT)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
              'bytes':p.stat().st_size} for p in paths}
  mpath = ROOT/'data/artifact_manifest.json'
  mpath.write_text(json.dumps({'status':'Research note; PRL novelty criterion not met.',
                              'files':manifest},indent=2)+'\n')
  archive = ROOT/'critical-noise-research.zip'
  with ZipFile(archive,'w',ZIP_DEFLATED) as z:
    for p in paths+[mpath]:
      z.write(p,Path('prl-critical-noise')/p.relative_to(ROOT))
  with ZipFile(archive) as z:
    assert z.testzip() is None
    for rel, rec in manifest.items():
      assert hashlib.sha256(z.read('prl-critical-noise/'+rel)).hexdigest()==rec['sha256']
  print(f'{archive}: {archive.stat().st_size} bytes, {len(paths)+1} files; hashes verified')

if __name__=='__main__':
  main()
