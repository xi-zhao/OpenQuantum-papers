"""Excited-mode convergence under coherence-order truncation.

The Liouvillian couples q=k_row-k_col only to q and q+-1. Restricting
|q|<=qmax is an algebraic diagnostic, not a completely positive model.
No truncated model is claimed to predict the gap unless the relevant real
pair is present and converged. All candidates are saved.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.sparse.linalg import eigs
from biorthogonal_generator import generators

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/coherence_mechanism'


def run(N):
  start=time.monotonic()
  ref=json.loads((ROOT/f'data/biorthogonal_audit/N{N}_float64_s17_t-0.76.json').read_text())
  L0,L1,s=generators(N)
  L=(L0+ref['detuning']*L1).tocsc()
  q=s['ne'][s['row']]-s['ne'][s['col']]
  C=L.tocoo()
  assert max(abs(q[C.row]-q[C.col]))==1
  records=[]
  for qmax in sorted(set([0,1,2,3,4,6,8,10,12,14,16,20,24,N])):
    if qmax>N: continue
    ids=np.flatnonzero(abs(q)<=qmax)
    A=L[ids][:,ids].tocsc()
    rng=np.random.default_rng(247+qmax)
    w,R=eigs(A,k=6,sigma=-.76,tol=2e-13,
             v0=rng.normal(size=len(ids))+1j*rng.normal(size=len(ids)))
    order=np.argsort(abs(w+.76)); w=w[order];R=R[:,order]
    real=[i for i,z in enumerate(w) if abs(z.imag)<1e-7]
    pair=np.sort(w[real[:2]].real) if len(real)>=2 else None
    records.append(dict(qmax=qmax,dimension=len(ids),
      candidates=[[float(z.real),float(z.imag)] for z in w],
      pair=None if pair is None else pair.tolist(),
      halfgap=None if pair is None else float(np.diff(pair)[0]/2),
      relative_gap_error=None if pair is None else float(np.diff(pair)[0]/(2*ref['g'])-1),
      residual_max=float(max(np.linalg.norm(A@R[:,i]-w[i]*R[:,i])/np.linalg.norm(R[:,i]) for i in range(6)))))
    print(N,qmax,'pair',pair,'gap_error',records[-1]['relative_gap_error'],flush=True)
  DATA.mkdir(exist_ok=True)
  out=dict(N=N,detuning=ref['detuning'],reference_g=ref['g'],
           reference_pair=ref['lam'],records=records,seconds=time.monotonic()-start)
  (DATA/f'coherence_N{N}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
  parser=argparse.ArgumentParser(description=__doc__)
  parser.add_argument('N',type=int)
  run(parser.parse_args().N)
