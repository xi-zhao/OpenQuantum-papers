"""Signed boundary-flux groups and a deliberately simple partition test.

Both are diagnostics. A hard operator-coordinate cut changes local modes and
is not a physical quantum channel. No failed partition is used as a coupling
prediction. Original complex coefficients are retained before grouping.
"""
from pathlib import Path
import json
import numpy as np
from scipy.sparse.linalg import eigs
from biorthogonal_generator import generators
from biorthogonal_audit import encode
from biorthogonal_flux import accurate_sum

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/coherence_mechanism'
AUDIT=ROOT/'data/biorthogonal_audit'


def groups(t,ids):
  return [dict(label=int(k),signed=encode([accurate_sum(t[ids==k])])[0],
               absolute=float(sum(abs(t[ids==k])))) for k in np.unique(ids)]


def phases(N):
  tag=f'N{N}_float64_s17_t-0.76'
  ref=json.loads((AUDIT/f'{tag}.json').read_text())
  z=np.load(AUDIT/f'{tag}.npz');R,W=z['R'],z['W']
  L0,L1,s=generators(N);L=L0+ref['detuning']*L1;C=L.tocoo()
  nr=s['ne'][s['row']];nc=s['ne'][s['col']]
  mask=((nr+nc)/(2*N)<=.157).astype(int)
  cross=mask[C.row]-mask[C.col];ix=np.flatnonzero(cross)
  i,j=C.row[ix],C.col[ix]
  t=W[i,1]*C.data[ix]*R[j,0]*cross[ix]
  out=dict(N=N,cut=.157,flux=encode([accurate_sum(t)])[0],absolute=float(sum(abs(t))))
  for label,ids in [('signed_q',nr[i]-nc[i]),('twoj',s['twoj'][s['row'][i]]),
                    ('direction',cross[ix]),('deltaq',nr[i]-nc[i]-nr[j]+nc[j])]:
    g=groups(t,ids)
    out[label]=dict(groups=g,grouped_absolute=float(sum(abs(complex(*x['signed'])) for x in g)))
  return out


def partition(N,cut):
  ref=json.loads((AUDIT/f'N{N}_float64_s17_t-0.76.json').read_text())
  L0,L1,s=generators(N);L=L0+ref['detuning']*L1
  coord=(s['ne'][s['row']]+s['ne'][s['col']])/(2*N)
  lo=np.flatnonzero(coord<=cut);hi=np.flatnonzero(coord>cut)
  modes=[]
  rng=np.random.default_rng(124)
  for ids in (hi,lo):
    A=L[ids][:,ids].tocsc()
    w,R=eigs(A,k=4,sigma=-.76,tol=2e-13,
             v0=rng.normal(size=len(ids))+1j*rng.normal(size=len(ids)))
    wl,W=eigs(A.T,k=4,sigma=-.76,tol=2e-13,
              v0=rng.normal(size=len(ids))+1j*rng.normal(size=len(ids)))
    ix=int(np.argmin(abs(w+.76)));iy=int(np.argmin(abs(wl-w[ix])))
    r=R[:,ix];v=W[:,iy];v/=sum(v*r)
    modes.append((w[ix],r,v,w))
  a,b=modes
  ab=sum(a[2]*(L[hi][:,lo]@b[1]));ba=sum(b[2]*(L[lo][:,hi]@a[1]))
  return dict(N=N,cut=cut,detuning=ref['detuning'],
    local_eigenvalues=encode([a[0],b[0]]),candidates=[encode(a[3]),encode(b[3])],
    product=encode([ab*ba])[0],local_sqrt_product=encode([np.sqrt(ab*ba)])[0],
    reference_g=ref['g'],interpretation='Unvalidated hard-cut local modes; not a gap prediction.')


def main():
  DATA.mkdir(exist_ok=True)
  result=dict(phase_groups=[phases(N) for N in range(24,57,4)],
              partition_tests=[partition(N,cut) for N in (24,40) for cut in (.157,.22,.26,.30)])
  (DATA/'phase_partition.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
  for x in result['phase_groups']:
    print(x['N'],'flux',x['flux'][0],'grouped_q',x['signed_q']['grouped_absolute'],
          'grouped_j',x['twoj']['grouped_absolute'],flush=True)
  for x in result['partition_tests']:
    print('partition',x['N'],x['cut'],x['local_eigenvalues'],x['product'],flush=True)


if __name__=='__main__':
  main()
