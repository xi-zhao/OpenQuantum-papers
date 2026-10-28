"""Boundary-current identity for the weak/strong two-mode couplings.

Chi selects operator entries with (k_row+k_col)/(2N) <= cut. It is a
coordinate mask, not a CPTP projection or a measured probability event.
J_ab = w_a.T [Chi,L] r_b
     = (e_b-e_a) w_a.T Chi r_b + V_ab(D_a-D_b),
where D_i=w_i.T Chi r_i. This holds for an invariant two-mode subspace.
"""
from pathlib import Path
import json
import math
import numpy as np
from biorthogonal_generator import generators
from biorthogonal_audit import mm, encode

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/biorthogonal_audit'


def accurate_sum(a):
  return complex(math.fsum(a.real),math.fsum(a.imag))


def main():
  records=[]
  cuts=[.157,.18,.22,.26,.30,.34,.38,.42]
  for N in range(24,57,4):
    tag=f'N{N}_float64_s17_t-0.76'
    d=json.loads((DATA/f'{tag}.json').read_text())
    z=np.load(DATA/f'{tag}.npz'); R,W=z['R'],z['W']
    L0,L1,s=generators(N); L=L0+d['detuning']*L1; coo=L.tocoo()
    H=mm(W.T,L@R)
    ncoord=(s['ne'][s['row']]+s['ne'][s['col']])/(2*N)
    for cut in cuts:
      mask=(ncoord<=cut).astype(float)
      overlap=mm(W.T,R*mask[:,None])
      edge_mask=mask[coo.row]-mask[coo.col]
      at=np.flatnonzero(edge_mask)
      js,ks=coo.row[at],coo.col[at]
      out=dict(N=N,cut=cut,actual_cut=float(max(ncoord[mask>0])),
               self_overlap=np.diag(overlap).real.tolist())
      for name,a,b in [('ab',0,1),('ba',1,0)]:
        terms=W[js,a]*coo.data[at]*R[ks,b]*edge_mask[at]
        flux=accurate_sum(terms)
        denom=overlap[a,a]-overlap[b,b]
        corrected=flux-(H[b,b]-H[a,a])*overlap[a,b]
        prediction=corrected/denom
        out[name]=dict(flux=encode([flux])[0],absolute_sum=float(sum(abs(terms))),
                       cancellation=float(sum(abs(terms))/abs(flux)),
                       denominator=encode([denom])[0],coupling=encode([prediction])[0],
                       relative_error=float(abs(prediction/H[a,b]-1)))
      records.append(out)
  fit={}
  for cut in cuts:
    r=[x for x in records if x['cut']==cut]
    n=np.array([x['N'] for x in r])
    vals=dict(flux=[abs(complex(*x['ba']['flux'])) for x in r],
              envelope=[x['ba']['absolute_sum'] for x in r],
              cancellation_fraction=[1/x['ba']['cancellation'] for x in r],
              denominator=[abs(complex(*x['ba']['denominator'])) for x in r])
    fit[str(cut)]={key:float(-np.polyfit(n,np.log(val),1)[0]) for key,val in vals.items()}
  (DATA/'boundary_flux.json').write_text(json.dumps(dict(records=records,fits=fit),indent=2)+'\n')
  print(json.dumps(fit,indent=2))
  for r in records:
    if r['N'] in (40,56):
      print(r['N'],r['cut'],'denom',r['ba']['denominator'][0],
            'flux',r['ba']['flux'][0],'cancel',r['ba']['cancellation'],
            'error',r['ba']['relative_error'])


if __name__=='__main__':
  main()
