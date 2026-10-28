"""Analyze audited couplings, physical norms, and signed contractions.

All exponent fits below are finite-window descriptive slopes, not estimates
with statistical sampling errors and not proofs of an N->infinity limit.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigvalsh
from scipy.sparse.linalg import eigs
from biorthogonal_generator import generators, populations
from biorthogonal_audit import mm, bins, encode

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/biorthogonal_audit'
NF=0.15692704039732738


def physical_norm(vec,s,left=False):
  # PIQS rho_j is the multiplicity-summed block. Its physical trace norm
  # equals the sum of block trace norms; the dual observable norm is max
  # over the blocks. W is the transpose dual; transposition preserves norms.
  matrix=np.zeros((s['d'],s['d']),complex)
  matrix[s['row'],s['col']]=vec
  norms=[]
  for j in range(s['N'],-1,-2):
    at=np.flatnonzero(s['twoj']==j)
    block=matrix[np.ix_(at,at)]
    assert np.linalg.norm(block-block.conj().T)<1e-7*max(1,np.linalg.norm(block))
    ev=eigvalsh((block+block.conj().T)/2)
    norms.append(float(max(abs(ev)) if left else sum(abs(ev))))
  return max(norms) if left else sum(norms)


def fitted_slope(N,values):
  fit=np.polyfit(N,np.log(np.abs(values)),1)
  residual=np.log(np.abs(values))-np.polyval(fit,N)
  return dict(alpha=float(-fit[0]),intercept=float(fit[1]),
              max_log_residual=float(max(abs(residual))),
              rms_log_residual=float(np.sqrt(np.mean(residual**2))))


def point_analysis(N):
  tag=f'N{N}_float64_s17_t-0.76'
  d=json.loads((DATA/f'{tag}.json').read_text())
  z=np.load(DATA/f'{tag}.npz'); R,W=z['R'],z['W']
  L0,L1,s=generators(N); L=(L0+d['detuning']*L1).tocsc()
  H=mm(W.T,L@R)
  right_norm=[physical_norm(R[:,i],s) for i in range(2)]
  left_norm=[physical_norm(W[:,i],s,left=True) for i in range(2)]
  Htrace=H*np.array(right_norm)[:,None]/np.array(right_norm)[None,:]
  product=H[0,1]*H[1,0]
  c_ratio=float(np.sqrt(abs(H[1,0]/H[0,1])))
  # Scale c_b/c_a by c_ratio, equalizing the offdiagonal magnitudes.
  Hb=H.copy(); Hb[0,1]*=c_ratio; Hb[1,0]/=c_ratio
  out=dict(N=N,detuning=d['detuning'],g=d['g'],
           Vab=float(H[0,1].real),Vba=float(H[1,0].real),
           trace_norm_gauge=[float(Htrace[0,1].real),float(Htrace[1,0].real)],
           right_trace_norm=right_norm,left_operator_norm=left_norm,
           balanced_gauge=[float(Hb[0,1].real),float(Hb[1,0].real)],
           balanced_scale_ratio=c_ratio,
           invariant_product_error=float(abs(Htrace[0,1]*Htrace[1,0]/product-1)))
  nr=(s['ne'][s['row']]+s['ne'][s['col']])/(2*N)
  group=s['ne'][s['row']]+s['ne'][s['col']]
  out['contractions']={}
  for name,a,b in [('ab',0,1),('ba',1,0)]:
    raw=W[:,a]*(L@R[:,b])
    # Subtract a scalar eigenvalue reference; offdiagonal contraction is unchanged.
    residual=(L@R[:,b])-H[b,b]*R[:,b]
    shifted=W[:,a]*residual
    expected=H[a,b]*R[:,a]
    identity_error=float(np.linalg.norm(residual-expected)/np.linalg.norm(expected))
    out['contractions'][name]=dict(
      raw_signed_sum=encode([sum(raw)])[0],raw_abs_sum=float(sum(abs(raw))),
      shifted_signed_sum=encode([sum(shifted)])[0],
      shifted_abs_sum=float(sum(abs(shifted))),residual_identity_relative_error=identity_error,
      raw_population=float(sum(raw[s['row']==s['col']]).real),
      raw_coherence=float(sum(raw[s['row']!=s['col']]).real),
      raw_absolute_mean_n=float(sum(nr*abs(raw))/sum(abs(raw))),
      shifted_absolute_mean_n=float(sum(nr*abs(shifted))/sum(abs(shifted))),
      raw_low_n_fraction=float(sum(abs(raw[nr<.22]))/sum(abs(raw))),
      shifted_low_n_fraction=float(sum(abs(shifted[nr<.22]))/sum(abs(shifted))),
      excitation_sum_raw=encode(bins(raw,group,2*N+1)),
      excitation_sum_shifted=encode(bins(shifted,group,2*N+1)))
  profiles={}
  for name,vec in [('ra',R[:,0]),('rb',R[:,1]),('la',W[:,0]),('lb',W[:,1])]:
    mask=s['row']==s['col']; k=s['ne'][s['row'][mask]]
    maxdiag=np.zeros(N+1)
    np.maximum.at(maxdiag,k,abs(vec[mask]))
    pop=populations(vec,s)
    profiles[name]=dict(pop_at_fold=float(np.interp(NF,np.arange(N+1)/N,pop.real)),
                        maxdiag_at_fold=float(np.interp(NF,np.arange(N+1)/N,maxdiag)),
                        population_peak_n=float(np.argmax(abs(pop))/N),
                        maxdiag=maxdiag.tolist())
  out['fold_profiles']=profiles
  return out


def nearby_spectrum(N):
  tag=f'N{N}_float64_s17_t-0.76'
  d=json.loads((DATA/f'{tag}.json').read_text()); H,H1=np.array(d['H']),np.array(d['H1'])
  L0,L1,s=generators(N)
  width=2*d['g']/abs(H1[0,0]-H1[1,1])
  rng=np.random.default_rng(32)
  out=[]
  for mult in (-1.,-.5,.5,1.):
    delta=mult*width
    w=eigs(L0+(d['detuning']+delta)*L1,k=4,sigma=float(np.trace(H)/2),
           return_eigenvectors=False,tol=2e-14,v0=rng.normal(size=L0.shape[0]))
    select=np.argsort(abs(w-np.trace(H)/2))[:2]
    actual=np.sort(w[select].real)
    prediction=np.sort(np.linalg.eigvals(H+delta*H1).real)
    out.append(dict(offset=delta,actual=actual.tolist(),prediction=prediction.tolist(),
                    relative_halfgap_error=float((np.diff(prediction)[0]/np.diff(actual)[0])-1),
                    absolute_eigenvalue_error=float(max(abs(actual-prediction)))))
  return dict(N=N,width=width,points=out)


def main():
  records=[point_analysis(N) for N in range(24,57,4)]
  fits={}
  for lo in (24,32,40):
    r=[x for x in records if x['N']>=lo]; sizes=np.array([x['N'] for x in r])
    quantities=dict(Vab=[x['Vab'] for x in r],Vba=[x['Vba'] for x in r],g=[x['g'] for x in r],
      trace_Vab=[x['trace_norm_gauge'][0] for x in r],trace_Vba=[x['trace_norm_gauge'][1] for x in r],
      balanced_Vab=[x['balanced_gauge'][0] for x in r],balanced_Vba=[x['balanced_gauge'][1] for x in r],
      ba_raw_envelope=[x['contractions']['ba']['raw_abs_sum'] for x in r],
      ba_cancellation_fraction=[x['Vba']/x['contractions']['ba']['raw_abs_sum'] for x in r],
      ra_population_fold=[x['fold_profiles']['ra']['pop_at_fold'] for x in r],
      lb_maxdiag_fold=[x['fold_profiles']['lb']['maxdiag_at_fold'] for x in r],
      ra_maxdiag_fold=[x['fold_profiles']['ra']['maxdiag_at_fold'] for x in r])
    fits[f'{lo}-56']={key:fitted_slope(sizes,values) for key,values in quantities.items()}
    fits[f'{lo}-56']['sum_rule_error']=abs((fits[f'{lo}-56']['Vab']['alpha']+
      fits[f'{lo}-56']['Vba']['alpha'])/2-fits[f'{lo}-56']['g']['alpha'])
  local=[dict(N1=a['N'],N2=b['N'],**{k:float(-np.log(b[k]/a[k])/(b['N']-a['N']))
         for k in ['Vab','Vba','g']}) for a,b in zip(records[:-1],records[1:])]
  neighbor=[nearby_spectrum(N) for N in (40,56)]
  summary=dict(records=records,fits=fits,local_exponents=local,nearby_spectrum=neighbor,
               interpretation='Descriptive finite-N fits; no independent WKB coupling prediction.')
  (DATA/'analysis.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
  print('FITS',json.dumps(fits,indent=2))
  print('NEARBY',json.dumps(neighbor,indent=2))
  for r in records:
    x=r['contractions']['ba']
    print('N',r['N'],'trace_norms',r['right_trace_norm'],'left_norms',r['left_operator_norm'],
          'weak diagonal/coherence',x['raw_population'],x['raw_coherence'],
          'residual cancellation',x['shifted_abs_sum']/r['Vba'],
          'mean_n',x['raw_absolute_mean_n'],'low_n_fraction',x['raw_low_n_fraction'])


if __name__=='__main__':
  main()
