"""Analytic spectral-parameter derivatives of the coherence continued fraction.

Sigma'= -U G (I-Sigma') G D
Sigma''= U [2G(I-Sigma')G(I-Sigma')G + G Sigma'' G]D.
The repeated Sigma symbols on the right belong to the next coherence level.
"""
from pathlib import Path
import json
import argparse
import numpy as np
from scipy.linalg import solve,eig
from coherence_fraction import Fraction,DATA,nearest_real_index
from biorthogonal_audit import mm


def derivatives(model,z,depth=None):
  if not np.isreal(z):
    raise ValueError('Conjugate-sector derivatives require real z')
  top=model.qmax if depth is None else depth
  if not isinstance(top,(int,np.integer)) or not 0<=top<=model.qmax:
    raise ValueError('depth must be an available nonnegative coherence order')
  S=np.zeros_like(model.A[top]);S1=S.copy();S2=S.copy()
  for k in range(top,0,-1):
    X=z*np.eye(len(model.A[k]))-model.A[k]-S
    D=model.down[k-1];U=model.up[k-1]
    Gd=solve(X,D,assume_a='gen')
    E=np.eye(len(X))-S1
    G1=solve(X,mm(E,Gd),assume_a='gen')
    G2=solve(X,2*mm(E,G1)+mm(S2,Gd),assume_a='gen')
    S=mm(U,Gd);S1=-mm(U,G1);S2=mm(U,G2)
  return model.A[0].real+2*S.real,2*S1.real,2*S2.real


def branch_derivatives(model,z):
  K,K1,K2=derivatives(model,z)
  w,vl,vr=eig(K,left=True,right=True)
  index=nearest_real_index(w,z);mu=w[index]
  r=vr[:,index].real;l=vl[:,index].real
  l/=np.dot(l,r)
  mu1=float(np.dot(l,mm(K1,r[:,None])[:,0]))
  n=len(r)
  augmented=np.zeros((n+1,n+1))
  augmented[:n,:n]=K-mu.real*np.eye(n)
  augmented[:n,n]=r;augmented[n,:n]=l
  rhs=np.r_[-mm(K1-mu1*np.eye(n),r[:,None])[:,0],0.]
  r1=solve(augmented,rhs)[:n]
  mu2=float(np.dot(l,mm(K2,r[:,None])[:,0]+2*mm(K1,r1[:,None])[:,0]))
  f=float(mu.real-z)
  out=dict(z=z,mu=float(mu.real),f=f,mu_prime=mu1,f_prime=mu1-1,
           f_second=mu2,quadratic_g=float(np.sqrt(-2*f/mu2)) if f<0 and mu2>0 else None,
           mode_residual=float(np.linalg.norm(mm(K,r[:,None])[:,0]-mu.real*r)),
           differentiated_residual=float(np.linalg.norm(mm(K-mu.real*np.eye(n),r1[:,None])[:,0]+
                       mm(K1-mu1*np.eye(n),r[:,None])[:,0])))
  return out


def main(N):
  source=json.loads((DATA/f'fraction_N{N}_q{N}.json').read_text())
  model=Fraction(N,source['detuning'])
  z=source['center'];out=branch_derivatives(model,z)
  dz=-out['f_prime']/out['f_second']
  if abs(dz)>1e-12:
    out=branch_derivatives(model,z+dz)
  out.update(N=N,reference_g=source['reference_g'],
              relative_prediction_error=out['quadratic_g']/source['reference_g']-1)
  (DATA/f'fraction_derivatives_N{N}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
  print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':
  p=argparse.ArgumentParser(description=__doc__);p.add_argument('N',type=int)
  main(p.parse_args().N)
