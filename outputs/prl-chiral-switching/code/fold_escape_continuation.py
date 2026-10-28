"""Continue a previously converged high-to-saddle instanton toward the fold.

This tests a new exploratory candidate sigma_infinity = S_high->fold/2.
It uses no Liouvillian gaps or eigenvectors. Detuning-dependent saddle roots
and the classical reduced WKB Hamiltonian are the only dynamical inputs.
The initial angles come from the earlier D=3.8 switching-action calculation.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.linalg import orthogonal_procrustes
from wkbr import F,JAC,H,fixed
from mf import roots

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/coherence_mechanism'
OM=2.0; V=10.0; DF=3.599636890457168


def mul(a,b):
  return np.einsum('ij,j->i',a,b,optimize=False)


def basis(y,d,stable,previous=None):
  w,r=np.linalg.eig(JAC(y,d,OM,V)); columns=[]
  for i,z in enumerate(w):
    if (z.real<0)!=stable: continue
    if abs(z.imag)<1e-9: columns.append(r[:,i].real)
    elif z.imag>0: columns.extend([r[:,i].real,r[:,i].imag])
  U=np.linalg.qr(np.column_stack(columns))[0]
  assert U.shape==(6,3)
  if previous is not None:
    Q,_=orthogonal_procrustes(U,previous)
    U=np.einsum('ij,jk->ik',U,Q,optimize=False)
  return U


def sph(a):
  return np.array([np.cos(a[0]),np.sin(a[0])*np.cos(a[1]),np.sin(a[0])*np.sin(a[1])])


def ray(a,Y,U,d,section,sign,eps=1e-6):
  ds=np.array([[-np.sin(a[0]),0.],
    [np.cos(a[0])*np.cos(a[1]),-np.sin(a[0])*np.sin(a[1])],
    [np.cos(a[0])*np.sin(a[1]),np.sin(a[0])*np.cos(a[1])]])
  sensitivity=eps*np.einsum('ij,jk->ik',U,ds,optimize=False)
  initial=np.r_[Y+eps*mul(U,sph(a)),0.,sensitivity.ravel()]
  def fun(t,z):
    f=F(z[:6],d,OM,V)
    jac=JAC(z[:6],d,OM,V)
    tangent=np.einsum('ij,jk->ik',jac,z[7:].reshape(6,2),optimize=False)
    return sign*np.r_[f,np.dot(z[3:6],f[:3]),tangent.ravel()]
  def hit(t,z): return z[0]-section
  hit.terminal=True
  def blow(t,z): return 15-max(abs(z[:6]))
  blow.terminal=True
  solution=solve_ivp(fun,(0.,500.),initial,rtol=2e-11,atol=2e-13,
                    events=(hit,blow),max_step=1.)
  if len(solution.t_events[0]):
    end=solution.y_events[0][0]
    tangent=end[7:].reshape(6,2)
    velocity=sign*F(end[:6],d,OM,V)
    tangent-=np.outer(velocity,tangent[0]/velocity[0])
    return end[:7],float(solution.t_events[0][0]),tangent
  return None


def run(eps=1e-6):
  DATA.mkdir(exist_ok=True)
  schedule=[3.8,3.78,3.76,3.74,3.72,3.70,3.68,3.66,3.64,3.625,3.615,3.608,3.604,3.602,3.601,3.6003,3.59995,3.59980,3.59972]
  x=np.array([0.5033821051895226,2.722432367131719,1.6047795886584613,1.6688348067092307])
  UA=US=None; records=[]
  for d in schedule:
    start=time.monotonic()
    r=roots(d,OM,V)
    if len(r)!=3: raise ValueError('Not in bistable region')
    A=fixed(d,OM,V,r[2]); S=fixed(d,OM,V,r[1])
    UA=basis(A,d,False,UA); US=basis(S,d,True,US)
    section=(r[1]+r[2])/2
    cache={}
    def evaluate(y):
      key=tuple(y)
      if key not in cache:
        f=ray(y[:2],A,UA,d,section,1,eps); b=ray(y[2:],S,US,d,section,-1,eps)
        cache.clear()
        cache[key]=(f,b)
      return cache[key]
    def mismatch(y):
      f,b=evaluate(y)
      return np.full(5,10.) if f is None or b is None else f[0][1:6]-b[0][1:6]
    def jacobian(y):
      f,b=evaluate(y)
      return np.zeros((5,4)) if f is None or b is None else np.column_stack([f[2][1:6],-b[2][1:6]])
    fit=least_squares(mismatch,x,jac=jacobian,xtol=1e-13,ftol=1e-13,gtol=1e-13,
                      max_nfev=140)
    x=fit.x
    f=ray(x[:2],A,UA,d,section,1,eps); b=ray(x[2:],S,US,d,section,-1,eps)
    miss=float(np.linalg.norm(mismatch(x)))
    record=dict(detuning=d,distance_from_fold=d-DF,eps=eps,miss=miss,
      ok=bool(f is not None and b is not None and miss<1e-6),
      action=None if f is None or b is None else float(f[0][6]-b[0][6]),
      section=float(section),angles=x.tolist(),
      energy_error=None if f is None or b is None else float(max(abs(H(f[0][:6],d,OM,V)),abs(H(b[0][:6],d,OM,V)))),
      travel_times=None if f is None or b is None else [f[1],b[1]],seconds=time.monotonic()-start)
    records.append(record)
    (DATA/f'fold_escape_eps{eps:g}.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
    print(json.dumps(record),flush=True)
    if not record['ok']:
      break


if __name__=='__main__':
  parser=argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--eps',type=float,default=1e-6)
  run(parser.parse_args().eps)
