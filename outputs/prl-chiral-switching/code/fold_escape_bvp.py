"""Finite-time boundary-value continuation of the exact WKB Hamiltonian.

Fix q(0)=q_high and q(T)=q_saddle; solve for p(t). The large-T limit is
checked by endpoint momenta and Hamiltonian energy, not assumed. Seed is a
previous converged switching instanton. The Hamiltonian is indefinite, so
this computes a continued stationary path, not a proven global minimum.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.integrate import solve_ivp,solve_bvp,simpson
from wkbr import F,JAC,H,fixed
from mf import roots
from fold_escape_continuation import basis,sph,mul

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/coherence_mechanism'
OM=2.;V=10.;DF=3.599636890457168


def seed():
  d=3.8;r=roots(d,OM,V);A=fixed(d,OM,V,r[-1]);S=fixed(d,OM,V,r[1])
  angles=[.5033821965491999,2.722432671333336,1.6047795887480882,1.668834806702011]
  section=(r[1]+r[2])/2
  def segment(Y,U,a,sign):
    initial=Y+1e-6*mul(U,sph(a))
    def event(t,y):return y[0]-section
    event.terminal=True
    z=solve_ivp(lambda t,y:sign*F(y,d,OM,V),(0,100),initial,events=event,
                rtol=2e-11,atol=2e-13,dense_output=True,max_step=.5)
    if not len(z.t_events[0]):raise RuntimeError('Seed failed to reach section')
    return z
  f=segment(A,basis(A,d,False),angles[:2],1)
  b=segment(S,basis(S,d,True),angles[2:],-1)
  tf,tb=f.t[-1],b.t[-1]
  t=np.linspace(0,tf+tb,500)
  y=np.empty((6,len(t)))
  before=t<=tf
  y[:,before]=f.sol(t[before]);y[:,~before]=b.sol(tf+tb-t[~before])
  return t,y


def run(scale=1.,manifold=False,tol=2e-7):
  DATA.mkdir(exist_ok=True)
  t0,y0=seed();T=(t0[-1])*scale
  grid=np.linspace(0,1,500)
  # Increase T by padding near the endpoints, retaining the central path.
  core=np.clip(grid*T-(T-t0[-1])/2,0,t0[-1])
  y=np.vstack([np.interp(core,t0,row) for row in y0])
  records=[]
  tag=f'manifold{int(manifold)}_scale{scale:g}_tol{tol:g}'
  schedule=[3.8,3.78,3.76,3.74,3.72,3.70,3.68,3.66,3.64,3.625,3.615,3.608,3.604,3.602,3.601,3.6003,3.59995,3.59980,3.59972]
  for d in schedule:
    start=time.monotonic();r=roots(d,OM,V)
    A=fixed(d,OM,V,r[-1]);S=fixed(d,OM,V,r[1])
    if manifold:
      ua=basis(A,d,False);us=basis(S,d,True)
      qa=np.linalg.qr(ua,mode='complete')[0][:,3:]
      qs=np.linalg.qr(us,mode='complete')[0][:,3:]
      def bc(ya,yb):return np.r_[mul(qa.T,ya-A),mul(qs.T,yb-S)]
    else:
      def bc(ya,yb):return np.r_[ya[:3]-A[:3],yb[:3]-S[:3]]
    solution=solve_bvp(lambda t,z:T*F(z,d,OM,V),bc,grid,y,
                       fun_jac=lambda t,z:T*JAC(z,d,OM,V),
                       tol=tol,max_nodes=15000)
    sample=np.linspace(0,1,3001);path=solution.sol(sample)
    energy=H(path,d,OM,V);velocity=F(path,d,OM,V)[:3]
    action=simpson(np.sum(path[3:]*velocity,axis=0),x=sample*T)
    record=dict(detuning=d,mu=d-DF,T=T,scale=scale,manifold=manifold,requested_tolerance=tol,ok=bool(solution.success),
      status=int(solution.status),message=solution.message,action=float(action),
      energy_max=float(max(abs(energy))),energy_mean=float(np.mean(energy)),
      momentum_end_norms=[float(np.linalg.norm(path[3:,0])),float(np.linalg.norm(path[3:,-1]))],
      collocation_residual_max=float(max(solution.rms_residuals)),
      nodes=len(solution.x),seconds=time.monotonic()-start)
    records.append(record)
    (DATA/f'fold_bvp_{tag}.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
    print(json.dumps(record),flush=True)
    if not solution.success:break
    grid=np.linspace(0,1,600);y=solution.sol(grid)
    np.savez_compressed(DATA/f'fold_bvp_last_{tag}.npz',t=sample*T,y=path,detuning=d)


if __name__=='__main__':
  p=argparse.ArgumentParser(description=__doc__);p.add_argument('--scale',type=float,default=1.)
  p.add_argument('--manifold',action='store_true')
  p.add_argument('--tol',type=float,default=2e-7)
  args=p.parse_args()
  run(args.scale,args.manifold,args.tol)
