"""Microscopic Dicke first-click survival, without a soft-mode approximation.

H = a^dag a + J_z + (a+a^dag) J_x / sqrt(N), J=N/2.
The monitored jump is sqrt(2*kappa)*a. Initial state: ground state of H.
All calculations remain in the permutation-symmetric spin sector.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import eigsh, expm_multiply
from scipy.integrate import simpson

ROOT = Path(__file__).resolve().parents[1]


def matrices(N, cutoff):
  a = diags(np.sqrt(np.arange(1, cutoff)), 1, shape=(cutoff,cutoff), format='csr')
  number = diags(np.arange(cutoff), format='csr')
  k = np.arange(N+1)
  jz = diags(k-N/2, format='csr')
  jp = diags(np.sqrt((N-np.arange(N))*(np.arange(N)+1)), -1,
             shape=(N+1,N+1), format='csr')
  jx = (jp+jp.T)/2
  H = (kron(number, eye(N+1)) + kron(eye(cutoff), jz)
       + kron(a+a.T, jx)/np.sqrt(N)).tocsr()
  return H, np.repeat(np.arange(cutoff),N+1)


def simulate(N, d, cutoff=36, ntime=241):
  started = time.monotonic()
  H, number = matrices(N,cutoff)
  energy, eigvec = eigsh(H,k=2,which='SA',tol=1e-12,
                         v0=np.random.default_rng(817).normal(size=H.shape[0]))
  idx = np.argsort(energy)
  energy, eigvec = energy[idx], eigvec[:,idx]
  psi = eigvec[:,0].astype(complex)
  ell = (4*N)**(1/6)
  kappa = 2*d/ell**4
  generator = -1j*(H-energy[0]*eye(H.shape[0]))-kappa*diags(number)
  tau_end = 7 / max(1.,np.sqrt(d))
  # All observables computed from exact no-jump evolution, no normalization.
  states = expm_multiply(generator*ell**2,psi,start=0,stop=tau_end,num=ntime)
  survival = np.sum(abs(states)**2,axis=1)
  tau = np.linspace(0,tau_end,ntime)
  photon = np.sum(abs(states)**2*number[None,:],axis=1)
  boundary = np.sum(abs(states[:,-2*(N+1):])**2,axis=1)
  return dict(N=N,d=d,cutoff=cutoff,kappa=kappa,ell=ell,
    gap_scaled=float((energy[1]-energy[0])*ell**2),
    number_ground=float(np.sum(abs(psi)**2*number)),
    purity_initial_slope_scaled=float(-4*kappa*np.sum(abs(psi)**2*number)*ell**2),
    maximum_cutoff_weight=float(np.max(boundary)),
    tau=tau.tolist(),survival=survival.tolist(),
    nojump_photon_numerator=photon.tolist(),
    truncated_mean=float(simpson(survival,x=tau)),
    seconds=time.monotonic()-started)


def main():
  p=argparse.ArgumentParser()
  p.add_argument('--N',nargs='+',type=int,default=[16,32,64,128])
  p.add_argument('--d',nargs='+',type=float,default=[.3,1.,3.])
  p.add_argument('--cutoff',type=int,default=36)
  p.add_argument('--output',default='dicke.json')
  args=p.parse_args()
  results=[]
  for d in args.d:
    for N in args.N:
      r=simulate(N,d,args.cutoff)
      results.append(r)
      print({k:v for k,v in r.items() if k not in ('tau','survival','nojump_photon_numerator')},flush=True)
      (ROOT/'data'/args.output).write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
  main()
