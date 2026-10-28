"""Deterministic first-emission statistics of the critical soft mode.

Convention: d rho / d tau = -i[h,rho] + d D[Q]rho,
h=P^2/2+Q^4/4.  The first-click survival is the norm under
h_eff=h-i*d*Q^2/2. No random trajectory sampling is used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad, simpson
from scipy.linalg import eigh_tridiagonal, solve_banded
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]


def ground_grid(n=1201, xmax=9.0, state=0):
  q = np.linspace(-xmax, xmax, n)
  dx = q[1] - q[0]
  diagonal = 1 / dx**2 + q**4 / 4
  off = np.full(n - 1, -0.5 / dx**2)
  energy, vec = eigh_tridiagonal(
    diagonal, off, select="i", select_range=(0, max(1, state))
  )
  psi = vec[:, state].astype(complex)
  return q, diagonal, off, energy, psi


def boundary_shape(u):
  u = np.asarray(u)
  return np.sqrt(2 * np.sqrt(2) * np.pi) / np.sqrt(
    np.sinh(np.sqrt(2) * u) + np.sin(np.sqrt(2) * u)
  )


def analytic_constants(n=2401):
  q, diagonal, off, energy, psi = ground_grid(n)
  dx = q[1] - q[0]
  p0 = float(abs(psi[n // 2])**2 / dx)
  q2 = float(np.sum(abs(psi)**2 * q**2))
  s50 = brentq(lambda s: np.sum(abs(psi)**2 * np.exp(-s*q*q)) - .5, .01, 30)
  cmean = p0 * quad(lambda u: boundary_shape(u), 0, 40, epsabs=1e-11)[0]
  csecond = 2 * p0 * quad(lambda u: u*boundary_shape(u), 0, 40, epsabs=1e-11)[0]
  return dict(E0=float(energy[0]), gap=float(energy[1]-energy[0]),
              Q2=q2, density_at_zero=p0, median_coefficient=float(s50),
              mean_coefficient=cmean, second_moment_coefficient=csecond)


def evolve(d, n=1201, dt_factor=.008, xmax=9., tail=1e-11, state=0):
  """Crank--Nicolson on a centered grid with geometric time steps.

  Each step obeys dt<=dt_factor/max(1,sqrt(d)); a small starting
  step resolves the rapid initial survival decrease for d>>1.
  Norm integration is checked by independent step/grid refinement.
  """
  q, diagonal, off, energy, psi = ground_grid(n, xmax, state)
  diagonal = diagonal - energy[state] - .5j*d*q*q
  dtmax = dt_factor / max(1., np.sqrt(d))
  tmax = 36 / (d * .456) if d < 1 else 40 / np.sqrt(d)
  times = [0.]
  survival = [1.]
  emission = [float(d*np.sum(q*q*abs(psi)**2))]
  t = 0.
  step = 0
  started = time.monotonic()
  # Save every point for quadrature; file output is decimated later.
  while t < tmax and survival[-1] > tail:
    dt = min(dtmax, max(dt_factor/d, .015*t), tmax-t)
    ab = np.empty((3, n), dtype=complex)
    ab[0, 0] = 0
    ab[0, 1:] = .5j*dt*off
    ab[1] = 1 + .5j*dt*diagonal
    ab[2, :-1] = .5j*dt*off
    ab[2, -1] = 0
    rhs = (1-.5j*dt*diagonal)*psi
    rhs[:-1] -= .5j*dt*off*psi[1:]
    rhs[1:] -= .5j*dt*off*psi[:-1]
    psi = solve_banded((1, 1), ab, rhs, overwrite_ab=True,
                       overwrite_b=True, check_finite=False)
    t += dt
    times.append(t)
    survival.append(float(np.vdot(psi, psi).real))
    emission.append(float(d*np.sum(q*q*abs(psi)**2)))
    step += 1
  times = np.asarray(times)
  survival = np.asarray(survival)
  mean = float(simpson(survival, x=times))
  second = float(simpson(2*times*survival, x=times))
  # Integrate the emission density instead of t^(q-1) S for q<1,
  # avoiding the integrable singularity at the time origin.
  density = np.asarray(emission)
  orders = [.25, .5, .75, 1., 1.5, 2., 3.]
  moments = {str(q): float(simpson(times**q*density, x=times))
             for q in orders}
  median = float(np.interp(.5, survival[::-1], times[::-1]))
  idx = np.unique(np.r_[0, np.searchsorted(times, np.geomspace(
    times[1], times[-1], 420)), len(times)-1])
  idx = idx[idx < len(times)]
  return dict(d=float(d), grid=n, dt_factor=dt_factor, xmax=xmax,
              state=state, moments=moments,
              mean=mean, median=median, second_moment=second,
              survival_end=float(survival[-1]), t_end=float(times[-1]),
              maximum_survival_increase=float(np.max(np.diff(survival))),
              emission_probability=float(simpson(density, x=times)),
              mean_flux_relative_difference=float(moments['1.0']/mean-1),
              steps=step, seconds=time.monotonic()-started,
              tau=times[idx].tolist(), survival=survival[idx].tolist())


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--d', type=float, nargs='+', default=[.1, 1, 10, 100, 1000, 10000])
  parser.add_argument('--grid', type=int, default=1201)
  parser.add_argument('--dt-factor', type=float, default=.008)
  parser.add_argument('--output', default='soft_mode.json')
  parser.add_argument('--state', type=int, choices=[0, 1], default=0)
  args = parser.parse_args()
  results = []
  for d in args.d:
    result = evolve(d, n=args.grid, dt_factor=args.dt_factor, state=args.state)
    results.append(result)
    print({k: v for k, v in result.items() if k not in ('tau','survival')}, flush=True)
  path = ROOT / 'data' / args.output
  path.write_text(json.dumps(dict(constants=analytic_constants(), results=results), indent=2)+'\n')


if __name__ == '__main__':
  main()
