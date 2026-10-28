"""Mean-field fold checks and physical-state preparation conventions."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from scipy.optimize import root
from pra_modes import fixed, jacobian

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/pra_spectral'


def main():
  n, u, v, delta, omega, interaction, lam = sp.symbols('n u v delta omega V lam', real=True)
  drive = delta-interaction*n
  field = sp.Matrix([-n-omega*v, -drive*v-u/2,
                     drive*u-omega*(1-2*n)/2-v/2])
  jac = field.jacobian([n, u, v])
  char = sp.expand((lam*sp.eye(3)-jac).det())
  stationary = {u:2*n*(delta-interaction*n)/omega, v:-n/omega}
  reduced = sp.factor(char.subs(stationary))
  target = lam*((lam+1)**2+omega**2/2)
  difference = sp.Poly(sp.expand(reduced-target),lam)
  F = n*((delta-interaction*n)**2+sp.Rational(1,4)+omega**2/2)-omega**2/4
  Fn = sp.diff(F,n)
  relations = sp.groebner([F,Fn],delta,interaction,n,omega,order='grevlex')
  residuals = [sp.factor(relations.reduce(c)[1]) for c in difference.all_coeffs()]
  if any(r != 0 for r in residuals):
    raise RuntimeError(f'Fold polynomial identity failed: {residuals}')
  def equations(x):
    density, detuning = x; shift = detuning-10*density
    return [density*(shift**2+.25+2)-1,
            shift**2+.25+2-20*density*shift]
  folds=[]
  for guess in ([.157,3.60],[.414,4.70]):
    solution=root(equations,guess,tol=1e-12)
    density,detuning=solution.x
    error=float(np.linalg.norm(equations(solution.x)))
    if error>1e-10:raise RuntimeError('Fold root failed')
    state=fixed(detuning,density)
    values=np.linalg.eigvals(jacobian(detuning,state))
    folds.append(dict(population=float(density),detuning=float(detuning),
      state=state.tolist(),eigenvalues=[[float(z.real),float(z.imag)] for z in values],
      residual=error))
  out=dict(fold_identity=str(reduced),difference_remainders=[str(r) for r in residuals],
           folds=folds,parameters=dict(omega=2.,interaction=10.,gamma=1.))
  DATA.mkdir(exist_ok=True)
  (DATA/'meanfield.json').write_text(json.dumps(out,indent=2)+'\n')
  print(json.dumps(out,indent=2))


if __name__=='__main__':main()
