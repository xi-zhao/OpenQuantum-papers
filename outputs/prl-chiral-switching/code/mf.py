import numpy as np
# Mean-field Rydberg 2-level: H_i = -Delta n + Om/2 sx + V <n> n ; L = sqrt(g) s-, g=1
# steady state: n = (Om^2/4)/((Delta - V n)^2 + 1/4 + Om^2/2)
def roots(Delta,Om,V):
    # n[(D - V n)^2 + 1/4 + Om^2/2] - Om^2/4 = 0 -> V^2 n^3 - 2 D V n^2 + (D^2+1/4+Om^2/2) n - Om^2/4
    r=np.roots([V**2,-2*Delta*V,Delta**2+0.25+Om**2/2,-Om**2/4])
    return np.sort(r[np.abs(r.imag)<1e-9].real)
def cusp(V):
    # cubic a n^3+b n^2+c n+d with triple root n0: b=-3a n0, c=3a n0^2, d=-a n0^3
    # -2DV=-3V^2 n0 -> D=3V n0/2 ; D^2+1/4+Om^2/2=3V^2 n0^2 ; Om^2/4=V^2 n0^3
    # => 9V^2n0^2/4+1/4+2V^2n0^3=3V^2n0^2 -> 2V^2 n0^3 - 3/4 V^2 n0^2 + 1/4=0
    r=np.roots([2*V**2,-0.75*V**2,0,0.25]); r=r[(abs(r.imag)<1e-9)&(r.real>0)].real
    out=[]
    for n0 in r:
        out.append((n0,1.5*V*n0,2*V*n0**1.5))
    return out
if __name__=="__main__":
    for V in [5,8,10,15,20]:
        print(V,cusp(V))
