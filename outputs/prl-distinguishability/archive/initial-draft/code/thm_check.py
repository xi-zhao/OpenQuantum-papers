import numpy as np
from math import comb, log2
rows=[l.split() for l in open("sim10.dat")]
import subprocess
for n in [4,6,8]:
    rows+=[l.split() for l in subprocess.run(["python3","sim.py",str(n),"40"],capture_output=True,text=True).stdout.split("\n") if l]
ok1=ok2=True
for n,r,m,s,_ in rows:
    n=int(n); r=float(r); m=float(m); s=float(s); q=3*r/4; d=2**n
    ks=range(n+1); N=[comb(n,k)*3**k for k in ks]; p=[(1-q)**(n-k)*(q/3)**k for k in ks]
    sk=[-log2(x) for x in p]
    best_up=min(sum(N[k]*p[k] for k in ks if sk[k]<n+a)+2**(-a/2)/np.sqrt(2) for a in np.linspace(0,3*n,600))
    best_lo=max(sum(N[k]*p[k] for k in ks if sk[k]<n-a)-2**(-a) for a in np.linspace(0,n,600))
    part=sum(min(N[k]*p[k],np.sqrt(d*N[k]/2)*p[k]) for k in ks)
    ok1&= best_lo<=m+3*s and m-3*s<=best_up
    ok2&= m-3*s<=part
    print(n,r,round(best_lo,3),m,round(part,3),round(best_up,3))
print("Thm1",ok1,"partition",ok2)
