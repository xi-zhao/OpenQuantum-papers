import numpy as np
from math import comb
def bound(n,r):
    q=3*r/4; d=2**n
    ks=range(n+1); N=[comb(n,k)*3**k for k in ks]; p=[(1-q)**(n-k)*(q/3)**k for k in ks]
    best=9
    for tau in sorted(set(p))+[0]:
        B=sum(N[k]*p[k] for k in ks if p[k]>tau); G=sum(N[k]*p[k]**2 for k in ks if p[k]<=tau)
        best=min(best,B+np.sqrt(d/2*G))
    return min(best,1)
for f in ["sim10.dat"]:pass
import subprocess
rows=[l.split() for l in open("sim10.dat")]
for n in [4,6,8]:
    out=subprocess.run(["python3","sim.py",str(n),"40"],capture_output=True,text=True).stdout.split("\n")
    rows+= [l.split() for l in out if l]
ok=True
for n,r,m,s,_ in rows:
    b=bound(int(n),float(r)); ok&= b>=float(m)-3*float(s)
    print(n,r,m,round(b,4))
print("bound holds:",ok)
def lower(n,r):
    q=3*r/4; d=2**n; best=0
    ks=range(n+1); N=[comb(n,k)*3**k for k in ks]; p=[(1-q)**(n-k)*(q/3)**k for k in ks]
    for K in range(n+1):  # good = k<=K
        best=max(best,sum(N[k]*p[k] for k in ks if k<=K)-sum(N[k] for k in ks if k<=K)/d)
    return best
print("lower bound check:",all(lower(int(n),float(r))<=float(m)+3*float(s) for n,r,m,s,_ in rows))
for n,r,m,s,_ in rows[:9]: print(n,r,round(lower(int(n),float(r)),4),m,round(bound(int(n),float(r)),4))
