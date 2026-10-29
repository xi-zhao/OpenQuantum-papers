import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, binom
from theory2 import lnT
from exponent import E, qc
L=np.log2
rc=4*qc/3
def Hq(q): return -(1-q)*L(1-q)-q*L(q/3)
def sig(q):
    s0,s1=-L(1-q),-L(q/3); m=(1-q)*s0+q*s1; return np.sqrt((1-q)*s0**2+q*s1**2-m*m)
# exact sims
sim={}
for f,ns in [("sim10.dat",[10])]:
    for l in open(f): n,r,m,s,_=l.split(); sim.setdefault(int(n),[]).append((float(r),float(m),float(s)))
import subprocess
for n in [6,8]:
    for l in subprocess.run(["python3","sim.py",str(n),"40"],capture_output=True,text=True).stdout.split("\n"):
        if l: a=l.split(); sim.setdefault(n,[]).append((float(a[1]),float(a[2]),float(a[3])))
fig,ax=plt.subplots(1,2,figsize=(7,2.7))
rs=np.linspace(0.02,0.7,60); cols=plt.cm.viridis(np.linspace(0,.85,6))
for c,n in zip(cols,[6,8,10,40,160,640]):
    ax[0].plot(rs,[np.exp(lnT(n,r,6001)) for r in rs],color=c,lw=1.2,label=f"$n={n}$")
    if n in sim:
        a=np.array(sim[n]); ax[0].errorbar(a[:,0],a[:,1],a[:,2],fmt='o',ms=3,color=c,mfc='none')
ax[0].axvline(rc,ls='--',c='gray',lw=.8); ax[0].axvline(0.4227,ls=':',c='gray',lw=.8)
ax[0].set_xlabel("depolarizing strength $r$"); ax[0].set_ylabel(r"$\overline{T}$"); ax[0].legend(fontsize=6.5,frameon=False)
ax[0].text(0.03,0.04,"(a)",transform=ax[0].transAxes)
# collapse
for c,n in zip(cols[2:],[10,40,160,640]):
    rr=np.linspace(0.12,0.42,40); q=3*rr/4
    x=np.sqrt(n)*(1-Hq(q))/sig(q)
    ax[1].plot(x,[np.exp(lnT(n,r,6001)) for r in rr],color=c,lw=1.2,label=f"$n={n}$")
xx=np.linspace(-4,4,200); ax[1].plot(xx,norm.cdf(xx),'k--',lw=1,label=r"$\Phi(x)$")
ax[1].set_xlim(-4,4); ax[1].set_xlabel(r"$x=\sqrt{n}\,[1-H(\mathbf{p})]/\sigma$"); ax[1].legend(fontsize=6.5,frameon=False)
ax[1].text(0.03,0.9,"(b)",transform=ax[1].transAxes)
plt.tight_layout(); plt.savefig("../fig1.pdf")
# exponent figure
fig,ax=plt.subplots(figsize=(3.4,2.5))
r=np.linspace(0,1,400); q=3*r/4
ax.plot(r,[E(x) for x in r],'k',lw=1.4,label=r"$E(r)$ (this work)")
H2=-L((1-q)**2+q*q/3); ax.plot(r,np.maximum(0,(H2-1)/2),'C1--',lw=1.1,label=r"$\frac{1}{2}[H_2-1]_+$")
# lower-bound exponent from projector bound
from scipy.optimize import brentq
D=lambda a,b: a*L(a/b)+(1-a)*L((1-a)/(1-b))
h=lambda f: -f*L(f)-(1-f)*L(1-f)+f*L(3)
def Elb(rv):
    qq=3*rv/4
    if qq<=qc: return 0
    t=brentq(lambda t: D(t,qq)+h(t)-1,1e-9,qc); return D(t,qq)
m=(r>rc+1e-3)&(r<0.655); ax.plot(r[m],[Elb(x) for x in r[m]],'C0:',lw=1.2,label="lower-bound exponent")
for n,mk in [(200,'s'),(800,'o')]:
    rp=np.array([.35,.4,.5,.6,.7,.8,.9]); ax.plot(rp,[-lnT(n,x,30001)/np.log(2)/n for x in rp],mk,ms=3.5,mfc='none',c='C2',label=f"det. eq. $n={n}$")
ax.axvline(rc,ls='--',c='gray',lw=.7); ax.axvline(0.6075,ls='-.',c='gray',lw=.7)
ax.set_xlabel("$r$"); ax.set_ylabel(r"$-\frac{1}{n}\log_2\overline{T}$"); ax.set_ylim(-0.02,0.9); ax.legend(fontsize=6.5,frameon=False,loc='upper left')
plt.tight_layout(); plt.savefig("../fig2.pdf")
print("ok")
