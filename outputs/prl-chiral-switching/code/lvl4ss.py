import numpy as np
from lvl import ops, liou, ss_lin
from lvl4 import P4, LS
from scipy.optimize import fsolve
P=P4
def Hlin(Dc,Omw,n3,n4):
    D1=Dc-P['alpha']*n3-P['beta']*n4; D2=-P['gam']*n3-P['eps']*n4
    return np.array([[0,P['Op']/2,0,0],[P['Op']/2,0,P['Oc']/2,0],[0,P['Oc']/2,-D1,Omw/2],[0,0,Omw/2,-D2]],complex)
def G(x,Dc,Omw):
    r=ss_lin(Hlin(Dc,Omw,x[0],x[1]),LS); return np.array([r[2,2].real-x[0],r[3,3].real-x[1]])
def all_roots(Dc,Omw,grid=12):
    sols=[]
    for a in np.linspace(0,0.3,grid):
        for b in np.linspace(0,0.06,grid):
            x,info,ier,msg=fsolve(G,[a,b],args=(Dc,Omw),full_output=True,xtol=1e-12)
            if ier==1 and np.max(abs(G(x,Dc,Omw)))<1e-10 and x.min()>-1e-9 and not any(np.hypot(*(x-s))<1e-6 for s in sols): sols.append(x)
    return sorted(sols,key=lambda s:s[0]+s[1])
def track(path,x0,M=2000):
    """quasi-static branch following: continuation, jump to nearest remaining root at folds"""
    x=np.array(x0); out=[]
    for k in range(M+1):
        Dc,Omw=path(k/M)
        y,info,ier,msg=fsolve(G,x,args=(Dc,Omw),full_output=True,xtol=1e-12)
        if ier!=1 or np.max(abs(G(y,Dc,Omw)))>1e-9 or np.hypot(*(y-x))>0.02:
            rs=all_roots(Dc,Omw,8); 
            # branch lost: pick root reached by relaxing = nearest stable root in n (fold jump)
            y=min(rs,key=lambda s: abs((s[0]+s[1])-(x[0]+x[1])) if np.hypot(*(s-x))>0.005 else 1e9) if rs else x
        x=y; out.append(x[0]+x[1])
    return np.array(out),x
