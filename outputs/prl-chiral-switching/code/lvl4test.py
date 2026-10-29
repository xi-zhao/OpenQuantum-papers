import numpy as np
from lvl4 import *
g=np.zeros((4,4),complex); g[0,0]=1
for name in ['black','cyan','magenta']:
    O0,Ot,D0,dl,cut=LOOPS[name]; Ds,Os=D0,O0-Ot
    # prepare both branches by adiabatic ramps in Delta from far red / far blue detuning
    def ramp(Da):
        p=lambda t:(Da+(Ds-Da)*min(t/400,1),Os); S=solve_ivp(rhs,[0,500],vec(g),args=(p,),rtol=1e-8,atol=1e-10,method='LSODA'); return S.y[:,-1]
    yA=ramp(-110.); yB=ramp(-20.)
    print(name,'start',(Ds,Os),'n from left ramp %.4f, right ramp %.4f'%(n(yA),n(yB)),flush=True)
    for T in [500.]:
        for sgn,lab in [(-1,'CW'),(1,'CCW')]:
            pass
    for T in [500.]:
        for sgn in (1,-1):
            lab='CW' if sgn>0 else 'CCW'
            for y0,lb in [(yA,'A'),(yB,'B')]:
                S=solve_ivp(rhs,[0,T],y0,args=(mkpath(name,T,sgn),),rtol=1e-8,atol=1e-10,method='LSODA')
                print('  ',lab,'init',lb,'%.4f -> %.4f'%(n(y0),n(S.y[:,-1])),flush=True)
