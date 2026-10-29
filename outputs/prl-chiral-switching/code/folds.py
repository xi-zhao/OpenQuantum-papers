import numpy as np
V=10.
def fold_curve(n_grid=np.linspace(1e-4,0.4999,20000)):
    pts=[]
    for n in n_grid:
        # u^2(1+ n/(1/2-n)) - 2Vn u + 1/4 + n/4/(1/2-n) = 0
        a=1+n/(0.5-n); b=-2*V*n; c=0.25+0.25*n/(0.5-n)
        disc=b*b-4*a*c
        if disc<0: continue
        for u in [(-b+np.sqrt(disc))/(2*a),(-b-np.sqrt(disc))/(2*a)]:
            Om2=n*(u*u+0.25)/(0.25-n/2)
            pts.append((n,u+V*n,np.sqrt(Om2)))
    return np.array(pts)
