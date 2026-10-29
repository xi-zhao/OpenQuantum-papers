import numpy as np, json, sys
from multiprocessing import Pool
from gmamv import gmam, fixed, roots, Om, V
nf=json.load(open('../data/normalform.json'))
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
def val(D):   # hi -> saddle inside wedge (validation vs exact s_hl)
    r=roots(D,Om,V); A=fixed(D,Om,V,r[2])[:3]; Z=fixed(D,Om,V,r[1])[:3]
    P=np.array([fixed(D,Om,V,x)[:3] for x in np.linspace(r[2],r[1],40)]); out=[]
    for m in (40,80,160): s,P=gmam(D,A,Z,P,m); out.append(s)
    return ['val',D,out]
def ghost(N):  # pre-registered: hi -> fold MF state q_f at D = Re Delta_EP(N)
    D=EP[N][1]; r=roots(D,Om,V); A=fixed(D,Om,V,r[-1])[:3]; Z=fixed(nf['Df'],Om,V,nf['nf'])[:3]
    P=np.array([A+(Z-A)*s for s in np.linspace(0,1,40)]); out=[]
    for m in (40,80,160): s,P=gmam(D,A,Z,P,m); out.append(s)
    return ['ghost',N,out]
if __name__=="__main__":
    jobs=[(val,D) for D in (3.62,3.7,3.8,4.0,4.2)]+[(ghost,N) for N in (24,32,40,48,56)]
    with Pool(9) as P:
        for r in P.imap_unordered(lambda a: a[0](a[1]),jobs) if False else P.starmap(lambda f,x: f(x),jobs) if False else []: pass
