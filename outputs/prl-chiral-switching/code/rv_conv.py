import numpy as np, sys
from chir import job
from multiprocessing import Pool
if __name__=="__main__":
    args=[(20,'A',1000.0,s,'hi',st) for st in (300,600,1200) for s in (1,-1)]
    with Pool(6) as P: res=P.map(job,args)
    for r in res: print(r['steps'],r['sgn'],'%.6f'%r['nf'],flush=True)
