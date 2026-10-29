from chir import job
from multiprocessing import Pool
if __name__=="__main__":
    args=[(20,'A',1e4,s,'hi',st) for st in (600,1200) for s in (1,-1)]+[(20,'B',31.6227766,s,'hi',st) for st in (600,2400) for s in (1,-1)]
    with Pool(8) as P: res=P.map(job,args)
    for r in res: print(r['loop'],r['T'],r['steps'],r['sgn'],'%.6f'%r['nf'],flush=True)
