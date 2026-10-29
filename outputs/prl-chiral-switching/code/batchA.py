import numpy as np, json, sys, os
from multiprocessing import Pool
from chir import job
if __name__=="__main__":
    Ts=list(np.logspace(1.5,4.5,7))
    args=[]
    for N,steps,TT in [(10,600,Ts),(20,600,Ts),(30,600,Ts)]:
        for loop in 'A':
            for T in TT:
                for s in (1,-1): args.append((N,loop,T,s,'hi',steps))
    out='../data/batch.jsonl'
    done=set()
    if os.path.exists(out):
        for l in open(out): r=json.loads(l); done.add((r['N'],r['loop'],round(r['T'],6),r['sgn']))
    args=[a for a in args if (a[0],a[1],round(a[2],6),a[3]) not in done]
    print(len(args),'jobs',flush=True)
    with Pool(9) as P:
        for r in P.imap_unordered(job,args):
            with open(out,'a') as f: f.write(json.dumps(r)+'\n')
            print(r['N'],r['loop'],round(r['T'],1),r['sgn'],round(r['nf'],4),flush=True)
