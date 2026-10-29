import numpy as np, time
from dyn import ems
from mf import roots
for N in [10,20,30]:
    t=time.time(); g,rss,(hi,lo)=ems(N,4.1,2.0,10.)
    from ex import obs_vec
    tr,nv=obs_vec(N)
    print(N,g,(nv@hi).real/N,(nv@lo).real/N,(nv@rss).real/N, (tr@hi).real, time.time()-t)
print(roots(4.1,2.0,10.))
