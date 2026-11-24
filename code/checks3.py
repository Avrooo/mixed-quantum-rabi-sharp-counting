import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from cert import *
import time
for n in (2,3,4):
    for t in (F(1,50),F(1,20)):
        for e in (F(-10**2),F(-10**4),F(-10**6)):
            t0=time.time(); r=certify(n,t,e)
            print(f"n={n} t={t} eps={float(e):.0e}: B_n={((n+1)**2)//4} certified count={r.get('count')} status={r.get('status')} [{time.time()-t0:.1f}s]",flush=True)
