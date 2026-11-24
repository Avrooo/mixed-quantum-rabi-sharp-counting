import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from exc_eps import solve_t
from gforms import inf_forms, binres
from fractions import Fraction as F
from math import gcd
def farey(n):
    S=set()
    for M in range(2,n+3):
        for k in range(1,M):
            if 2*k<M and gcd(k,M)==1: S.add(F(k,M))
    return sorted(S)
out=[]
for n in range(1,7):
    for lam in farey(n):
        t=solve_t(lam); (d1,G1),(d0,G0)=inf_forms(n,t)
        m1=max(abs(float(c)) for c in G1.values()); m0=max(abs(float(c)) for c in G0.values())
        res=binres(G1,d1,G0,d0); scale=m1**d0*m0**d1
        # leading/trailing coefficients (points [1:0] and [0:1])
        print(f"n={n} lambda+={lam} t={float(t):.10f}: G1 nonzero:{bool(G1)} G0 nonzero:{bool(G0)} relative|Res|={abs(float(res))/scale:.1e}  G1(1,0)!=0:{G1.get((d1,0),0)!=0} G0(1,0)!=0:{G0.get((d0,0),0)!=0} G1(0,1)!=0:{G1.get((0,d1),0)!=0} G0(0,1)!=0:{G0.get((0,d0),0)!=0}",flush=True)
