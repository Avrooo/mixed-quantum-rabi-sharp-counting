"""Exact check of Prop. n=1 (ii): for t in (-1/2,0) the positivity conditions Delta^2=A-B eps>0, g1^2=C-D eps>0 are incompatible."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from ratfunc import RF,T
import ratrec
A,B,C,D=ratrec.ABCD(T)
import audit_symbolic_n1 as S   # re-runs symbolic check (cheap)
def sg(x): return (x>0)-(x<0)
tp=None
print("t        sgnQ sgnB sgnD sgnA sgnC  feasible? (exists eps: A-B eps>0, C-D eps>0)")
def feasible(a,b,c,d):
    # intervals for eps: b*eps<a ; d*eps<c
    lo,hi=-float('inf'),float('inf')
    for (k,v) in ((b,a),(d,c)):
        if k>0: hi=min(hi,v/k)
        elif k<0: lo=max(lo,v/k)
        elif v<=0: return False
    return lo<hi
bad=[]
N=2000
for i in range(1,N):
    t=F(-1,2)+F(i,2*N)
    if t>=0: break
    a,b,c,d=A(t),B(t),C(t),D(t)
    if feasible(a,b,c,d): bad.append(float(t))
print("feasible rational t sample in (-1/2,0):",len(bad),"of",N//1,"; first few:",bad[:5])
# exact: find sign pattern at roots
Wn=S.W.n
import math
def fl(p,x): return sum(float(c)*x**i for i,c in enumerate(p))
# bisect roots of Q and W in (-1/2,0)
def bis(p,lo,hi):
    for _ in range(80):
        mid=(lo+hi)/2
        if fl(p,lo)*fl(p,mid)<=0: hi=mid
        else: lo=mid
    return lo
print("root of Q:",bis([1,6,4,24,16],-.5+1e-9,-1e-9) , " root of W:", bis(Wn,-.5+1e-9,-1e-9))
for t in (F(-49,100),F(-3,10),F(-1,5),F(-17,100),F(-1,6),F(-1,10),F(-1,100)):
    a,b,c,d=A(t),B(t),C(t),D(t); print(float(t),sg(float(RF([1,6,4,24,16])(t))),sg(b),sg(d),sg(a),sg(c),feasible(float(a),float(b),float(c),float(d)))
