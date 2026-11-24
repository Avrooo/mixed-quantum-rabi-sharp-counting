"""Exact: A numerator NA(t) has exactly one root in (0,1/2), located to 1e-12; V>0, Q monotone on (-1/2,0), Q(-1/2)=-3, Q(0)=1."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from ratfunc import pdivmod,pgcd,trim
NA=[F(c) for c in (9,36,-24,-432,-96,-1728,-384,2304,2304)]
def seq(P):
    s=[trim(P),trim([i*c for i,c in enumerate(P)][1:])]
    while len(s[-1])>1:
        r=pdivmod(s[-2],s[-1])[1]
        if not r: break
        s.append([-c for c in r])
    return s
def ev(p,x): return sum((c*x**i for i,c in enumerate(p)),F(0))
S=seq(NA); g=pgcd(NA,[i*c for i,c in enumerate(NA)][1:]); print("gcd(NA,NA') degree:",len(g)-1,"(squarefree)" if len(g)==1 else "")
var=lambda x:sum(1 for a,b in zip([ev(p,x) for p in S if ev(p,x)!=0],[ev(p,x) for p in S if ev(p,x)!=0][1:]) if a*b<0)
print("exact number of roots of NA in (0,1/2):",var(F(0))-var(F(1,2)))
lo,hi=F(0),F(1,2)
for _ in range(60):
    mid=(lo+hi)/2
    if ev(NA,lo)*ev(NA,mid)<=0: hi=mid
    else: lo=mid
print("root t0 = %.12f  (g2/omega = %.5f)"%(float(lo),float(2*lo/(1+4*lo*lo))))
Qp=lambda t:6+8*t+72*t**2+64*t**3
print("Q(-1/2)=",1+6*F(-1,2)+4*F(1,4)+24*F(-1,8)+16*F(1,16)," Q(0)=1;  min of Q' on a fine grid of (-1/2,0):",float(min(Qp(F(-i,2000)) for i in range(0,1001))))
V=lambda t:16*t**4+16*t**3+24*t**2+4*t+1
print("min V on grid (-1/2,1/2):",float(min(V(F(i,2000)) for i in range(-1000,1001))),"; lower bound 16t^2+4t+1 has discriminant",16-64)
