"""f(t)=m^2 q /(r^2 (q-8t)) is strictly increasing on (-1/2,0)  <=> lambda_+ = (1-sqrt f)/2 strictly decreasing from 1/2 to 0.
f' numerator computed exactly; sign checked by exact Sturm count (no root in (-1/2,0)) and sign at one point."""
from fractions import Fraction as F
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from audit_divisibility import mul,add,sc,trim,sturm_count
t=[F(0),F(1)]; one=[F(1)]
m=add(one,sc(t,2)); r=add(one,sc(t,-2)); q=add(one,sc(mul(t,t),4)); d=add(q,sc(t,-8))
def der(p): return [i*p[i] for i in range(1,len(p))]
Nn=mul(mul(m,m),q); Dn=mul(mul(r,r),d)
num=trim(add(mul(der(Nn),Dn),sc(mul(Nn,der(Dn)),-1)))
print("numerator of f' (low->high):",[str(c) for c in num])
cnt=sturm_count(num,F(-1,2)+F(1,10**12),F(0)-F(1,10**12)) if False else None
# strip trivial factors (1+2t),(1-2t)? just count roots strictly inside
def ev(p,x):
    s=F(0)
    for c in reversed(p): s=s*x+c
    return s
cnt=sturm_count(num,F(-1,2)+F(1,10**15),F(-1,10**15))
print("distinct real roots of numerator in (-1/2,0):",cnt,"; sign at t=-1/4:",ev(num,F(-1,4))>0, "(f' sign = numerator/(denominator^2>0))")
