import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from fractions import Fraction as F
from math import isqrt
from gforms import inf_forms, binres
D=10**80
tstar=F(D-isqrt(3*D*D)//2*1,D)  # 1 - sqrt(3)/2 approx
s3=F(isqrt(3*D*D),D); tstar=1-s3/2
q=1+4*tstar**2
print("q-8t at t*:",float(q-8*tstar))
def rel(n,t):
    (d1,G1),(d0,G0)=inf_forms(n,t)
    m1=max(abs(float(c)) for c in G1.values()); m0=max(abs(float(c)) for c in G0.values())
    return abs(float(binres(G1,d1,G0,d0)))/(m1**d0*m0**d1) if True else 0
for n in range(1,7):
    print(n,"rel|Res| at t*:",f"{rel(n,tstar):.2e}","  at t*+1e-3:",f"{rel(n,tstar+F(1,1000)):.2e}","  at t*+1e-6:",f"{rel(n,tstar+F(1,10**6)):.2e}",flush=True)
