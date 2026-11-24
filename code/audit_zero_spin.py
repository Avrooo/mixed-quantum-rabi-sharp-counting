"""Zero squared-spin eigenvalue (k=0) and the line w=0.  Exact.
 (a) m even, U=2(m+2): A^_{m+1}(2s,c) is odd in c, its root c_j=0 (j=s/2, k=0) is simple and is exactly the removed parity factor:
     P_{m+1}(2s,w) = C * prod_{k>0, k=m mod 2}((s^2-k^2)w-8k^2) with C!=0 and P_{m+1}(2s,0)!=0   (m<=24).
 (b) the stripped pair has NO common zero on the line w=0 for any U  (n<=24): gcd_U(P_{n+1}(U,0),(U-2n-4)P_n(U,0)) = 1."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from audit_limit_bezout import build, ev
from ratfunc import pgcd, trim, pmul
P=build(26); ok=True
for m in range(1,25):
    s=m+2; U=F(2*s)
    # polynomial in w: evaluate P_{m+1}(2s,w) at w=0..m+1 and compare with product (leading coefficient fit)
    # product polynomial
    prod=[F(1)]
    for k in range(1,m+1):
        if (k-m)%2==0: prod=pmul(prod,[F(-8*k*k),F(s*s-k*k)])
    # P_{m+1}(2s,w) as polynomial: coefficients from dict
    co={}
    for (i,j),v in P[m+1].items(): co[j]=co.get(j,0)+v*U**i
    pol=trim([co.get(j,F(0)) for j in range(max(co)+1)])
    ratio=pol[-1]/prod[-1]
    same=len(pol)==len(prod) and all(pol[i]==ratio*prod[i] for i in range(len(pol)))
    nz0=pol[0]!=0
    ok&=same and nz0
    if m<=8 or m%8==0: print("m=%2d: P_{m+1}(2s,w)=C*prod_{k>0}((s^2-k^2)w-8k^2): %s ; deg=%d=ceil(m/2) : %s ; P(2s,0)!=0 (k=0 is not a root of the stripped polynomial): %s"%(m,same,len(pol)-1,len(pol)-1==(m+1)//2,nz0))
for n in range(1,25):
    a=[]; 
    co={}
    for (i,j),v in P[n+1].items():
        if j==0: co[i]=co.get(i,0)+v
    a=trim([co.get(i,F(0)) for i in range(max(co)+1)])
    co={}
    for (i,j),v in P[n].items():
        if j==0: co[i]=co.get(i,0)+v
    b=trim([co.get(i,F(0)) for i in range(max(co)+1)]); b=pmul(b,[F(-(2*n+4)),F(1)])
    g=pgcd(a,b); ok&=len(g)==1
print("(b) no common zero of the stripped pair on w=0, n<=24:",ok)
print("ALL PASS" if ok else "FAILURE")
