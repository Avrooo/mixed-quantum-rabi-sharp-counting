"""Exact audit of the limit system (Theorem 7) in the form used by the Bezout-based proof:
 stripped polynomials  P_i(U,w):  A^_i = c^{i mod 2} P_i,  P_i=-(1/2)P_{i-1}+k_i P_{i-2} (i odd), -(w/2)P_{i-1}+k_i P_{i-2} (i even), k_i=(2i-U)/(4i).
 Stripped pair: Phi1=P_{n+1}, Phi0=(U-2n-4)P_n.   Check, exactly, for n<=24:
  (i)   total degrees d1,d0;     (ii) the top-degree forms have no common root in P^1 (resultant != 0);
  (iii) the B_n points of Z_n are zeros of both polynomials (rational arithmetic), pairwise distinct, w>0;
  (iv)  Jacobian d(Phi1,Phi0)/d(U,w) != 0 at each (consistent with Bezout: total multiplicity B_n = number of points)."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from gforms import binres
def padd(a,b):
    r=dict(a)
    for k,v in b.items():
        r[k]=r.get(k,0)+v
        if r[k]==0: del r[k]
    return r
def pscal(a,c): return {k:v*c for k,v in a.items()} if c!=0 else {}
def mulU(a): return {(i+1,j):v for (i,j),v in a.items()}
def mulw(a): return {(i,j+1):v for (i,j),v in a.items()}
def ev(a,U,w): return sum((v*U**i*w**j for (i,j),v in a.items()),F(0))
def evd(a,U,w,which):
    s=F(0)
    for (i,j),v in a.items():
        if which=='U' and i>0: s+=v*i*U**(i-1)*w**j
        if which=='w' and j>0: s+=v*j*U**i*w**(j-1)
    return s
def build(N):
    P=[{(0,0):F(1)}]
    for i in range(1,N+1):
        # k_i P_{i-2} = ((2i-U)/(4i)) P_{i-2}
        prev2=P[i-2] if i>=2 else {}
        kt=padd(pscal(prev2,F(2*i,4*i)),pscal(mulU(prev2),F(-1,4*i)))
        first=pscal(P[i-1],F(-1,2)) if i%2==1 else pscal(mulw(P[i-1]),F(-1,2))
        P.append(padd(first,kt))
    return P
P=build(26)
allok=True
for n in range(1,25):
    d1=(n+1)//2; d0=n//2+1; Bn=d1*d0
    Phi1=P[n+1]; Phi0=padd(mulU(P[n]),pscal(P[n],F(-(2*n+4))))
    deg=lambda a:max(i+j for (i,j) in a)
    e1=deg(Phi1)==d1; e0=deg(Phi0)==d0
    G1={k:v for k,v in Phi1.items() if sum(k)==d1}; G0={k:v for k,v in Phi0.items() if sum(k)==d0}
    R=binres(G1,d1,G0,d0)
    pts=[(F(2*(m+2)),F(8*k*k,(m+2)**2-k*k)) for m in range(1,n+1) for k in range(m%2 if m%2 else 2,m+1,2) if k>0]
    # k ≡ m mod 2, 0<k<=m
    pts=[(F(2*(m+2)),F(8*k*k,(m+2)**2-k*k)) for m in range(1,n+1) for k in range(1,m+1) if (k-m)%2==0]
    zero=all(ev(Phi1,U,w)==0 and ev(Phi0,U,w)==0 for U,w in pts)
    distinct=len(set(pts))==len(pts)==Bn and all(w>0 for U,w in pts)
    jac=all(evd(Phi1,U,w,'U')*evd(Phi0,U,w,'w')-evd(Phi1,U,w,'w')*evd(Phi0,U,w,'U')!=0 for U,w in pts)
    ok=e1 and e0 and R!=0 and zero and distinct and jac
    allok&=ok
    print("n=%2d  deg=(%d,%d) exact:%s  Res(top forms)!=0:%s  B_n=%d points all zeros:%s distinct,w>0:%s  Jacobians!=0:%s"%(n,d1,d0,e1 and e0,R!=0,Bn,zero,distinct,jac),flush=True)
print("ALL PASS" if allok else "FAILURE")
