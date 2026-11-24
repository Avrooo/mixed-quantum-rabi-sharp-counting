"""Exact audit: the resultant numerator N_n(t) (interpolated exactly in V4) must be divisible by the
quartic  P_{k/M}(t) = (M-2k)^2 r^2 (q-8t) - M^2 m^2 q   for every reduced k/M<1/2, M<=n+2
(this encodes lambda_+(t)=k/M, lambda_+ = 1/2 - (m/2r) sqrt(q/(q-8t)) ), and N_n has no other factor with a root in (-1/2,1/2)."""
import pickle
from fractions import Fraction as F
from math import gcd
def trim(p):
    p=list(p)
    while p and p[-1]==0: p.pop()
    return p
def mul(a,b):
    r=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def add(a,b):
    n=max(len(a),len(b)); return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)]
def sc(a,c): return [c*x for x in a]
def divmod_(a,b):
    a=trim(a); b=trim(b); q=[F(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and a:
        c=a[-1]/b[-1]; k=len(a)-len(b); q[k]=c
        a=trim([a[i]-(c*b[i-k] if 0<=i-k<len(b) else 0) for i in range(len(a))])
    return trim(q),a
def Pk(k,M):
    t=[F(0),F(1)]; one=[F(1)]
    q=add(one,sc(mul(t,t),4)); r=add(one,sc(t,-2)); m=add(one,sc(t,2))
    lhs=sc(mul(mul(r,r),add(q,sc(t,-8))),(M-2*k)**2)
    rhs=sc(mul(mul(m,m),q),M*M)
    return trim(add(lhs,sc(rhs,-1)))
def farey(n):
    return sorted({F(k,M) for M in range(2,n+3) for k in range(1,M) if 2*k<M and gcd(k,M)==1})
def sturm_count(P,lo,hi):
    def deriv(p): return [i*p[i] for i in range(1,len(p))]
    seq=[trim(P),trim(deriv(P))]
    while len(seq[-1])>1:
        _,rem=divmod_(seq[-2],seq[-1]); 
        if not rem: break
        seq.append(sc(rem,F(-1)))
    def ev(p,x): 
        s=F(0)
        for c in reversed(p): s=s*x+c
        return s
    def var(x):
        v=[ev(p,x) for p in seq]; v=[a for a in v if a!=0]
        return sum(1 for i in range(len(v)-1) if v[i]*v[i+1]<0)
    return var(lo)-var(hi)

allok=True
print("n | multiplicities of each quartic P_{k/M} in N_n | residual after removing all P and q-powers")
for n in range(1,7):
    a,c1,c2,(b,N)=pickle.load(open(f"data/resN_{n}.pkl","rb")); cof=trim(N)
    mult={}
    for fr in farey(n):
        P=Pk(fr.numerator,fr.denominator); k=0
        while True:
            qq,rem=divmod_(cof,P)
            if rem: break
            cof=qq; k+=1
        mult[str(fr)]=k; allok&=(k>=1)
    qk=0
    while True:
        qq,rem=divmod_(cof,[F(1),F(0),F(4)])
        if rem: break
        cof=qq; qk+=1
    print(f"n={n}: mult {mult}; q^{qk}; residual degree {len(cof)-1}")
    allok&=(len(cof)==1)
print("EXACT DIVISIBILITY AUDIT PASSED (N_n = const * prod P_{k/M}^{mult} * q^j, no other factor):",allok)
