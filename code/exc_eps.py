import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys,math; from decimal import Decimal, getcontext
getcontext().prec=80
from fractions import Fraction as F
from crit import fit, eliminants, real_roots, swap, peval
import geom
def lam_dec(t):
    q=1+4*t*t; r=1-2*t; w2=1-2*t+4*t*t
    s0=(1-8*t/q).sqrt()
    return w2*(s0-1)/(s0*(r*r-(r*q/2)*(s0-1)))
def solve_t(target,lo=Decimal("-0.4999"),hi=Decimal("-0.0000001")):
    f=lambda t:lam_dec(t)-Decimal(target.numerator)/Decimal(target.denominator)
    for _ in range(260):
        mid=(lo+hi)/2
        if f(lo)*f(mid)<=0: hi=mid
        else: lo=mid
    return F(str((lo+hi)/2)).limit_denominator(10**45)
def pad(R,k): return R[k] if k<len(R) else F(0)
def analyse(n,target):
    t=solve_t(target); c1,c0,mk=eliminants(n,t)
    # coefficient polynomials in eps of Res_v (eliminating v -> polynomial in u): list of coefficient fits
    m_u=(lambda e: mk(e)[1])
    deg=max(len(m_u(F(k))) for k in range(0,3))
    polys=[fit(lambda e,i=i: pad(m_u(e),i),maxd=60) for i in range(deg)]
    out=[]
    for i,p in enumerate(polys):
        if p is None or len(p)<2: out.append((i,"const/zero" if p is None or len(p)<2 else "")); continue
        out.append((i,len(p)-1))
    # common real roots: roots of the first nonconstant poly, tested on the others
    cand=None
    for p in polys:
        if p is not None and len(p)>1: cand=p;break
    res=[]
    if cand is not None:
        L=F(10**4)
        for r in real_roots(cand,-L,L):
            rr=F(r).limit_denominator(10**14)
            vals=[abs(float(peval(p,rr))) for p in polys if p is not None]
            sc=[max(abs(float(c))for c in p) for p in polys if p is not None]
            res.append((float(rr),max(v/s for v,s in zip(vals,sc))))
    return t,[len(p)-1 if p else None for p in polys],res
if __name__=="__main__":
    for n,target in ((1,F(1,3)),(2,F(1,3)),(2,F(1,4))):
        t,degs,res=analyse(n,target)
        print(f"n={n} lambda+={target} t_j={float(t):.12f}  deg of res-coeff polys in eps: {degs}")
        print("   candidate eps (root of first coeff poly) and max relative size of the other coeff polys there:",[(round(a,6),f"{b:.1e}") for a,b in res],flush=True)
