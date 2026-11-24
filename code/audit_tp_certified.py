"""n=1, t=t_p: RIGOROUS (exact rational interval arithmetic) proof that the degenerate line has no admissible point.
Steps: (1) symbolically in t: the (u,v) determinant of the 2x2 linear system equals Q(t) times a nonvanishing factor,
(2) bracket the unique root t_p of Q in (-1/2,0) by exact rationals [tl,th] with Q(tl)<0<Q(th) (Q increasing there),
(3) evaluate the rational functions a1(t), b1(t), R(t) with interval arithmetic over that bracket; the true value at t_p is
contained in the enclosure. Claim: a1>0, b1>0, R<0 on the bracket => a1 u + b1 v = R has no solution with u,v>0."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from ratfunc import RF,T,pmul,pdivmod
import ratrec
c1,c0=ratrec.reduced_coeffs(1,T)          # dict of RF coefficients, symbolic in t
g=lambda c,k: c.get(k,RF(0)) if not isinstance(c.get(k,0),(int,F)) else RF(c.get(k,0))
row=lambda c:(g(c,(1,0,0)),g(c,(0,0,1)),g(c,(0,1,0)),g(c,(0,0,0)))
a1,b1,e1,k1=row(c1); a0,b0,e0,k0=row(c0)
det=a1*b0-b1*a0
Qp=RF([1,6,4,24,16])
ratio=det/Qp
print("det/Q(t) = RF with numerator deg %d, denominator deg %d"%(len(ratio.n)-1,len(ratio.d)-1))
# numerator of det must contain Q: ratio has no pole at roots of Q  <=> (det/Q).d is coprime to Q
from ratfunc import pgcd
print("det vanishes at roots of Q (Q | numerator of det):", len(pgcd(det.n,[F(1),F(6),F(4),F(24),F(16)]))==5)
print("gcd(den of det/Q, Q) degree:",len(pgcd(ratio.d,[F(1),F(6),F(4),F(24),F(16)]))-1,"(0 means no pole at roots of Q)")
# line quantities as functions of t (valid where a1 != 0)
lam=a0/a1
eps=-(k0-lam*k1)/(e0-lam*e1)
R=-(k1+e1*eps)
# --- interval arithmetic
class I:
    def __init__(s,lo,hi=None): s.lo=F(lo); s.hi=F(hi if hi is not None else lo)
    def __add__(s,o): o=I.c(o); return I(s.lo+o.lo,s.hi+o.hi)
    def __mul__(s,o):
        o=I.c(o); p=[s.lo*o.lo,s.lo*o.hi,s.hi*o.lo,s.hi*o.hi]; return I(min(p),max(p))
    @staticmethod
    def c(o): return o if isinstance(o,I) else I(o)
def ipoly(p,x):
    r=I(0)
    for c in reversed(p): r=r*x+I(c)
    return r
def ieval(rf,x):
    N=ipoly(rf.n,x); D=ipoly(rf.d,x)
    assert D.lo>0 or D.hi<0, "denominator interval contains 0"
    inv=I(1/D.hi,1/D.lo) if D.lo>0 else I(1/D.lo,1/D.hi)
    return N*inv
Qf=lambda t: 1+6*t+4*t**2+24*t**3+16*t**4
# Q' > 0 on (-1/2,0): Q' = 6+8t+72t^2+64t^3 >= 6+8t+40t^2 >0 (disc 64-960<0) -- so unique root; bisect with exact rationals
lo,hi=F(-1,2),F(0)
assert Qf(lo)<0<Qf(hi)
for _ in range(100):
    mid=(lo+hi)/2
    if Qf(mid)<0: lo=mid
    else: hi=mid
assert Qf(lo)<0<Qf(hi)
X=I(lo,hi); print("t_p bracket width: %.3e   [%.15f, %.15f]"%(float(hi-lo),float(lo),float(hi)))
A1,B1,RR=ieval(a1,X),ieval(b1,X),ieval(R,X)
print("a1 in [%.6f,%.6f]  b1 in [%.6f,%.6f]  R in [%.6f,%.6f]"%tuple(float(v) for v in (A1.lo,A1.hi,B1.lo,B1.hi,RR.lo,RR.hi)))
ok=A1.lo>0 and B1.lo>0 and RR.hi<0
print("CERTIFIED: a1>0, b1>0, R<0 on the whole bracket:",ok)
assert ok
den=e0-lam*e1
Dn=ieval(den,X); print("e0 - (a0/a1) e1 in [%.6f,%.6f]  (nonzero on bracket: %s)"%(float(Dn.lo),float(Dn.hi),Dn.lo>0 or Dn.hi<0))
assert Dn.lo>0 or Dn.hi<0
