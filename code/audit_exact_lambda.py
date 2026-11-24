import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
"""Independent exact audit of the coprimality theorem (V4 / Theorem 3).
For rational t with q(q-8t) a rational square, lambda_+(t) is rational, so the predicted
exceptional condition [reduced denominator of lambda_+ <= n+2] and the predicted common root
chi = u/v = -4 m w2 M'/r^3 can be tested in EXACT arithmetic (no floats, no rootfinding)."""
import sys; from fractions import Fraction as F
from math import isqrt
from gforms import inf_forms, binres
def fsqrt(x):
    if x<0: return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return F(a,b) if a*a==x.numerator and b*b==x.denominator else None
def params(t):
    q=1+4*t*t; r=1-2*t; m=1+2*t; w2=1-2*t+4*t*t
    return q,r,m,w2
def lam_plus(t):
    q,r,m,w2=params(t); s=fsqrt(q/(q-8*t)) if q-8*t!=0 else None
    return None if s is None else F(1,2)-m/(2*r)*s
def evalform(G,u,v): return sum(c*u**i*v**k for (i,k),c in G.items())
found=[]
for b in range(2,60):
    for a in range(1,b):
        t=F(-a,b)
        if abs(t)>=F(1,2): continue
        L=lam_plus(t)
        if L is not None and 0<L<F(1,2) and L.denominator<=9: found.append((t,L))
found.sort(key=lambda z:(z[1].denominator,z[1]))
print("rational t<0 with rational lambda_+ (denominator<=9):",[(str(t),str(L)) for t,L in found])
ok=True
for t,L in found:
    q,r,m,w2=params(t); Mp=L.denominator
    chi=-4*m*w2*Mp/r**3
    for n in range(1,10):
        (d1,G1),(d0,G0)=inf_forms(n,t)
        res=binres(G1,d1,G0,d0)
        pred = Mp<=n+2
        # common root (chi:1) check
        v1=evalform(G1,chi,F(1)); v0=evalform(G0,chi,F(1))
        ok &= ((res==0)==pred)
        if pred: ok &= (v1==0 and v0==0)
        print(f"t={t} lambda+={L} M'={Mp} n={n}: Res==0:{res==0} predicted:{pred}  G1(chi,1)==0:{v1==0} G0(chi,1)==0:{v0==0}")
print("ALL PREDICTIONS CONFIRMED EXACTLY:",ok)
# positive t: lambda_+ never a positive rational => no exceptional t. Test rational t>0 with q(q-8t) square.
pos=[]
for b in range(2,200):
    for a in range(1,b):
        t=F(a,b)
        if t>=F(1,2): continue
        L=lam_plus(t)
        if L is not None: pos.append((t,L))
print("rational t in (0,1/2) with rational lambda_+:",[(str(t),str(L)) for t,L in pos][:10],"; all lambda_+<0:",all(L<0 for t,L in pos))
