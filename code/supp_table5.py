"""Full-precision (40-digit) parameters of Table V: refines (Delta,g1) by Newton iteration in 50-digit decimal arithmetic."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from decimal import Decimal as Dc, getcontext
import ratrec, geom
getcontext().prec=60
t=F(1,10)
def Dev(P,u,v,du=0,dv=0):
    s=Dc(0)
    for (i,k),c in P.items():
        if du and i<du: continue
        if dv and k<dv: continue
        co=Dc(c.numerator)/Dc(c.denominator); a_=i; b_=k
        if du: co*=a_; a_-=1
        if dv: co*=b_; b_-=1
        s+=co*u**a_*v**b_
    return s
def refine(n,eps,D0,g0):
    c1,c0=ratrec.reduced_coeffs(n,t); P=geom.fix_eps(c1,eps); Q=geom.fix_eps(c0,eps)
    u=Dc(D0)**2; v=Dc(g0)**2
    for _ in range(100):
        f1=Dev(P,u,v); f0=Dev(Q,u,v)
        a11=Dev(P,u,v,du=1); a12=Dev(P,u,v,dv=1); a21=Dev(Q,u,v,du=1); a22=Dev(Q,u,v,dv=1)
        det=a11*a22-a12*a21; du=(-f1*a22+a12*f0)/det; dv=(-a11*f0+a21*f1)/det; u+=du; v+=dv
        if abs(du)+abs(dv)<Dc(10)**-45: break
    return u,v,abs(f1)+abs(f0)
tq=Dc(1)+4*Dc('0.1')**2; r=1-2*Dc('0.1'); m=1+2*Dc('0.1'); Om=r*m/tq; g2=2*Dc('0.1')/tq
states=[(1,F(0),2.5702,0.27392),(1,F(1,2),2.0288,0.25141),(2,F(0),2.4522,0.24538),(3,F(0),4.4847,0.2166),(4,F(1,2),4.7092,0.3860),(5,F(1,2),5.8404,0.1785),(6,F(0),4.1409,0.1797)]
print("n,eps,Delta,g1,E,residual_of_constraints")
for n,e,D0,g0 in states:
    u,v,res=refine(n,e,D0,g0); eps=Dc(e.numerator)/Dc(e.denominator)
    E=(n+Dc('2.5'))*Om-Dc('0.5')-v/(1-2*g2)-eps
    print("%d,%s,%.40s,%.40s,%.40s,%.2e"%(n,e,u.sqrt(),v.sqrt(),E,float(res)))
