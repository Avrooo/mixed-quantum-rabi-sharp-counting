"""n=1, t=t_p (root of Q in (-1/2,0)): the system degenerates to one line in (u,v); check its signs at 60+ digits."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import ratrec
def Qf(t): return 1+6*t+4*t**2+24*t**3+16*t**4
lo,hi=F(-1,2),F(0)
for _ in range(220):
    mid=(lo+hi)/2
    if Qf(lo)*Qf(mid)<=0: hi=mid
    else: lo=mid
tp=(lo+hi)/2
print("t_p =",float(tp)," g2 =",float(2*tp/(1+4*tp**2)))
c1,c0=ratrec.reduced_coeffs(1,tp)
g=lambda c,k: c.get(k,F(0))
# keys (i_u, i_eps, i_v)
row=lambda c:(g(c,(1,0,0)),g(c,(0,0,1)),g(c,(0,1,0)),g(c,(0,0,0)))
a1,b1,e1,k1=row(c1); a0,b0,e0,k0=row(c0)
det=a1*b0-b1*a0
print("det of (u,v) block at t_p: %.3e   (scale %.3e)"%(float(det),float(abs(a1*b0))))
lam=a0/a1; print("row0/row1 ratio in (u,v) block consistent:",float(abs(b0-lam*b1)))
eps=-(k0-lam*k1)/(e0-lam*e1); print("eps* =",float(eps))
A_=a1; B_=b1; R=-(k1+e1*eps)
print("line: %.6f u + %.6f v = %.6f"%(float(A_),float(B_),float(R)))
print("both coefficients > 0:",A_>0 and B_>0," RHS < 0:",R<0,"  => no (u,v)>0 solution:", (A_>0 and B_>0 and R<0))
print("margin |R| =%.4f vs perturbation ~1e-60 (t_p known to 66 digits)"%float(abs(R)))
