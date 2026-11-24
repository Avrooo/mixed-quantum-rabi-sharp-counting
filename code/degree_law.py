import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
"""Exact verification of the weighted-degree law and the leading-coefficient formulas (proof ingredients)."""
import sys; from fractions import Fraction as F
import exact as ex
out=[]
def wt(k): return k[0]*1+k[1]*2+k[2]*1   # key=(Delta exp, eps exp, g exp): weights Delta=1, eps=2, g=1
allok=True
for t in (F(1,10),F(3,10),F(2,5),F(-1,10),F(-3,10),F(1,50)):
    q,r,m,Om,g2,ka=ex.setup(t)
    for n in range(1,15):
        R1,R0,E=ex.constraints_qp(n,t)
        w1=max(wt(k) for k in R1.d); w0=max(wt(k) for k in R0.d)
        ok=(w1==n+2 and w0==n+3)
        # leading monomial coefficients predicted by the linearised recursion
        rho=-1/(2*t)
        if n%2==0:
            pr1=(1/(4*t))*rho**((n-2)//2)*(n+2)*q*m/(2*t*r*r); k1=(1,n//2,1)
            pr0=-(1/(2*t))*rho**(n//2); k0=(1,n//2+1,0)
        else:
            pr1=-(1/(2*t))*rho**((n-1)//2); k1=(1,(n+1)//2,0)
            pr0=(1/(2*t))*rho**((n-1)//2)*(q/r)*(F(n+1)/(4*t)+F(n+2)/r); k0=(1,(n+1)//2,1)
        c1=R1.d.get(k1,0); c0=R0.d.get(k0,0)
        ok=ok and c1==pr1 and c0==pr0 and wt(k1)==n+2 and wt(k0)==n+3
        # parity: only monomials Delta^(odd) g^(parity) appear in the right pattern
        par1=all((k[0]%2==1) and (k[2]%2==(n+1)%2) for k in R1.d); par0=all((k[0]%2==1) and (k[2]%2==n%2) for k in R0.d)
        ok=ok and par1 and par0
        allok&=ok
        if not ok: print("FAIL",t,n,w1,w0,c1==pr1,c0==pr0,par1,par0)
print("weighted degree (Delta=1,g1=1,eps=2): max weight R1=n+2, R0=n+3; leading-monomial coefficients equal the closed forms; parity pattern holds -> all pass:",allok,"(n=1..14, 6 values of alpha)")
# where does the odd-n R0 leading monomial vanish for alpha<0?
for n in (1,3,5,7):
    print("n=%d: R0 key monomial coefficient vanishes at alpha=%s"%(n,-F(n+1,2*(n+3))))
