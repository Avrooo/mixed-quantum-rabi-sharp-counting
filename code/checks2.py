import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from topw import top
from limitpoly import chain, peval2
from fractions import Fraction as F
A=chain(14)
def lim(i,U,c): return sum(v*U**a*c**b for (a,b),v in A[i].items())
# C7: tilde values from the exact top-weight recursion converge to the limit recurrence as tau->0
U=F(9,4); c=F(7,5); x=F(3,2)   # x^2=U
for tau in (F(1,100),F(1,1000),F(1,10000)):
    t=tau*tau; g=tau*c
    worst=0
    for n in (3,6):
        R1,R0,Ai,Bi=top(n,t,x,g,F(-1))
        A0=Ai[0]
        for i in range(0,n+3):
            rat=Ai[i]/A0*(tau**i)
            worst=max(worst,float(abs(rat/lim(i,U,c)-1)) if lim(i,U,c)!=0 else 0)
    print(f"C7 tau={float(tau):.0e}: max rel dev of A_i/A_0 * tau^i from limit recurrence (i<=n+2, n=3,6): {worst:.2e}")
