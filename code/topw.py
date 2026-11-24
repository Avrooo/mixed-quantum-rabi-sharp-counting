import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
# top-weight (leading-form) recursion, exact, derived by hand from exact.core (weights: x=g=1, eps=2)
import sys; from fractions import Fraction as F
import exact as ex
from limit_sys import leading_forms
def top(n,t,x,g,e):
    q,r,m,Om,g2,ka=ex.setup(t); p=4*t
    beta=ka*g; h=g+2*beta*g2; s=g*beta+g2*beta*beta; sig=s+e
    A={-2:F(0),-1:F(0)}; B={}
    B[0]=F(1); A[0]=-x/p
    for i in range(1,n+3):
        B[i]=(x*A[i-2]-(n+3-i)*h*B[i-1])/(i*Om)
        A[i]=-(2*beta*A[i-1]+2*sig*A[i-2]+x*B[i])/p
    R1=-p*A[n+1]; R0=2*sig*A[n]+x*B[n+2]
    return R1,R0,A,B
def ev(form,u,e,v): return sum(c*u**i*e**j*v**k for (i,j,k),c in form.items())
if __name__=="__main__":
    import random; random.seed(1)
    for t in (F(1,10),F(-1,7),F(2,5)):
        for n in (1,2,3,4,5,6,7):
            F1,F0=leading_forms(n,t)
            rat1=set();rat0=set()
            for _ in range(5):
                x=F(random.randint(1,9),random.randint(1,5)); g=F(random.randint(1,9),random.randint(1,5)); e=F(random.randint(-9,9),random.randint(1,4))
                R1,R0,_,_=top(n,t,x,g,e)
                par1=(n+1)%2; par0=n%2
                f1=ev(F1,x*x,e,g*g); f0=ev(F0,x*x,e,g*g)
                rat1.add(R1/(x*g**par1)/f1); rat0.add(R0/(x*g**par0)/f0)
            print(f"t={t} n={n}: constant ratio R1_top/(Delta g^par F1): {len(rat1)==1}, R0: {len(rat0)==1}")
