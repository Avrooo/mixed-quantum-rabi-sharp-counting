import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys,random; from fractions import Fraction as F
import exact as ex
from topw import top
from gforms import inf_forms, binres
random.seed(7)
# (1) generating-function ODE / 4-term recurrence for the top-weight series A_i (g=1, x^2=chi, eps arbitrary)
def series(n,t,x,g,e,K):
    q,r,m,Om,g2,ka=ex.setup(t); p=4*t; beta=ka*g; h=g+2*beta*g2; s=g*beta+g2*beta*beta; sig=s+e; N=n+3
    A={-2:F(0),-1:F(0)}; B={0:F(1)}; A[0]=-x/p
    for i in range(1,K):
        B[i]=(x*A[i-2]-(N-i)*h*B[i-1])/(i*Om); A[i]=-(2*beta*A[i-1]+2*sig*A[i-2]+x*B[i])/p
    L0=Om*p; L1=2*beta*Om-h*p; L2=2*sig*Om-2*beta*h; L3=-2*sig*h
    M0=2*beta*Om+h*(N-1)*p; M1=4*sig*Om-2*beta*h+2*beta*h*(N-1)+x*x; M2=2*sig*h*(N-3)
    return A,(L0,L1,L2,L3),(M0,M1,M2),N,sig,h
ok=True; cnt=0
for t in (F(1,10),F(-1,7),F(3,10),F(-2,5),F(1,9)):
    for n in (1,2,3,5,8):
        x=F(random.randint(1,9),random.randint(1,4)); g=F(random.randint(1,9),random.randint(1,4)); e=F(random.randint(-9,9),3)
        K=n+9; A,L,M,N,sig,h=series(n,t,x,g,e,K)
        for i in range(0,K-2):
            val=L[0]*(i+1)*A[i+1]+(L[1]*i+M[0])*A[i]+(L[2]*(i-1)+M[1])*A[i-1]+(L[3]*(i-2)+M[2])*A[i-2]
            cnt+=1
            if val!=0: ok=False
print(f"S1 four-term recurrence from L A'+M A=0 holds exactly for the top-weight series ({cnt} checks, 5 t, n in 1,2,3,5,8):",ok)
# (2) quadratic-irrational identities in Q(sqrt D), D=1-8t/q : lambda+ + lambda- = 1,  Q(y0)=4 w2/q, and (b^2 - a c) = r^2 m/(2q)
class QE:
    def __init__(s,a,b,D): s.a=F(a); s.b=F(b); s.D=D
    def __add__(s,o): o=s._c(o); return QE(s.a+o.a,s.b+o.b,s.D)
    __radd__=__add__
    def __sub__(s,o): o=s._c(o); return QE(s.a-o.a,s.b-o.b,s.D)
    def __rsub__(s,o): return s._c(o)-s
    def __mul__(s,o): o=s._c(o); return QE(s.a*o.a+s.b*o.b*s.D, s.a*o.b+s.b*o.a, s.D)
    __rmul__=__mul__
    def inv(s): n=s.a*s.a-s.b*s.b*s.D; return QE(s.a/n,-s.b/n,s.D)
    def __truediv__(s,o): o=s._c(o); return s*o.inv()
    def _c(s,o): return o if isinstance(o,QE) else QE(o,0,s.D)
    def iszero(s): return s.a==0 and s.b==0
ok2=True
for t in (F(1,10),F(-1,7),F(3,10),F(-2,5),F(1,9),F(-1,3),F(2,5),F(-1,11)):
    q=1+4*t*t; r=1-2*t; m=1+2*t; w2=1-2*t+4*t*t; D=1-8*t/q
    s0=QE(0,1,D); yp=(r/2)*(s0-1); ym=(-r/2)*(s0+1); y0=QE(r*r/q,0,D)
    lp=2*w2*yp/(q*(yp-ym)*(y0-yp)); lm=2*w2*ym/(q*(ym-yp)*(y0-ym))
    Qy0=4*t+(2*q/r)*y0+(2*q/(r*r))*y0*y0
    c1=(lp+lm-1).iszero(); c2=(Qy0-QE(4*w2/q,0,D)).iszero()
    # roots check: yp+ym=-r, yp*ym=2 t r^2/q
    c3=((yp+ym)+r).iszero() and (yp*ym-QE(2*t*r*r/q,0,D)).iszero()
    # closed formula lambda+ = w2 (s0-1)/(s0 (r^2 - (r q/2)(s0-1)))
    c4=(lp-(w2*(s0-1))/(s0*(r*r-(r*q/2)*(s0-1)))).iszero()
    if not(c1 and c2 and c3 and c4): ok2=False
    if D<0:
        a=-r/2; b2=-D*r*r/4; c=r*r/q+r/2
        if b2-a*c!=r*r*m/(2*q): ok2=False
print("S2 exact identities in Q(sqrt D): lambda_+ + lambda_- = 1, Q(y0)=4 w2/q, root relations, closed form of lambda_+, b^2-ac = r^2 m/(2q):",ok2)
# (3) exact: for t in (0,1/2) rational, n<=12: G1,G0 nonzero, exact degrees, binary resultant != 0, and G_i(1,0),G_i(0,1) not both zero
ok3=True; cnt3=0; sgn=set()
ts=[F(k,40) for k in range(1,20)]+[F(13,100),F(67,500),F(1,1000),F(499,1000)]
for t in ts:
    for n in range(1,13):
        (d1,G1),(d0,G0)=inf_forms(n,t)
        res=binres(G1,d1,G0,d0); cnt3+=1
        if not G1 or not G0 or res==0: ok3=False; print("FAIL",t,n)
        a1=G1.get((d1,0),0); a0=G0.get((d0,0),0); b1=G1.get((0,d1),0); b0=G0.get((0,d0),0)
        if (a1==0 and a0==0) or (b1==0 and b0==0): ok3=False; print("FAIL inf pts",t,n)
print(f"S3 exact: for 23 rational t in (0,1/2) and n=1..12 ({cnt3} cases): leading forms nonzero of exact degrees (d1,d0), binary resultant != 0, points [1:0],[0:1] not common roots:",ok3)
# (4) negative t: Res != 0 at rational t not in Lambda_n (sample) 
neg=[F(-1,10)+F(1,997),F(-3,10),F(-2,5),F(-1,50),F(-11,40)]
ok4=all(binres(inf_forms(n,t)[0][1],inf_forms(n,t)[0][0],inf_forms(n,t)[1][1],inf_forms(n,t)[1][0])!=0 for t in neg for n in range(1,9))
print("S4 exact: Res != 0 at 5 generic negative rational t (n<=8):",ok4)
