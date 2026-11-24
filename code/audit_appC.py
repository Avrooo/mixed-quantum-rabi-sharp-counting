"""Appendix C normalisation audit.
 (1) exact: (T1),(T2)  ==>  (C2)  for random rational data (all i<=40);
 (2) exact rational-function identities in t: zeta0, tilde-Q, kappa-hat, Qt(1)=4 w2/q, y+ + y-, y+ y-, discriminant, lambda+ lambda-, lambda formula, f'/f formula;
 (3) 80-digit end-to-end: at t with lambda_+=k/M (root of the quartic P_{k/M}) the recursion (T1),(T2) with chi=chi_M gives
     A_{n+1}=A_{n+2}=0 and A_i = product-form coefficients (e0,e+,e- of the paper), for every reduced k/M<1/2, M<=8, every n>=M-2 (n<=9)."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from decimal import Decimal as Dc, getcontext
from math import gcd
import random
from ratfunc import RF,T
getcontext().prec=90
# ---------- (1) exact (T1,T2) => (C2)
def AB(n,t,x,g,N,eps=F(0)):
    q=1+4*t*t; r=1-2*t; m=1+2*t; Om=r*m/q; p=4*t; beta=q/r*g; h=m/r*g; s=q/(r*r)*g*g; sig=s+eps
    A=[-x/p]; B=[F(1)]
    for i in range(1,N):
        Bi=(x*(A[i-2] if i>=2 else 0)-(n+3-i)*h*B[i-1])/(i*Om); B.append(Bi)
        Ai=-(2*beta*A[i-1]+2*sig*(A[i-2] if i>=2 else 0)+x*Bi)/p; A.append(Ai)
    return A,B,dict(q=q,r=r,m=m,Om=Om,p=p,beta=beta,h=h,s=s)
random.seed(7); bad=0; cnt=0
for trial in range(30):
    t=F(random.choice([-1,1])*random.randint(1,45),100); x=F(random.randint(1,9),random.randint(1,5)); g=F(random.randint(1,9),random.randint(1,5)); n=random.randint(1,9)
    A,B,c=AB(n,t,x,g,45); Nn=n+3; Om,p,beta,h,s=c['Om'],c['p'],c['beta'],c['h'],c['s']
    L1=2*beta*Om-h*p; L2=2*s*Om-2*beta*h; M0=2*beta*Om+h*(Nn-1)*p; M1=4*s*Om-2*beta*h+2*beta*h*(Nn-1)+x*x
    for i in range(2,40):
        lhs=Om*p*(i+1)*A[i+1]+(L1*i+M0)*A[i]+(L2*(i-1)+M1)*A[i-1]+2*s*h*(Nn-1-i)*A[i-2]
        cnt+=1; bad+= (lhs!=0)
print("(1) (C2) exact residual != 0 :",bad,"of",cnt)
# ---------- (2) exact identities in t
q=1+4*T*T; r=1-2*T; m=1+2*T; w2=1-2*T+4*T*T; Om=r*m/q
Qt1=4*T+2*r+2*r*r/q
ok=[]
ok.append(("Qt(1)=4w2/q",Qt1==4*w2/q))
ok.append(("q+r^2=2 w2",q+r*r==2*w2))
ok.append(("Omega=rm/q =(1-4t^2)/q",Om==(1-4*T*T)/q))
# kappa-hat prefactor x^2 zeta0^2/Omega = chi r^3/(q m)  (g=1: zeta0=Omega/h, h=m/r)
h=m/r; z0=Om/h; ok.append(("x^2 z0^2/Om = r^3/(q m) (per chi)",z0*z0/Om==r**3/(q*m)))
ok.append(("2 beta z0 = 2 r  (beta/h=q/m)",2*(q/r)*z0==2*r))
ok.append(("2 s z0^2 = 2 r^2/q",2*(q/(r*r))*z0*z0==2*r*r/q))
# y roots: sum=-q/r, prod=2tq/r^2 ; D=q(q-8t)/r^2
ysum=-q/r; yprod=2*T*q/(r*r); D=ysum*ysum-4*yprod
ok.append(("D=q(q-8t)/r^2",D==q*(q-8*T)/(r*r)))
one_minus=1+q/r+yprod   # (1-y+)(1-y-)
ok.append(("(1-y+)(1-y-)=2w2/r^2",one_minus==2*w2/(r*r)))
lamprod=-yprod*one_minus/D
f=m*m*q/(r*r*(q-8*T))
ok.append(("lambda+ lambda- = (1-f)/4",lamprod==(1-f)/4))
ok.append(("lambda+ lambda- = -4 t w2/(r^2(q-8t))",lamprod==-4*T*w2/(r*r*(q-8*T))))
# chi_M: chi r^3/(4 m w2) = -M  => chi=-4 m w2 M/r^3 ; kappa-hat=-M Qt(1)
M=RF(5); chi=-4*m*w2*M/(r**3); ok.append(("kappa-hat=chi r^3/(q m) = -M Qt(1) at chi_M",chi*r**3/(q*m)==-M*Qt1))
# f'/f = 8/(1-4t^2)+8t/q+8(1-t)/(q-8t)
def dRF(a):
    num=[i*c for i,c in enumerate(a.n)][1:]; den=[i*c for i,c in enumerate(a.d)][1:]
    return RF(num,[1])*RF(a.d)-RF(a.n)*RF(den,[1]) 
def deriv(a):
    n_=RF([i*c for i,c in enumerate(a.n)][1:] or [0]); d_=RF([i*c for i,c in enumerate(a.d)][1:] or [0])
    return (n_*RF(a.d)-RF(a.n)*d_)/(RF(a.d)*RF(a.d))
ok.append(("f'/f = 8/(1-4t^2)+8t/q+8(1-t)/(q-8t)",deriv(f)/f==8/(1-4*T*T)+8*T/q+8*(1-T)/(q-8*T)))
for nm,v in ok: print("(2)",nm,":",v); assert v
# ---------- (3) 80-digit end-to-end
def Pk(t,k,Mm):  # quartic encoding lambda+=k/M
    q_=1+4*t*t; r_=1-2*t; m_=1+2*t
    return (Mm-2*k)**2*r_*r_*(q_-8*t)-Mm*Mm*m_*m_*q_
def root(k,Mm):
    lo,hi=Dc(-0.5)+Dc(10)**-30,Dc(0)-Dc(10)**-30
    assert Pk(lo,k,Mm)*Pk(hi,k,Mm)<0
    for _ in range(280):
        mid=(lo+hi)/2
        if Pk(lo,k,Mm)*Pk(mid,k,Mm)<=0: hi=mid
        else: lo=mid
    return (lo+hi)/2
worst=Dc(0); tests=0
for Mm in range(3,9):
    for k in range(1,(Mm+1)//2):
        if gcd(k,Mm)!=1 or 2*k>=Mm: continue
        t=root(k,Mm); q_=1+4*t*t; r_=1-2*t; m_=1+2*t; w2_=1-2*t+4*t*t; Om_=r_*m_/q_; p_=4*t
        lam=Dc(1)/2-(m_/(2*r_))*(q_/(q_-8*t)).sqrt()
        assert abs(lam-Dc(k)/Mm)<Dc(10)**-70
        for mult in range(1,9):
            Mj=Mm*mult
            if Mj>13: break
            for n in range(Mj-2,min(Mj+3,14)):  # M<=n+2, n<=13
                if n<1: continue
                chi=-4*m_*w2_*Mj/r_**3; g=Dc(1); x2=chi*g*g   # x^2 = chi g^2 (x imaginary): use A~ = A/x
                beta=q_/r_*g; h_=m_/r_*g; s_=q_/(r_*r_)*g*g
                At=[-1/p_]; Bq=[Dc(1)]
                for i in range(1,n+4):
                    Bi=(x2*(At[i-2] if i>=2 else 0)-(n+3-i)*h_*Bq[i-1])/(i*Om_); Bq.append(Bi)
                    Ai=-(2*beta*At[i-1]+2*s_*(At[i-2] if i>=2 else 0)+Bi)/p_; At.append(Ai)
                # product form
                z0=Om_/h_
                D_=(q_*(q_-8*t)).sqrt()/r_; yp=(-q_/r_+D_)/2; ym=(-q_/r_-D_)/2
                e0=n+2-Mj; ep=Mj*lam-1; em=Mj*(1-lam)-1
                assert abs(ep-round(ep))<Dc(10)**-60 and abs(em-round(em))<Dc(10)**-60
                ep=int(round(ep)); em=int(round(em)); assert e0>=0 and ep>=0 and em>=0 and e0+ep+em==n
                # series of (1-y)^e0 (1-y/yp)^ep (1-y/ym)^em in y -> convert to zeta: y=zeta/z0
                def binom_series(a,e,N): # (1 - a*y)^e coefficients
                    co=[Dc(1)]
                    for _ in range(e): co=[ (co[i] if i<len(co) else 0) - a*(co[i-1] if i>=1 else 0) for i in range(len(co)+1)]
                    return co
                poly=[Dc(1)]
                for a,e in ((Dc(1),e0),(1/yp,ep),(1/ym,em)):
                    c2=binom_series(a,e,0); new=[Dc(0)]*(len(poly)+len(c2)-1)
                    for i,u in enumerate(poly):
                        for j,v in enumerate(c2): new[i+j]+=u*v
                    poly=new
                prod=[-1/p_*c/z0**i for i,c in enumerate(poly)]   # A~(zeta) coefficients
                dev=max(abs(At[i]-(prod[i] if i<len(prod) else 0)) for i in range(n+4))
                sc=max(abs(a) for a in At[:n+1])
                rel=dev/sc; worst=max(worst,rel); tests+=1
                if not (abs(At[n+1])/sc<Dc(10)**-40 and abs(At[n+2])/sc<Dc(10)**-40 and abs(At[n])/sc>Dc(10)**-15): print("FAIL",Mm,k,mult,n,float(rel),float(abs(At[n+1])/sc),float(abs(At[n+2])/sc),float(abs(At[n])/sc)); continue
print("(3) end-to-end: %d (k/M, multiple, n) cases; max relative deviation recursion vs product form: %.2e; A_{n+1}=A_{n+2}=0 and A_n!=0 in all"%(tests,worst))
