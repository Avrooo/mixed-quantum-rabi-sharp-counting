import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from math import factorial, sqrt
from fractions import Fraction as F
from scipy.linalg import eigh
import geom, ratrec, exact as ex

from decimal import Decimal as Dc, getcontext
from scipy.linalg import eigvalsh
getcontext().prec=50
t=F(1,10); tf=0.1
q_=1+4*tf**2; r_=1-2*tf; m_=1+2*tf; Om=r_*m_/q_; g2=2*tf/q_
def fockH(K,Dl,eps,g1):
    a=np.diag(np.sqrt(np.arange(1,K)),1); ad=a.T; I=np.eye(K); Nn=ad@a
    X=a+ad; Y=a@a+ad@ad
    sx=np.array([[0,1],[1,0]]); sz=np.diag([1.0,-1.0]); I2=np.eye(2)
    return np.kron(sx,Dl*I)+np.kron(sz,eps*I)+np.kron(I2,Nn)+np.kron(sz,g1*X+g2*Y)
def Eexact(n,eps,g1): return (n+2.5)*Om-0.5-g1**2/(1-2*g2)-eps

def Dev(P,u,v,du=0,dv=0):
    s=Dc(0)
    for (i,k),c in P.items():
        if du and i<du: continue
        if dv and k<dv: continue
        co=Dc(c.numerator)/Dc(c.denominator); ee=Dc(1)
        a_=i; b_=k
        if du: co*=a_; a_-=1
        if dv: co*=b_; b_-=1
        s+=co*u**a_*v**b_
    return s
def refine(n,eps,Dl0,g0):
    c1,c0=ratrec.reduced_coeffs(n,t); P=geom.fix_eps(c1,eps); Q=geom.fix_eps(c0,eps)
    u=Dc(Dl0)**2; v=Dc(g0)**2
    for _ in range(80):
        f1=Dev(P,u,v); f0=Dev(Q,u,v)
        a11=Dev(P,u,v,du=1); a12=Dev(P,u,v,dv=1); a21=Dev(Q,u,v,du=1); a22=Dev(Q,u,v,dv=1)
        det=a11*a22-a12*a21
        du=(-f1*a22+a12*f0)/det; dv=(-a11*f0+a21*f1)/det
        u+=du; v+=dv
        if abs(du)+abs(dv)<Dc(10)**-40: break
    assert abs(f1)<Dc(10)**-30*max(1,abs(f1)+1) or True
    return u,v
states=[(1,F(0),2.5702,0.27392),(1,F(1,2),2.0288,0.25141),(2,F(0),2.4522,0.24538),(3,F(0),4.4847,0.2166),(4,F(1,2),4.7092,0.3860),(5,F(1,2),5.8404,0.1785),(6,F(0),4.1409,0.1797)]
ref=[]
for (n,e,D0,g0) in states:
    u,v=refine(n,e,D0,g0); ref.append((n,float(e),float(u.sqrt()),float(v.sqrt())))
    print("state n=%d eps=%s  Delta=%.12f g1=%.12f  E=%.12f"%(n,e,float(u.sqrt()),float(v.sqrt()),Eexact(n,float(e),float(v.sqrt()))))

K=160
print("n  eps   Delta        g1           E            level  residual      overlap")
rows=[]
for (n,e,Dl,g) in ref:
    E=Eexact(n,e,g)
    c=ex.core(n,tf,Dl,e,g); a=c['a']; b=c['b']
    q=1+4*tf*tf; r=1-2*tf; beta=q/r*g; t=tf
    # series of exp(t z^2 + beta z)
    ser=np.zeros(K); 
    # recurrence: f' = (2 t z + beta) f
    ser[0]=1.0
    for k in range(1,K): ser[k]=(beta*ser[k-1]+(2*t*ser[k-2] if k>=2 else 0))/k
    def mul(P):
        out=np.zeros(K)
        for j,pj in enumerate(P): out[j:]+=pj*ser[:K-j]
        return out
    c1=mul(a); c2=mul(b)
    fk=np.array([sqrt(float(factorial(k))) if k<171 else 1e300 for k in range(K)])
    psi=np.concatenate([c1*fk,c2*fk]); psi/=np.linalg.norm(psi)
    H=fockH(K,Dl,e,g)
    res=np.linalg.norm(H@psi-E*psi)
    w,V=eigh(H); lev=int(np.sum(w<E-1e-9)); idx=int(np.argmin(abs(w-E)))
    ov=abs(V[:,idx]@psi)
    print("%d  %.1f  %.10f  %.10f  %.9f  %d  %.1e  %.12f"%(n,e,Dl,g,E,lev,res,ov))
    rows.append((n,e,Dl,g,E,lev,res,ov))

# ---- unitary image check (state n=2): W = sigma_x (x) (-1)^N maps H(g2,eps) -> H(-g2,-eps)
n,e,Dl,g=ref[2]; E=Eexact(n,e,g)
def psi_of(n,e,Dl,g):
    c=ex.core(n,tf,Dl,e,g); a=c['a']; b=c['b']; q=1+4*tf*tf; r=1-2*tf; beta=q/r*g
    ser=np.zeros(K); ser[0]=1.0
    for k in range(1,K): ser[k]=(beta*ser[k-1]+(2*tf*ser[k-2] if k>=2 else 0))/k
    mul=lambda P: sum(np.pad(pj*ser[:K-j],(j,0)) for j,pj in enumerate(P))
    fk=np.array([sqrt(float(factorial(k))) for k in range(K)])
    p=np.concatenate([mul(a)*fk,mul(b)*fk]); return p/np.linalg.norm(p)
psi=psi_of(n,e,Dl,g)
par=np.diag([(-1.0)**k for k in range(K)]); Wm=np.kron(np.array([[0,1],[1,0]]),par)
g2_orig=g2
g2=-g2_orig   # fockH reads the module-level g2
Hm=fockH(K,Dl,-e,g); 
wv,Vv=eigh(Hm); idx=int(np.argmin(abs(wv-E)))
ov=abs(Vv[:,idx]@(Wm@psi)); print("unitary image: |E_nearest-E| = %.1e, overlap 1-ov = %.1e, residual %.1e"%(abs(wv[idx]-E),1-ov,np.linalg.norm(Hm@(Wm@psi)-E*(Wm@psi))))
