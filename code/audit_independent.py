"""Independent check of the attainment/counting claims using raw operators only.
For selected (n,t,eps): solve the constraint system (geom.real_solutions), build the analytic Gaussian-polynomial spinor in a
Fock basis, and test ||(H-E)psi|| directly with the Hamiltonian matrix (no band formulas) and against Fock diagonalization.
Exits nonzero if any check fails."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as F
from scipy.linalg import eigh
from scipy.special import gammaln
import ratrec, geom, exact as ex
fails=0
def run(n,t,eps,K=160,tol=1e-7):
    global fails
    tf=float(t); q=1+4*tf*tf; r=1-2*tf; g2=2*tf/q; Om=r*(1+2*tf)/q; B=((n+1)**2)//4
    c1,c0=ratrec.reduced_coeffs(n,t); sols=geom.real_solutions(c1,c0,eps) or []
    a=np.diag(np.sqrt(np.arange(1,K)),1); ad=a.T; I=np.eye(K); N=ad@a
    sx=np.array([[0,1],[1,0.]]); sz=np.diag([1.,-1.]); worst=0.0; wo=1.0; ok=len(sols)==B
    for (u,v) in sols:
        Dl,g=np.sqrt(u),np.sqrt(v); e=float(eps)
        H=np.kron(sx,Dl*I)+np.kron(sz,e*I)+np.kron(np.eye(2),N)+np.kron(sz,g*(a+ad)+g2*(a@a+ad@ad))
        E=(n+2.5)*Om-0.5-g**2/(1-2*g2)-e
        c=ex.core(n,tf,Dl,e,g); beta=q/r*g; ser=np.zeros(K); ser[0]=1
        for k in range(1,K): ser[k]=(beta*ser[k-1]+(2*tf*ser[k-2] if k>=2 else 0))/k
        mul=lambda P: sum(np.pad(pj*ser[:K-j],(j,0)) for j,pj in enumerate(P))[:K]
        fk=np.exp(0.5*gammaln(np.arange(K)+1)); p=np.concatenate([mul(c['a'])*fk,mul(c['b'])*fk]); p/=np.linalg.norm(p)
        res=np.linalg.norm(H@p-E*p)/max(1,abs(E)); w,V=eigh(H); i=int(np.argmin(abs(w-E))); ov=abs(V[:,i]@p)
        worst=max(worst,res); wo=min(wo,ov); ok&=(res<tol and abs(w[i]-E)<1e-6*max(1,abs(E)) and ov>1-1e-6)
    print("n=%d t=%s eps=%s: solutions %d (B_n=%d)  worst rel. residual %.1e  min overlap %.9f  %s"%(n,t,eps,len(sols),B,worst,wo,"PASS" if ok else "FAIL"))
    if not ok: fails+=1
for n in (2,3,4): run(n,F(1,20),F(-100))
for n in (2,3): run(n,F(1,50),F(-100))
run(1,F(1,10),F(1,2)); run(2,F(1,10),F(0))
sys.exit(1 if fails else 0)
