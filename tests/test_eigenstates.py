"""Analytic Gaussian-polynomial spinor vs Fock-space diagonalization for the Table V state n=1, eps=1/2 (t=1/10)."""
import numpy as np
from scipy.linalg import eigh
import exact as ex
tf=0.1; q=1+4*tf**2; r=1-2*tf; g2=2*tf/q; Om=r*(1+2*tf)/q
Dl,g,e,n=2.0288407037594594,0.25141065329047163,0.5,1
def fockH(K):
    a=np.diag(np.sqrt(np.arange(1,K)),1); ad=a.T; I=np.eye(K); N=ad@a
    sx=np.array([[0,1],[1,0.]]); sz=np.diag([1.,-1.])
    return np.kron(sx,Dl*I)+np.kron(sz,e*I)+np.kron(np.eye(2),N)+np.kron(sz,g*(a+ad)+g2*(a@a+ad@ad))
def spinor(K):
    c=ex.core(n,tf,Dl,e,g); beta=q/r*g; ser=np.zeros(K); ser[0]=1
    for k in range(1,K): ser[k]=(beta*ser[k-1]+(2*tf*ser[k-2] if k>=2 else 0))/k
    mul=lambda P: sum(np.pad(pj*ser[:K-j],(j,0)) for j,pj in enumerate(P))[:K]
    from scipy.special import gammaln
    fk=np.exp(0.5*gammaln(np.arange(K)+1)); p=np.concatenate([mul(c['a'])*fk,mul(c['b'])*fk]); return p/np.linalg.norm(p)
def test_energy_and_overlap():
    K=64; E=(n+2.5)*Om-0.5-g**2/(1-2*g2)-e; p=spinor(K); H=fockH(K)
    assert np.linalg.norm(H@p-E*p)<1e-12
    w,V=eigh(H); i=int(np.argmin(abs(w-E))); assert abs(w[i]-E)<1e-12 and abs(abs(V[:,i]@p)-1)<1e-12
