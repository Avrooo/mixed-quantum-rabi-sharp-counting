"""Check: at U=2s (s=m+2), c=c_j the limiting polynomials A_i of recurrence (L) equal
A_i = binom(m,i) K_i(j-1; p=j/s, m) (-1/b)^i,  K_i = 2F1(-i,-x;-N;1/p) (Krawtchouk), b<0 the negative root of eta^2-c*eta-2."""
from math import comb, factorial
from decimal import Decimal as D, getcontext
getcontext().prec=80
def sqrt(x): return D(x).sqrt()
def poch(a,k):
    r=D(1)
    for i in range(k): r*=D(a)+i
    return r
def kraw(i,x,p,N):
    return sum(poch(-i,k)*poch(-x,k)/(poch(-N,k)*factorial(k))*D(p)**(-k) for k in range(0,min(i,x)+1))
def A_rec(U,c,nmax):
    A=[D(1),-(c/2)]  # A_0=1, A_1=-(c/2)A_0 (A_{-1}=0)
    for i in range(2,nmax+1):
        A.append(-(c/2)*A[i-1]+(D(2*i-U)/D(4*i))*A[i-2])
    return A
worst=0; term=0
for m in range(1,13):
    s=m+2
    for j in range(1,m+2):
        c=sqrt(2)*(2*j-s)/sqrt(j*(s-j)); rho=sqrt(c*c+8); a=(c+rho)/2; b=(c-rho)/2
        assert abs(a/rho-D(j)/D(s))<D('1e-40')
        A=A_rec(2*s,c,m+4)
        term=max(term,max(abs(A[i]) for i in range(m+1,m+5)))      # must vanish for i>m
        for i in range(0,m+1):
            K=kraw(i,j-1,D(j)/D(s),m); pred=comb(m,i)*K*(-1/b)**i
            worst=max(worst,abs(A[i]-pred)/max(D(1),abs(A[i])))
print(f"max relative deviation recurrence vs Krawtchouk formula (m<=12, all j): {float(worst):.2e}")
print(f"max |A_i| for m<i<=m+4 (termination, expect ~0): {float(term):.2e}")
assert worst<D('1e-40') and term<D('1e-40')
print("KRAWTCHOUK IDENTIFICATION VERIFIED")
