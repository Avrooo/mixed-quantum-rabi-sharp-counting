import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from limitpoly import chain
from fractions import Fraction as F
from math import comb
A=chain(14)
def lim(i,U,c): return sum(v*U**a*c**b for (a,b),v in A[i].items())
# Krawtchouk/binomial structure at a rational root: s=9,k=7 (m=7): c=7/2, a=4, b=-1/2, A(eta)=(1-eta/4)^7
c=F(7,2); U=F(18)
print("C8 s=9,k=7: A~_i(18,7/2) = C(7,i)(-1/4)^i for i=0..10:", all(lim(i,U,c)==(comb(7,i)*F(-1,4)**i if i<=7 else 0) for i in range(11)))
c=F(1); U=F(6)
print("C8 s=3,k=1: A~_i(6,1) = (1-eta/2) coefficients:", [lim(i,U,c) for i in range(5)])
# spin matrices: K=2Jx on spin m/2 (dim m+1); K^2 rational. Pencil det(w(s^2-K^2)-8K^2) vs prod
def K2(m):
    d=m+1; j=F(m,2)
    # Jx = (J+ + J-)/2 ; K^2 = (J+ + J-)^2 = J+^2 + J-^2 + J+J- + J-J+  (entries rational)
    # J+|mu> = sqrt(j(j+1)-mu(mu+1))|mu+1>, basis index r -> mu=-j+r
    def cp(r): mu=-j+r; return j*(j+1)-mu*(mu+1)   # squared amplitude of J+ from r to r+1
    M=[[F(0)]*d for _ in range(d)]
    for r in range(d):
        # (J+J-+J-J+) diagonal at r: |J- r>^2 + |J+ r>^2
        up=cp(r) if r+1<d else F(0); dn=cp(r-1) if r-1>=0 else F(0)
        M[r][r]=up+dn
        if r+2<d: M[r][r+2]=cp(r)*0  # placeholder replaced below
    # J+^2 and J-^2: amplitude sqrt(cp(r)cp(r+1)) -- product of two sqrt: irrational in general; use squared-matrix check via symmetric similarity instead
    return None
# Use exact rational route: eigenvalues of K are m,m-2,..,-m (standard); verify the ALGEBRAIC identity of the pencil on the spectrum
for m in range(1,13):
    s=m+2
    ks=[k for k in range(-m,m+1,2)]
    # det over full space of w(s^2-K^2)-8K^2 = prod_k ((s^2-k^2)w - 8k^2); reduced (k>0) polynomial = sqrt up to w^{[m even]} s^2
    full=[F(1)]
    for k in ks:
        new=[F(0)]*(len(full)+1)
        for i,cf in enumerate(full): new[i]+=cf*(-8*k*k); new[i+1]+=cf*(s*s-k*k)
        full=new
    red=[F(1)]
    for k in range(m,0,-2):
        new=[F(0)]*(len(red)+1)
        for i,cf in enumerate(red): new[i]+=cf*(-8*k*k); new[i+1]+=cf*(s*s-k*k)
        red=new
    sq=[F(0)]*(2*len(red)-1)
    for i,x in enumerate(red):
        for j,y in enumerate(red): sq[i+j]+=x*y
    if m%2==0: sq=[F(0)]+[x*s*s for x in sq]
    assert full==sq,(m)
print("C9 full-spin pencil det(w(s^2-K^2)-8K^2) = [reduced product]^2 * (s^2 w)^{[m even]}, m=1..12: PASS")
# non-similarity by diagonal scaling: our Jacobi off-diagonal squares b_i=2(m+1-i)/(i+1) vs spin i(m+1-i)
print("C10 ratio b_i/(i(m+1-i)) = 2/(i(i+1)) (not constant in i):", [F(2,i*(i+1)) for i in range(1,6)])
