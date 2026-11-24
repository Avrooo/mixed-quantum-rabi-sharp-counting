import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from fractions import Fraction as F
from limit_sys import leading_forms
def inf_forms(n,t):
    """(u,v)-leading homogeneous parts G1,G0 (epsilon-free part of the top weighted forms) as dicts {(i,k):coef} u^i v^k"""
    out=[]
    for form in leading_forms(n,t):
        d=max(i+j+k for (i,j,k) in form)
        out.append((d,{(i,k):c for (i,j,k),c in form.items() if j==0 and i+k==d and c!=0}))
    return out
def binres(G1,d1,G0,d0):
    """resultant of binary forms G(u,v)=sum c_{i,k}u^i v^k (i+k=d) via Sylvester on coefficients of u^d.. v^d"""
    a=[G1.get((d1-k,k),F(0)) for k in range(d1+1)]   # coefficients in v (set u=1)
    b=[G0.get((d0-k,k),F(0)) for k in range(d0+1)]
    N=d1+d0
    if N==0: return F(1)
    rows=[]
    for i in range(d0): rows.append([F(0)]*i+a+[F(0)]*(N-i-d1-1))
    for i in range(d1): rows.append([F(0)]*i+b+[F(0)]*(N-i-d0-1))
    return det(rows)
def det(M):
    M=[r[:] for r in M]; n=len(M); d=F(1)
    for c in range(n):
        piv=None
        for r in range(c,n):
            if M[r][c]!=0: piv=r;break
        if piv is None: return F(0)
        if piv!=c: M[c],M[piv]=M[piv],M[c]; d=-d
        d*=M[c][c]
        for r in range(c+1,n):
            f=M[r][c]/M[c][c]
            if f!=0:
                for k in range(c,n): M[r][k]-=f*M[c][k]
    return d
if __name__=="__main__":
    for n in range(1,8):
        for t in (F(1,10),F(3,10),F(-1,10)):
            (d1,G1),(d0,G0)=inf_forms(n,t)
            print(n,t,"d1,d0",d1,d0,"nonzero G1,G0:",bool(G1),bool(G0),"Res",float(binres(G1,d1,G0,d0)) if G1 and G0 else "n/a")
