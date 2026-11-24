# exact t=0 limit polynomials A~_i(U,c): A~_0=1, A~_{-1}=0, A~_{i+1} = -(c/2) A~_i + ((2(i+1)-U)/(4(i+1))) A~_{i-1}
from fractions import Fraction as F
def padd(a,b,s=1):
    r=dict(a)
    for k,v in b.items():
        r[k]=r.get(k,0)+s*v
        if r[k]==0: del r[k]
    return r
def pmulc(a,c): return {k:v*c for k,v in a.items()} if c!=0 else {}
def pmulm(a,mono,c=F(1)): return {(k[0]+mono[0],k[1]+mono[1]):v*c for k,v in a.items()}
def chain(N):
    A=[{(0,0):F(1)}]; prev={}
    for i in range(0,N):
        nxt=padd(pmulm(A[i],(0,1),F(-1,2)), padd(pmulc(prev,F(2*(i+1),4*(i+1))), pmulm(prev,(1,0),F(-1,4*(i+1)))))
        prev=A[i]; A.append(nxt)
    return A   # A[i] dict {(powU,powC):coef}
def strip(poly,par):
    """parity strip: poly = c^par * P(U,w); return dict {(powU,powW):coef}"""
    out={}
    for (i,j),v in poly.items():
        assert (j-par)%2==0 and j>=par,(i,j,par)
        out[(i,(j-par)//2)]=v
    return out
def peval2(P,U,w):
    return sum(v*U**i*w**j for (i,j),v in P.items())
def dU(P): return {(i-1,j):v*i for (i,j),v in P.items() if i>0}
def dW(P): return {(i,j-1):v*j for (i,j),v in P.items() if j>0}
def univ_c(poly,U0):
    d={}
    for (i,j),v in poly.items(): d[j]=d.get(j,0)+v*U0**i
    n=max(d) if d else 0
    return [d.get(j,F(0)) for j in range(n+1)]
