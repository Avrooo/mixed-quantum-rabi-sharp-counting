import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from fractions import Fraction as F
import numpy as np
from ratrec import reduced_coeffs, trim, pdivmod, pgcd, peval, pmul
def fix_eps(c,e):
    # -> dict (i,k)->coef  for u^i v^k
    o={}
    for (i,j,k),v in c.items(): o[(i,k)]=o.get((i,k),0)+v*e**j
    return o
def det(M):
    M=[r[:] for r in M]; n=len(M); d=F(1)
    for c in range(n):
        p=next((r for r in range(c,n) if M[r][c]!=0),None)
        if p is None: return F(0)
        if p!=c: M[c],M[p]=M[p],M[c]; d=-d
        d*=M[c][c]
        for r in range(c+1,n):
            f=M[r][c]/M[c][c]
            if f: M[r]=[a-f*b for a,b in zip(M[r],M[c])]
    return d
def upoly(P,v,du):
    cs=[F(0)]*(du+1)
    for (i,k),c in P.items(): cs[i]+=c*v**k
    return cs  # low->high in u
def sylv(a,b):
    m=len(a)-1; n=len(b)-1; N=m+n; M=[]
    ar=a[::-1]; br=b[::-1]
    for i in range(n): M.append([F(0)]*i+ar+[F(0)]*(N-m-1-i))
    for i in range(m): M.append([F(0)]*i+br+[F(0)]*(N-n-1-i))
    return det(M)
def resultant_v(P,Q):
    du1=max(i for i,k in P); du0=max(i for i,k in Q)
    dv=du0*max(k for i,k in P)+du1*max(k for i,k in Q)
    pts=[F(j,3)+F(1,7) for j in range(dv+1)]
    vals=[sylv(upoly(P,x,du1),upoly(Q,x,du0)) for x in pts]
    # Lagrange -> coefficients
    coef=[F(0)]*(dv+1)
    for i,xi in enumerate(pts):
        num=[F(1)]; den=F(1)
        for j,xj in enumerate(pts):
            if j!=i: num=pmul(num,[-xj,F(1)]); den*=xi-xj
        for k,c in enumerate(num): coef[k]+=vals[i]*c/den
    return trim(coef)
def evalQ(P,u,v): return sum(float(c)*u**i*v**k for (i,k),c in P.items())
def real_solutions(c1,c0,e):
    P=fix_eps(c1,e); Q=fix_eps(c0,e); R=resultant_v(P,Q)
    if not R: return None
    rc=np.array([float(x) for x in R]); rc=rc/np.abs(rc).max()
    sols=[]
    for v in np.roots(rc[::-1]):
        if abs(v.imag)>1e-6*max(1,abs(v)) or v.real<=0: continue
        v=v.real
        pu=lambda PP: np.array([sum(float(c)*v**k for (i,k),c in PP.items() if i==ii) for ii in range(max(i for i,k in PP)+1)])
        a1=pu(P)
        for u in np.roots(a1[::-1]):
            if abs(u.imag)>1e-6*max(1,abs(u)) or u.real<=0: continue
            u=u.real
            # Newton polish in (u,v)
            x=np.array([u,v])
            for _ in range(30):
                f=np.array([evalQ(P,*x),evalQ(Q,*x)]); h=1e-7
                J=np.array([[ (evalQ(P,x[0]+h,x[1])-f[0])/h,(evalQ(P,x[0],x[1]+h)-f[0])/h],[(evalQ(Q,x[0]+h,x[1])-f[1])/h,(evalQ(Q,x[0],x[1]+h)-f[1])/h]])
                try: dx=np.linalg.solve(J,-f)
                except Exception: break
                x=x+dx
                if np.abs(dx).max()<1e-13: break
            sc=lambda PP:max(abs(c) for c in map(float,PP.values()))
            if x[0]>0 and x[1]>0 and abs(evalQ(P,*x))<1e-9*sc(P)*max(1,x.max())**3 and abs(evalQ(Q,*x))<1e-9*sc(Q)*max(1,x.max())**3:
                if all(abs(x[0]-s_[0])>1e-6*max(1,x[0]) or abs(x[1]-s_[1])>1e-6 for s_ in sols): sols.append((float(x[0]),float(x[1])))
    return sols
if __name__=="__main__":
    import sys
    t=F(1,10)
    for n in range(1,7):
        c1,c0=reduced_coeffs(n,t)
        d1=max(i+j+k for (i,j,k) in c1); d0=max(i+j+k for (i,j,k) in c0)
        print(f"n={n}: deg Q1={d1}, deg Q0={d0}, Bezout bound={d1*d0}")
        prev=None; row=[]
        for k in range(-60,61,2):
            e=F(k,20); s=real_solutions(c1,c0,e); m=len(s) if s is not None else -1
            row.append((float(e),m))
        print("   (eps,count of real solutions with Delta^2>0,g1^2>0):",row[::3])
        print("   max count over scan:",max(m for _,m in row))
