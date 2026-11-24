import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from fractions import Fraction as F
import exact as ex, pickle
def reduced_coeffs(n, t):
    R1, R0, E = ex.constraints_qp(n, t); res=[]
    for R in (R1, R0):
        mx=min(k[0] for k in R.d); mg=min(k[2] for k in R.d)
        res.append({((k[0]-mx)//2,k[1],(k[2]-mg)//2):v for k,v in R.d.items()})
    return res
def ABCD(t):
    c1,c0=reduced_coeffs(1,t)
    lin=lambda cf: dict(const=cf.get((0,0,0),0),u=cf.get((1,0,0),0),v=cf.get((0,0,1),0),e=cf.get((0,1,0),0))
    L1,L0=lin(c1),lin(c0); det=L1['u']*L0['v']-L1['v']*L0['u']
    sol=lambda r1,r0:((r1*L0['v']-L1['v']*r0)/det,(L1['u']*r0-r1*L0['u'])/det)
    u0,v0=sol(-L1['const'],-L0['const']); ue,ve=sol(-L1['e'],-L0['e'])
    return u0,-ue,v0,-ve
# polynomial utils (lists low->high)
def trim(p):
    while p and p[-1]==0: p=p[:-1]
    return p
def pmul(a,b):
    r=[F(0)]*(len(a)+len(b)-1) if a and b else []
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return trim(r)
def pdivmod(a,b):
    a=a[:]; q=[F(0)]*max(len(a)-len(b)+1,1)
    while len(a)>=len(b) and a:
        c=a[-1]/b[-1]; k=len(a)-len(b); q[k]=c
        for i,y in enumerate(b): a[i+k]-=c*y
        a=trim(a)
    return trim(q),a
def pgcd(a,b):
    while b: a,b=b,pdivmod(a,b)[1]
    return [c/a[-1] for c in a]
def peval(p,x):
    r=F(0)
    for c in reversed(p): r=r*x+c
    return r
def nullvec(M,ncol):
    M=[r[:] for r in M]; piv=[]; r=0
    for c in range(ncol):
        p=next((i for i in range(r,len(M)) if M[i][c]!=0),None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]; inv=1/M[r][c]; M[r]=[x*inv for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]; M[i]=[x-f*y for x,y in zip(M[i],M[r])]
        piv.append(c); r+=1
    free=[c for c in range(ncol) if c not in piv]
    if not free: return None
    fc=free[0]; v=[F(0)]*ncol; v[fc]=F(1)
    for i,c in enumerate(piv): v[c]=-M[i][fc]
    return v
def ratfit(vals,pts,maxd=45):
    for d in range(1,maxd):
        n=2*d+2
        if n+4>len(pts): return None
        rows=[[pts[i]**k for k in range(d+1)]+[-vals[i]*pts[i]**k for k in range(d+1)] for i in range(n)]
        v=nullvec(rows,2*d+2)
        if v is None: continue
        P=trim(v[:d+1]); Q=trim(v[d+1:])
        if not Q: continue
        if all(peval(P,pts[i])==vals[i]*peval(Q,pts[i]) for i in range(n,n+4)):
            g=pgcd(P,Q) if P else Q
            P=pdivmod(P,g)[0]; Q=pdivmod(Q,g)[0]
            lc=Q[-1]; return [c/lc for c in P],[c/lc for c in Q]
    return None
if __name__=="__main__":
    pts=[F(k,29) for k in range(-13,14) if k!=0]+[F(k,31)+F(1,1000) for k in range(-15,16)]+[F(3,5)+F(k,50) for k in range(0,12)]
    pts=[p for p in pts if abs(p)!=F(1,2) and p!=0]
    data=[ABCD(p) for p in pts]; print("evaluated",len(pts),"points")
    res={}
    for i,name in enumerate("ABCD"):
        r=ratfit([d[i] for d in data],pts)
        res[name]=r
        print(name, None if r is None else f"deg num {len(r[0])-1}, deg den {len(r[1])-1}")
    pickle.dump(res,open('ABCD_rat.pkl','wb'))
