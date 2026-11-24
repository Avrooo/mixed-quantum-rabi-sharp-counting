"""Exact audit of parity stripping, degrees, endpoint (1:0),(0:1) and projective common-root-freeness, from the FULL
constraint polynomials R1,R0 of exact.py (independent of the weight recursion).  n=1..12."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import exact as ex
from gforms import binres
ts_pos=[F(1,50),F(1,10),F(1,7),F(1,5),F(3,10),F(2,5),F(49,100)]
ts_neg=[F(-1,10),F(-1,7),F(-2,5)]
tot=ok=0; fails=[]
for n in range(1,13):
    d1=(n+1)//2; d0=n//2+1
    for t in ts_pos+ts_neg:
        R1,R0,_=ex.constraints_qp(n,t)
        res=[]
        for R,par,d,nm in ((R1,(n+1)%2,d1,'R1'),(R0,n%2,d0,'R0')):
            ks=list(R.d)
            c_parity = all(k[0]%2==1 and k[2]%2==par for k in ks)          # Delta odd, g1 parity (n+1)%2 / n%2
            c_min = min(k[0] for k in ks)==1 and min(k[2] for k in ks)==par   # exactly stripped: no extra u or v factor
            Q={((k[0]-1)//2,k[1],(k[2]-par)//2):v for k,v in R.d.items()}
            deg=max(a+b+c for (a,c,b) in [(k[0],k[1],k[2]) for k in Q])      # total degree in (u,eps,v)
            top={(a,b):v for (a,e,b),v in Q.items() if a+e+b==d and e==0}     # eps-free part of the top-degree layer
            maxuv=max(a+b for (a,e,b) in Q)                                   # degree in (u,v) at fixed eps (eps constant)
            ok_deg = (deg==d)
            res.append((c_parity,c_min,ok_deg,maxuv==d,top,d))
        (p1,m1,dg1,uv1,G1,_),(p0,m0,dg0,uv0,G0,_)=res
        good = p1 and m1 and p0 and m0 and dg1 and dg0 and uv1 and uv0 and bool(G1) and bool(G0)
        # endpoint coefficients: G(1:0)=coef u^d, G(0:1)=coef v^d
        ep = {('G1(1:0)',G1.get((d1,0),0)!=0) ,('G1(0:1)',G1.get((0,d1),0)!=0),('G0(1:0)',G0.get((d0,0),0)!=0),('G0(0:1)',G0.get((0,d0),0)!=0)}
        both_inf = not ((G1.get((d1,0),0)==0) and (G0.get((d0,0),0)==0))   # (1:0) not a common root
        both_0 = not ((G1.get((0,d1),0)==0) and (G0.get((0,d0),0)==0))      # (0:1) not a common root
        R=binres(G1,d1,G0,d0)
        tot+=1
        cond = good and both_inf and both_0 and (R!=0 if t>0 else True)
        if cond: ok+=1
        else: fails.append((n,str(t),dict(parity=p1 and p0,strip=m1 and m0,deg=dg1 and dg0,uvdeg=uv1 and uv0,Gnz=bool(G1) and bool(G0),inf=both_inf,zero=both_0,Res_nonzero=R!=0)))
    print("n=%2d done  (cumulative %d/%d)"%(n,ok,tot),flush=True)
print("RESULT: %d/%d (n,t) pairs pass: parity, exact stripping, exact degrees d1,d0, nonzero top forms, (1:0),(0:1) not common roots, Res!=0 for t>0"%(ok,tot))
for f in fails: print("FAIL",f)
