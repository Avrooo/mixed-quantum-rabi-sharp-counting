import os
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys,pickle,math,csv; from fractions import Fraction as F
from cert import sqfree, trim, peval, sturm_seq, var_at
from lam import lam_plus, roots
from math import gcd
def farey(n):
    S=set()
    for M in range(2,n+3):
        for k in range(1,M):
            if 2*k<M and gcd(k,M)==1: S.add(F(k,M))
    return sorted(S)
rows=[]
print("n | #real roots of Res-numerator in (-1/2,0) | in (0,1/2) | predicted #reduced fractions k/M<1/2, M<=n+2 | lambda+ at roots == predicted set?")
for n in range(1,7):
    a,c1,c2,(b,P)=pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"data",f"resN_{n}.pkl"),"rb")); sf=sqfree(trim(P))
    rneg=roots(sf,F(-1,2)+F(1,10**9),F(-1,10**9)); rpos=roots(sf,F(1,10**9),F(1,2)-F(1,10**9))
    pred=farey(n)
    lam=[lam_plus(x) for x in rneg]
    ok=len(rneg)==len(pred) and all(abs(l-float(p))<1e-6 for l,p in zip(sorted(lam),pred))
    print(n,len(rneg),len(rpos),len(pred),ok,flush=True)
    for x,l in zip(rneg,lam):
        fr=min(pred,key=lambda p:abs(float(p)-l))
        rows.append(dict(n=n,t_root=f"{x:.12f}",g2_over_omega=f"{2*x/(1+4*x*x):.9f}",lambda_plus=f"{l:.9f}",k_over_M=str(fr),
            chi_u_over_v_for_smallest_n0=f"{-4*(1+2*x)*(1-2*x+4*x*x)*fr.denominator/(1-2*x)**3:.6f}",common_root_at_infinity="yes (exact Sturm root of Res numerator)"))
w=csv.DictWriter(open("EXACT_EXCEPTIONAL_CASES.csv","w",newline=""),fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
