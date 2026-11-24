import os
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys,pickle,math; from fractions import Fraction as F
from cert import sturm_seq, var_at, sqfree, peval, trim
def lam_plus(t):
    q=1+4*t*t; r=1-2*t; m=1+2*t; w2=1-2*t+4*t*t
    s0=math.sqrt(1-8*t/q) if 1-8*t/q>=0 else None
    if s0 is None: return None
    return w2*(s0-1)/(s0*(r*r-(r*q/2)*(s0-1)))
def roots(P,lo,hi):
    seq=sturm_seq(P); out=[]
    def rec(a,b):
        c=var_at(seq,a)-var_at(seq,b)
        if c==0: return
        if c==1 and b-a<F(1,10**9): out.append(float((a+b)/2)); return
        mid=(a+b)/2
        if peval(seq[0],mid)==0: out.append(float(mid)); return
        rec(a,mid); rec(mid,b)
    rec(lo,hi); return sorted(out)
if __name__=="__main__":
    for t in (0.01,0.05,0.1,0.13,0.134,0.2,-0.01,-0.1,-0.3,-0.45): print(t,lam_plus(t))
    for n in (1,2,3,4):
        a,c1,c2,(b,P)=pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"data",f"resN_{n}.pkl"),"rb")); sf=sqfree(trim(P))
        rs=roots(sf,F(-1,2)+F(1,10**9),F(-1,10**9))
        print("n",n,"neg roots of Res numerator:",[round(x,6) for x in rs],"lambda+ there:",[round(lam_plus(x),6) for x in rs])
