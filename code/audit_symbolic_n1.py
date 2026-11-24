"""Exact SYMBOLIC derivation of the n=1 closed forms (no interpolation): substitute t -> indeterminate in the exact
recursion of exact.py (rational-function arithmetic over Q), then compare with the formulas printed in the paper."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from ratfunc import RF,T,pmul,padd,pdivmod,pgcd,trim
import ratrec
A,B,C,D=ratrec.ABCD(T)
t=T; q=1+4*t*t; Qp=RF([1,6,4,24,16])
pa=RF([9,36,-24,-432,-96,-1728,-384,2304,2304])/(q*q*Qp)
pb=6*(1-4*t*t)*(1+2*t)**2/Qp
pc=t*(1-2*t)**3*(1+2*t)*(5-8*t+20*t*t)/(2*q*q*Qp)
pd=t*(1-2*t)**4/(q*Qp)
for nm,x,y in (("A",A,pa),("B",B,pb),("C",C,pc),("D",D,pd)):
    print(nm,"symbolic == printed :",x==y)
    assert x==y
dt=A*D-B*C
print("AD-BC numerator (low->high):",[str(c) for c in dt.n]); print("denominator:",[str(c) for c in dt.d])
# factor out t(1-2t)^4/(q^3 Q^2)
pref=t*(1-2*t)**4/(q**3*Qp*Qp)
W=dt/pref
print("W = (AD-BC)/[t(1-2t)^4/(q^3 Q^2)] is polynomial?",W.d==[1], " deg",len(W.n)-1)
print("W coeffs low->high:",[str(c) for c in W.n])
# is W divisible by (1-2t)? by m=(1+2t)? factor repeated roots
def divmult(p,root_poly):
    k=0
    while True:
        qq,r=pdivmod(p,root_poly)
        if r: return k,p
        p=qq;k+=1
for nm,fac in (("(1-2t)",[F(1),F(-2)]),("(1+2t)",[F(1),F(2)]),("t",[F(0),F(1)]),("q",[F(1),F(0),F(4)])):
    print("mult of",nm,"in W:",divmult(W.n,fac)[0])
# strict sign of W on (0,1/2): evaluate, plus exact Sturm count of roots in (0,1/2)

def sturm(P):
    P=trim(P); seq=[P,trim([i*c for i,c in enumerate(P)][1:])]
    while len(seq[-1])>1:
        r=pdivmod(seq[-2],seq[-1])[1]
        if not r: break
        seq.append([-c for c in r])
    return seq
def nroots(P,lo,hi):
    g=pgcd(P,trim([i*c for i,c in enumerate(P)][1:]))
    if len(g)>1: P=pdivmod(P,g)[0]
    seq=sturm(P)
    def var(x):
        v=[sum((c*x**i for i,c in enumerate(p)),F(0)) for p in seq]; v=[a for a in v if a!=0]
        return sum(1 for a,b in zip(v,v[1:]) if a*b<0)
    return var(lo)-var(hi)
Wn=W.n
print("distinct real roots of W in (0,1/2):",nroots(Wn,F(0),F(1,2)), "; in (-1/2,0):",nroots(Wn,F(-1,2),F(0)))
print("W(0),W(1/4),W(1/2) =",[float(sum(c*x**i for i,c in enumerate(Wn))) for x in (F(0),F(1,4),F(1,2))])
print("distinct real roots of Q in (-1/2,1/2):",nroots([1,6,4,24,16],F(-1,2),F(1,2)))
print("sign of W on (0,1/2): W(0)=%s -> negative iff W<0"%Wn[0])
