"""Exact univariate rational functions over Q (lists low->high of Fractions); used for symbolic substitution."""
from fractions import Fraction as F
def trim(p):
    p=list(p)
    while p and p[-1]==0: p.pop()
    return p
def padd(a,b):
    n=max(len(a),len(b)); return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)])
def pneg(a): return [-c for c in a]
def pmul(a,b):
    if not a or not b: return []
    r=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return trim(r)
def pdivmod(a,b):
    a=trim(a); q=[F(0)]*max(len(a)-len(b)+1,1)
    while len(a)>=len(b) and a:
        c=a[-1]/b[-1]; k=len(a)-len(b); q[k]=c
        for i,y in enumerate(b): a[i+k]-=c*y
        a=trim(a)
    return trim(q),a
def pgcd(a,b):
    a,b=trim(a),trim(b)
    while b: a,b=b,pdivmod(a,b)[1]
    return a
class RF:
    def __init__(s,num,den=None,_n=True):
        num=trim([F(c) for c in (num if isinstance(num,(list,tuple)) else [num])])
        den=trim([F(c) for c in (den if den is not None else [1])])
        assert den
        if _n:
            if not num: den=[F(1)]
            else:
                g=pgcd(num,den)
                if len(g)>1: num=pdivmod(num,g)[0]; den=pdivmod(den,g)[0]
                c=den[-1]; num=[x/c for x in num]; den=[x/c for x in den]
        s.n,s.d=num,den
    @staticmethod
    def c(o): return o if isinstance(o,RF) else RF([o])
    def __add__(s,o):
        if not isinstance(o,(RF,int,F)): return NotImplemented
        o=RF.c(o); return RF(padd(pmul(s.n,o.d),pmul(o.n,s.d)),pmul(s.d,o.d))
    __radd__=__add__
    def __neg__(s): return RF(pneg(s.n),s.d,False)
    def __sub__(s,o):
        if not isinstance(o,(RF,int,F)): return NotImplemented
        return s+(-RF.c(o))
    def __rsub__(s,o): return RF.c(o)+(-s)
    def __mul__(s,o):
        if not isinstance(o,(RF,int,F)): return NotImplemented
        o=RF.c(o); return RF(pmul(s.n,o.n),pmul(s.d,o.d))
    __rmul__=__mul__
    def __truediv__(s,o):
        if not isinstance(o,(RF,int,F)): return NotImplemented
        o=RF.c(o); return RF(pmul(s.n,o.d),pmul(s.d,o.n))
    def __rtruediv__(s,o): return RF.c(o)/s
    def __eq__(s,o):
        o=RF.c(o) if isinstance(o,(RF,int,F)) else None
        return o is not None and s.n==o.n and s.d==o.d
    def __ne__(s,o): return not s.__eq__(o)
    def __hash__(s): return hash((tuple(s.n),tuple(s.d)))
    def __pow__(s,k):
        r=RF(1)
        for _ in range(k): r=r*s
        return r
    def __call__(s,x):
        f=lambda p: sum((c*x**i for i,c in enumerate(p)),F(0)); return f(s.n)/f(s.d)
T=RF([0,1])
