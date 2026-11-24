"""Check of Prop. 12 (spin spectral identity): det(w(S^2-K^2)-8K^2) = [prod_{k>0}((S^2-k^2)w-8k^2)]^2 (S^2 w)^{[l even]},
K=2J_x on spin l/2, S=l+2.  Also: eigenvalue multiplicities of K^2 in the two parity sectors of exp(i pi J_z).
Exact rational arithmetic (fractions) at several rational w."""
import numpy as np
from fractions import Fraction as F
def Jx_Jz(l):
    j=F(l,2); m=[j-i for i in range(l+1)]
    return m
def run(lmax=12):
    ok=True; out=[]
    for l in range(1,lmax+1):
        S=l+2
        ks=[l-2*i for i in range(l+1)]           # spectrum of K=2J_x: 2m
        # parity sector of exp(i pi J_z): J_x eigenvectors |k> map to |-k> under that parity; sector spans (|k>+-|-k>)
        pos=[k for k in ks if k>0]; zero=[k for k in ks if k==0]
        sector_dims=(len(pos)+len(zero),len(pos))     # sector containing k=0 (if l even), the other
        for w in (F(1,3),F(2),F(7,5),F(11,2)):
            lhs=1
            for k in ks: lhs*= w*(S*S-k*k)-8*k*k
            rhs=1
            for k in pos: rhs*=((S*S-k*k)*w-8*k*k)**2
            if l%2==0: rhs*=S*S*w
            ok&=(lhs==rhs)
        out.append((l,len(pos),len(zero),sector_dims,ok))
    return ok,out
ok,out=run()
print("l  #k>0  #k=0  sector dims (with k=0 / other)  identity ok")
for r in out: print(*r)
print("ALL OK" if ok else "FAIL")
