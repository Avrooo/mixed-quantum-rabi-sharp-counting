"""Prop. 14(ii) sign logic: for t<0, t != t_p, no (eps) makes Delta^2=A-B*eps and g1^2=C-D*eps both positive,
and sgn(AD-BC)=sgn(Q*) as in Eq. (ADBC).  Also the corrected implication: B>0>D and C/D<A/B  =>  AD-BC<0."""
from fractions import Fraction as F
import ratrec
def feasible(A,B,C,D):
    lo,hi=None,None            # need A-B e>0 and C-D e>0
    for (a,b) in ((A,B),(C,D)):
        if b>0: hi=a/b if hi is None else min(hi,a/b)
        elif b<0: lo=a/b if lo is None else max(lo,a/b)
        elif a<=0: return False
    return lo is None or hi is None or lo<hi
def test_no_negative_t_solution():
    tp=-0.168601502
    for k in range(-499,0,7):
        t=F(k,1000)
        if abs(float(t)-tp)<2e-3: continue
        A,B,C,D=ratrec.ABCD(t)
        assert not feasible(A,B,C,D),t
def test_sign_of_ADBC_matches_Qstar():
    for k in range(-499,0,11):
        t=F(k,1000)
        if abs(float(t)-(-0.168601502))<2e-3: continue
        A,B,C,D=ratrec.ABCD(t); Qs=1+6*t+4*t**2+24*t**3+16*t**4 if False else None
        # sign relation: sgn(B)=varsigma, sgn(D)=-varsigma, sgn(AD-BC)=varsigma
        s=1 if B>0 else -1
        assert (D>0)==(s<0) and ((A*D-B*C)>0)==(s>0),t
def test_corrected_implication():
    for (A,B,C,D) in ((F(1),F(2),F(-3),F(-1)),(F(5),F(1),F(-7),F(-2)),(F(1,3),F(3),F(-1),F(-4))):
        assert B>0>D
        if C/D<A/B: assert A*D-B*C<0
