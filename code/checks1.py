import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from limitpoly import *
import itertools
N=42; A=chain(N+8)
# C1: ODE identity (eta^2 - c eta - 2) A' = ((U/2-2) eta + c) A  coefficientwise, i up to 30
bad=0
for i in range(0,34):
    lhs=padd(padd(pmulc(A[i-1],F(i-1)) if i>=1 else {}, pmulm(A[i],(0,1),F(-i))), pmulc(A[i+1],F(-2*(i+1))))
    rhs=padd(padd(pmulm(A[i-1],(1,0),F(1,2)) if i>=1 else {}, pmulc(A[i-1],F(-2)) if i>=1 else {}), pmulm(A[i],(0,1)))
    if padd(lhs,rhs,-1): bad+=1
print("C1 ODE identity for coefficients i=0..33:", "PASS" if bad==0 else f"FAIL {bad}")
# helper: univariate in c at fixed U
def umul(a,b):
    r=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def utrim(p):
    while p and p[-1]==0: p=p[:-1]
    return p
def udivmod(a,b):
    a=utrim(a[:]); b=utrim(b); q=[F(0)]*max(len(a)-len(b)+1,1)
    while len(a)>=len(b) and a:
        c=a[-1]/b[-1]; k=len(a)-len(b); q[k]=c
        for i,y in enumerate(b): a[i+k]-=c*y
        a=utrim(a)
    return utrim(q),a
def ugcd(a,b):
    a=utrim(a); b=utrim(b)
    while b: a,b=b,udivmod(a,b)[1]
    return [x/a[-1] for x in a]
# C2: Atilde_{m+1}(2s,c) stripped = const * prod_k ((s^2-k^2) w - 8 k^2)
ok2=True; ok3=True; ok3b=True
for m in range(1,41):
    s=m+2; U0=F(2*s); par=(m+1)%2
    P=strip({k:v for k,v in {(0,j):c for j,c in enumerate(univ_c(A[m+1],U0)) if c!=0}.items()},par)
    Pw=[P.get((0,j),F(0)) for j in range(max(j for (_,j) in P)+1)]
    prod=[F(1)]
    for k in range(m,0,-2): prod=umul(prod,[F(-8*k*k),F(s*s-k*k)])
    prod=utrim(prod); Pw=utrim(Pw)
    if len(Pw)!=len(prod) or any(Pw[i]*prod[-1]!=prod[i]*Pw[-1] for i in range(len(prod))): ok2=False; print("C2 fail",m)
    # C3: A_i(2s,c) divisible by A_{m+1}(2s,c) for i=m+2..m+6 ; gcd(A_m,A_{m+1})=1
    a1=utrim(univ_c(A[m+1],U0))
    for i in range(m+2,m+7):
        if udivmod(utrim(univ_c(A[i],U0)),a1)[1]: ok3=False
    if len(ugcd(utrim(univ_c(A[m],U0)),a1))!=1: ok3b=False
print("C2 stripped A~_{m+1}(2(m+2),w) = const*prod_k((s^2-k^2)w-8k^2), m=1..40:", "PASS" if ok2 else "FAIL")
print("C3 A~_i(2s,c) divisible by A~_{m+1} for i=m+2..m+6 (m<=40):", "PASS" if ok3 else "FAIL", "| gcd(A~_m,A~_{m+1})=1:", "PASS" if ok3b else "FAIL")
# C5: classification of common zeros of (A_{n+1},A_n): at U0 not in {2,4,..,2(n+1)} gcd=1 ; at U0=2s (2<=s<=n+1) gcd degree = s-1
ok5=True
for n in range(1,31):
    g=ugcd(utrim(univ_c(A[n+1],F(1,3))),utrim(univ_c(A[n],F(1,3))))
    if len(g)!=1: ok5=False
    for s in range(2,n+2):
        g=ugcd(utrim(univ_c(A[n+1],F(2*s))),utrim(univ_c(A[n],F(2*s))))
        if len(g)-1!=s-1: ok5=False; print("C5 fail",n,s,len(g)-1)
print("C5 common zeros of (A~_{n+1},A~_n): none off U in {4,..,2(n+1)}; at U=2s exactly deg s-1 (n<=30):","PASS" if ok5 else "FAIL")
# C4: Jacobian of stripped system (P1,P0) at every indexed point, 1<=m<=n<=24
ok4=True; cnt=0; zero_found=[]
for n in range(1,25):
    P1=strip(A[n+1],(n+1)%2)
    Phi0=padd(pmulm(A[n],(1,0),F(1,n+2)),pmulc(A[n],F(-2)))   # (U-2n-4)/(n+2) * A_n
    P0=strip(Phi0,n%2)
    for m in range(1,n+1):
        s=m+2
        for k in range(m,0,-2):
            U0=F(2*s); w0=F(8*k*k,s*s-k*k)
            v1=peval2(P1,U0,w0); v0=peval2(P0,U0,w0)
            if v1!=0 or v0!=0: ok4=False; print("not a root",n,m,k,v1,v0)
            J=peval2(dU(P1),U0,w0)*peval2(dW(P0),U0,w0)-peval2(dW(P1),U0,w0)*peval2(dU(P0),U0,w0)
            cnt+=1
            if J==0: zero_found.append((n,m,k)); ok4=False
print(f"C4 all indexed points are common roots and Jacobian != 0 (n<=24, {cnt} points):","PASS" if ok4 else f"FAIL {zero_found[:5]}")
# degrees/top parts: P1,P0 total degree d1,d0 and top parts nonzero
okd=True
for n in range(1,31):
    P1=strip(A[n+1],(n+1)%2); Phi0=padd(pmulm(A[n],(1,0),F(1,n+2)),pmulc(A[n],F(-2))); P0=strip(Phi0,n%2)
    if max(i+j for (i,j) in P1)!=(n+1)//2 or max(i+j for (i,j) in P0)!=n//2+1: okd=False
print("C6 total degrees of stripped limit polynomials are d1=floor((n+1)/2), d0=floor(n/2)+1 (n<=30):","PASS" if okd else "FAIL")
