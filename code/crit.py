import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cert import *
from limit_sys import cert_polys
from ratrec import pmul
def lagr(pts,vals):
    n=len(pts); coef=[F(0)]*n
    for i,xi in enumerate(pts):
        num=[F(1)]; den=F(1)
        for j,xj in enumerate(pts):
            if j!=i: num=pmul(num,[-xj,F(1)]); den*=xi-xj
        for k,c in enumerate(num): coef[k]+=vals[i]*c/den
    return trim(coef)
def fit(f,maxd=120):
    pts=[F(k,3)-F(7,5) for k in range(maxd+6)]; cache={}
    def val(i):
        if i not in cache: cache[i]=f(pts[i])
        return cache[i]
    for d in range(1,maxd):
        p=lagr(pts[:d+1],[val(i) for i in range(d+1)])
        if all(peval(p,pts[d+1+j])==val(d+1+j) for j in range(3)): return p
    return None
def eliminants(n,t):
    c1,c0=reduced_coeffs(n,t)
    def mk(e):
        P=geom.fix_eps(c1,e); Q=geom.fix_eps(c0,e); return geom.resultant_v(P,Q),geom.resultant_v(swap(P),swap(Q))
    return c1,c0,mk
def run(n,t):
    c1,c0,mk=eliminants(n,t)
    def pad(R,k): return R[k] if k<len(R) else F(0)
    D=len(((n+1)**2)//4*[0])
    polys={}
    polys["u=0 (Delta->0)"]=fit(lambda e: pad(mk(e)[1],0))
    polys["v=0 (g1->0)"]=fit(lambda e: pad(mk(e)[0],0))
    polys["lead coeff v-elim (escape)"]=fit(lambda e: mk(e)[0][-1] if mk(e)[0] else F(0))
    polys["lead coeff u-elim (escape)"]=fit(lambda e: mk(e)[1][-1] if mk(e)[1] else F(0))
    def disc(e):
        Rv=mk(e)[0]; 
        from geom import sylv
        return sylv(Rv,deriv(Rv))
    polys["discriminant of v-elim (fold or v-collision)"]=fit(disc,maxd=160)
    return polys
def real_roots(p,lo=F(-60),hi=F(60)):
    if p is None or len(p)<2: return []
    seq=sturm_seq(p); 
    if peval(p,lo)==0 or peval(p,hi)==0: return ["boundary"]
    k=var_at(seq,lo)-var_at(seq,hi)
    if k==0: return []
    out=[]; 
    def rec(a,b):
        c=var_at(seq,a)-var_at(seq,b)
        if c==0: return
        if c==1 and b-a<F(1,10**8): out.append(float((a+b)/2)); return
        m=(a+b)/2
        if peval(seq[0],m)==0: out.append(float(m)); return
        rec(a,m); rec(m,b)
    rec(lo,hi); return out
if __name__=="__main__":
    t=F(sys.argv[1]); 
    for n in map(int,sys.argv[2:]):
        polys=run(n,t); allc=[]
        print(f"n={n}, alpha={t}")
        for k,p in polys.items():
            rr=real_roots(p) if p else None
            print(f"   {k}: deg={None if p is None else len(p)-1}, real roots in (-60,60): {rr if rr is None else [round(x,5) for x in rr]}",flush=True)
            if rr: allc+= [(x,k) for x in rr if isinstance(x,float)]
        allc.sort(); pts=sorted(set(round(x,6) for x,_ in allc))
        samples=[F(-100)]+[F(int(round((pts[i]+pts[i+1])/2*10**6)),10**6) for i in range(len(pts)-1)]+[F(100)] if pts else [F(0)]
        print("   critical eps (all mechanisms):",[(round(x,5),k.split()[0]) for x,k in allc])
        print("   certified count in each region (sample point -> count):",[(float(s),certify(n,t,s)['count'] if 'count' in certify(n,t,s) else 0) for s in samples],flush=True)
