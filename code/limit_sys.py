import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from cert import *
def leading_forms(n,t):
    c1,c0=reduced_coeffs(n,t)
    out=[]
    for c in (c1,c0):
        d=max(i+j+k for (i,j,k) in c)
        out.append({k:v for k,v in c.items() if sum(k)==d})
    return out
def at_eps(form,e):
    o={}
    for (i,j,k),v in form.items(): o[(i,k)]=o.get((i,k),0)+v*e**j
    return {k:v for k,v in o.items() if v!=0}
def cert_polys(P,Q,width=F(1,10**12)):
    Rv=geom.resultant_v(P,Q); Ru=geom.resultant_v(swap(P),swap(Q))
    if not Rv or not Ru: return dict(status="ZERO_RESULTANT")
    nv,sv=npos(Rv); nu,su=npos(Ru)
    info=dict(res_deg=len(Rv)-1,sqfree=(len(sqfree(Rv))==len(Rv) and len(sqfree(Ru))==len(Ru)),pos_v=nv,pos_u=nu,zero_root=(Rv[0]==0 or Ru[0]==0))
    if nv==0 or nu==0: info.update(count=0); return info
    Iv=isolate_pos(sv,width); Iu=isolate_pos(su,width)
    boxes=[(u,v) for u in Iu for v in Iv if box_contains_zero(P,u,v) and box_contains_zero(Q,u,v)]
    inj=len({b[0] for b in boxes})==len(boxes) and len({b[1] for b in boxes})==len(boxes)
    info.update(count=len(boxes),injective=inj); return info
if __name__=="__main__":
    ts=[F(1,50),F(1,20),F(1,10),F(1,5),F(3,10),F(2,5),F(9,20),F(-1,20),F(-1,10),F(-3,10),F(-2,5)]
    print("limit system F1(U,V,-1)=F0(U,V,-1)=0   (rows: n; entries: certified #positive roots / Bezout bound; '!' = not squarefree or ambiguous)")
    for t in ts:
        row=[]
        for n in range(1,9):
            F1,F0=leading_forms(n,t); r=cert_polys(at_eps(F1,F(-1)),at_eps(F0,F(-1)))
            row.append(f"{r.get('count','?')}/{((n+1)**2)//4}"+("" if r.get('sqfree',False) and r.get('injective',True) else "!"))
        print(f"  alpha={str(t):>6} g2={float(2*t/(1+4*t*t)):+.4f}: "," ".join(row))
