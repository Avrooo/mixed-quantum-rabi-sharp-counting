import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys,csv,time; from exc_eps import *
from exceptional import farey
def mx(p): return max(abs(float(c)) for c in p)
rows=[]
for n in map(int,sys.argv[1:]):
    for target in farey(n):
        t0=time.time()
        try:
            t=solve_t(target); c1,c0,mk=eliminants(n,t)
            m_u=(lambda e: mk(e)[1])
            deg=max(len(m_u(F(k))) for k in range(0,3))
            polys=[fit(lambda e,i=i: pad(m_u(e),i),maxd=80) for i in range(deg)]
            S=max(mx(p) for p in polys if p)
            live=[p for p in polys if p and mx(p)>1e-25*S]
            nz=len(live); zero=sum(1 for p in polys if (not p) or mx(p)<=1e-25*S)
            # common real roots among live polys (use the live poly of smallest positive degree)
            cand=sorted([p for p in live if len(p)>1],key=len)
            comm=[]
            if cand:
                L=F(10**4)
                for r in real_roots(cand[0],-L,L):
                    rr=F(r).limit_denominator(10**14)
                    worst=max(abs(float(peval(p,rr)))/mx(p) for p in live)
                    if worst<1e-8: comm.append((float(rr),worst))
            status="COMMON-COMPONENT real eps found" if comm else "no common real eps (no real common component)"
            rows.append(dict(n=n,lambda_plus=str(target),t=f"{float(t):.12f}",vanishing_coeffs=zero,live_coeffs=nz,common_real_eps=";".join(f"{a:.6f}" for a,_ in comm),status=status))
            print(rows[-1],f"[{time.time()-t0:.0f}s]",flush=True)
        except Exception as ex:
            print("ERR",n,target,repr(ex)[:80],flush=True)
w=csv.DictWriter(open("EXCEPTIONAL_EPS_SCAN.csv","w",newline=""),fieldnames=["n","lambda_plus","t","vanishing_coeffs","live_coeffs","common_real_eps","status"]); w.writeheader(); w.writerows(rows)
