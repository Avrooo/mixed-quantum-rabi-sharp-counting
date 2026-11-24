import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
"""Certified exact counting of admissible (Delta^2, g1^2) pairs. Pure Fraction arithmetic (no sympy available)."""
import sys; from fractions import Fraction as F
from ratrec import reduced_coeffs, trim, pdivmod, pmul, peval
import geom
def deriv(p): return trim([i*c for i, c in enumerate(p)][1:])
def pgcd_(a, b):
    while b: a, b = b, pdivmod(a, b)[1]
    return [c/a[-1] for c in a] if a else a
def sqfree(p):
    g = pgcd_(p, deriv(p)); return pdivmod(p, g)[0] if len(g) > 1 else p
def sturm_seq(p):
    p = sqfree(p); seq = [p, deriv(p)]
    while len(seq[-1]) > 1:
        r = pdivmod(seq[-2], seq[-1])[1]
        if not r: break
        seq.append([-c for c in r])
    return seq
def sign_at(p, x):
    v = peval(p, x); return (v > 0) - (v < 0)
def var_at_zero_plus(seq):
    s = []
    for q in seq:
        k = next(i for i, c in enumerate(q) if c != 0); s.append(1 if q[k] > 0 else -1)
    return sum(1 for i in range(len(s)-1) if s[i] != s[i+1])
def var_at_inf(seq):
    s = [1 if q[-1] > 0 else -1 for q in seq]; return sum(1 for i in range(len(s)-1) if s[i] != s[i+1])
def var_at(seq, x):
    s = [sign_at(q, x) for q in seq]; s = [v for v in s if v != 0]
    return sum(1 for i in range(len(s)-1) if s[i] != s[i+1])
def npos(p):
    """exact number of distinct roots in (0, inf) (a root at 0 is excluded, reported separately)"""
    seq = sturm_seq(p); return var_at_zero_plus(seq) - var_at_inf(seq), seq
def isolate_pos(seq, width):
    """disjoint rational intervals each holding exactly one positive root"""
    hi = F(1)
    while var_at(seq, hi) > 0: hi *= 2          # all positive roots below hi
    out = []
    def rec(a, b, va, vb):
        c = va - vb
        if c == 0: return
        if c == 1 and b - a <= width: out.append((a, b)); return
        m = (a+b)/2
        if peval(seq[0], m) == 0: m += width/8
        rec(a, m, va, var_at(seq, m)); rec(m, b, var_at(seq, m), vb)
    a = F(1, 10**40); rec(a, hi, var_at(seq, a), var_at(seq, hi)); return out
def swap(P): return {(k, i): c for (i, k), c in P.items()}
def box_contains_zero(P, ui, vi):
    lo = hi = F(0)
    for (i, k), c in P.items():
        ulo, uhi = ui[0]**i, ui[1]**i; vlo, vhi = vi[0]**k, vi[1]**k
        a, b = ulo*vlo, uhi*vhi
        lo += min(c*a, c*b); hi += max(c*a, c*b)
    return lo <= 0 <= hi
def certify(n, t, e, width=F(1, 10**12)):
    c1, c0 = reduced_coeffs(n, t); P = geom.fix_eps(c1, e); Q = geom.fix_eps(c0, e)
    Rv = geom.resultant_v(P, Q); Ru = geom.resultant_v(swap(P), swap(Q))
    d1 = max(i+k for (i, k) in P); d0 = max(i+k for (i, k) in Q)
    info = dict(n=n, t=str(t), eps=str(e), bezout=((n+1)**2)//4, degP=d1, degQ=d0)
    if not Rv or not Ru: info.update(status="ZERO_RESULTANT_possible_common_component"); return info
    nv, sv = npos(Rv); nu, su = npos(Ru)
    info.update(res_v_deg=len(Rv)-1, res_u_deg=len(Ru)-1, res_v_sqfree=(len(sqfree(Rv)) == len(Rv)), res_u_sqfree=(len(sqfree(Ru)) == len(Ru)), pos_v_roots=nv, pos_u_roots=nu)
    if nv == 0 or nu == 0: info.update(count=0, status="certified_empty_by_Sturm"); return info
    Iv = isolate_pos(sv, width); Iu = isolate_pos(su, width)
    # exclusion: keep boxes where both polynomials can vanish
    boxes = [(u, v) for u in Iu for v in Iv if box_contains_zero(P, u, v) and box_contains_zero(Q, u, v)]
    # every kept box must be the unique candidate for its u-interval and v-interval (injective pairing)
    us = [b[0] for b in boxes]; vs = [b[1] for b in boxes]
    inj = len(set(us)) == len(us) and len(set(vs)) == len(vs)
    info.update(count=len(boxes), pairing_injective=inj, status="certified" if inj else "AMBIGUOUS_pairing")
    info["u_roots"] = [float((b[0][0]+b[0][1])/2) for b in boxes]; info["v_roots"] = [float((b[1][0]+b[1][1])/2) for b in boxes]
    return info
if __name__ == "__main__":
    import time
    for n in range(1, 5):
        t0 = time.time(); r = certify(n, F(1, 10), F(0)); print(n, {k: r[k] for k in r if k not in ("u_roots", "v_roots")}, f"{time.time()-t0:.1f}s")
