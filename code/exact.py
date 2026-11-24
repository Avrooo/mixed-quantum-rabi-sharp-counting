"""Exact-arithmetic core for the biased mixed Rabi model (omega = 1).
H = D sx + e sz + a^dag a + sz [ g1 (a + a^dag) + g2 (a^2 + a^dag^2) ]
Rational uniformisation of the Gaussian branch:   alpha = t,
   q = 1+4t^2,  r = 1-2t,  m = 1+2t
   g2 = 2t/q,   Omega = (1-4t^2)/q  (= omega - 4 alpha g2),   kappa = q/r
   |t| < 1/2  : normalisable branch (alpha_s);   |t| > 1/2 : the other root alpha_l = 1/(4 alpha_s)  (Omega < 0)
   t -> -t    : the other spinor ordering of the same sign of g2 is NOT this map (see unitary equivalence in the report).
core() works with ints / Fractions / floats / QP polynomials in (x=Delta, e=eps, g=g1).
"""
from fractions import Fraction as F

class QP:
    def __init__(s, d=None): s.d = {k: v for k, v in (d or {}).items() if v != 0}
    @staticmethod
    def var(i):
        k = [0, 0, 0]; k[i] = 1; return QP({tuple(k): 1})
    def _c(s, o): return o if isinstance(o, QP) else QP({(0, 0, 0): o})
    def __add__(s, o):
        o = s._c(o); r = dict(s.d)
        for k, v in o.d.items(): r[k] = r.get(k, 0) + v
        return QP(r)
    __radd__ = __add__
    def __neg__(s): return QP({k: -v for k, v in s.d.items()})
    def __sub__(s, o): return s + (-s._c(o))
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if not isinstance(o, QP): return QP({k: v * o for k, v in s.d.items()})
        r = {}
        for k1, v1 in s.d.items():
            for k2, v2 in o.d.items():
                k = (k1[0]+k2[0], k1[1]+k2[1], k1[2]+k2[2]); r[k] = r.get(k, 0) + v1*v2
        return QP(r)
    __rmul__ = __mul__
    def __truediv__(s, a): return QP({k: v / a for k, v in s.d.items()})

def setup(t):
    q = 1 + 4*t*t; r = 1 - 2*t; m = 1 + 2*t
    return q, r, m, (1 - 4*t*t)/q, 2*t/q, q/r          # q, r, m, Omega, g2, kappa

def core(n, t, x, e, g):
    """closed-form recursion (omega=1). Returns dict(a,b,E,R1,R0)."""
    q, r, m, Om, g2, ka = setup(t)
    beta = ka*g; h = g + 2*beta*g2; s = g*beta + g2*beta*beta
    E = (n+2)*Om - 2*t*g2 - s - e
    p = 4*t
    d = lambda k: k*(1 + 4*t*g2) - (n+2)*Om + 4*t*g2 + 2*s + 2*e
    a = {j: 0 for j in range(-3, n+8)}; b = dict(a)
    b[n+2] = 1; a[n] = -x/p; b[n+1] = -(n+2)*h/Om
    def bj(j): return ((j+1)*h*b[j+1] + g2*(j+2)*(j+1)*b[j+2] - x*a[j])/(Om*(j-n-2))
    for j in range(n+2, 1, -1):
        if j <= n: b[j] = bj(j)
        if j != n+2:
            a[j-2] = -(2*beta*a[j-1] + d(j)*a[j] + (j+1)*h*a[j+1] + g2*(j+2)*(j+1)*a[j+2] + x*b[j])/p
    if n >= 1: b[1] = bj(1)
    if n >= 0: b[0] = bj(0)
    R1 = 2*beta*a[0] + d(1)*a[1] + 2*h*a[2] + 6*g2*a[3] + x*b[1]
    R0 = d(0)*a[0] + h*a[1] + 2*g2*a[2] + x*b[0]
    return dict(a=[a[j] for j in range(n+1)], b=[b[j] for j in range(n+3)], E=E, R1=R1, R0=R0)

def constraints_qp(n, t):
    X, EP, G = QP.var(0), QP.var(1), QP.var(2)
    c = core(n, t, X, EP, G)
    return c['R1'], c['R0'], c['E']

# ---------- independent residual from the raw operators (no band formulas) ----------
def matmul(A, B):
    K = len(A); return [[sum(A[i][k]*B[k][j] for k in range(K) if A[i][k] != 0) for j in range(K)] for i in range(K)]
def madd(A, B, s=1): return [[A[i][j] + s*B[i][j] for j in range(len(A))] for i in range(len(A))]
def mscal(A, c): return [[c*v for v in row] for row in A]
def raw_residual(n, t, x, e, g, a, b, E):
    """apply  M+ p1 + x p2  and  M- p2 + x p1  with M+- built from z, d/dz and the gauge shift d -> d + 2 t z + beta."""
    q, r, m, Om, g2, ka = setup(t); beta = ka*g
    K = n + 10; zero = 0
    Z = [[zero]*K for _ in range(K)]; D = [[zero]*K for _ in range(K)]; I = [[zero]*K for _ in range(K)]
    for k in range(K-1): Z[k+1][k] = 1
    for k in range(1, K): D[k-1][k] = k
    for k in range(K): I[k][k] = 1
    T = madd(madd(D, mscal(Z, 2*t)), mscal(I, beta))
    N = matmul(Z, T)
    A = madd(mscal(madd(Z, T), g), mscal(madd(matmul(Z, Z), matmul(T, T)), g2))
    Mp = madd(madd(N, A), mscal(I, e - E)); Mm = madd(madd(N, A, -1), mscal(I, -e - E))
    p1 = [a[j] if j <= n else 0 for j in range(K)]; p2 = [b[j] if j <= n+2 else 0 for j in range(K)]
    ap = lambda M, v: [sum(M[i][k]*v[k] for k in range(K)) for i in range(K)]
    r1 = [u + x*w for u, w in zip(ap(Mp, p1), p2)]
    r2 = [u + x*w for u, w in zip(ap(Mm, p2), p1)]
    return r1, r2
