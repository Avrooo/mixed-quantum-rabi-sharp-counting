import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys; from fractions import Fraction as F
import random
import exact as ex

random.seed(4)
ts = [F(1,10), F(3,10), F(-1,5), F(2,5), F(-2,5), F(5,2), F(-5,2)]   # last two: the other Gaussian root alpha_l = 1/(4 alpha_s)
bad = 0; tot = 0
for t in ts:
    for n in range(0, 7):
        for trial in range(2):
            x = F(random.randint(3, 40), random.randint(3, 17)); e = F(random.randint(-30, 30), random.randint(3, 17)); g = F(random.randint(1, 30), random.randint(3, 17))
            c = ex.core(n, t, x, e, g)
            r1, r2 = ex.raw_residual(n, t, x, e, g, c['a'], c['b'], c['E'])
            # prediction: r2 == 0 identically; r1[j]==0 for j>=2; r1[1]==R1; r1[0]==R0
            ok = all(v == 0 for v in r2) and all(v == 0 for v in r1[2:]) and r1[1] == c['R1'] and r1[0] == c['R0']
            tot += 1; bad += (not ok)
print(f"exact identity check (raw operators vs closed-form recursion): {tot-bad}/{tot} pass over t in {[str(t) for t in ts]}, n=0..6 (random rational Delta, eps, g1 NOT on the solution set)")
