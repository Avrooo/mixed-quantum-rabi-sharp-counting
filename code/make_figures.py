"""Generates all figures (vector PDF, black & white) and LaTeX tables of the manuscript."""
import os, csv, math
from fractions import Fraction as F
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"..","figures"); TAB=os.path.join(HERE,"..","tables")
os.makedirs(OUT,exist_ok=True); os.makedirs(TAB,exist_ok=True)
plt.rcParams.update({"font.family":"serif","font.size":9,"axes.linewidth":0.8,"pdf.fonttype":42,"axes.grid":False,
                     "mathtext.fontset":"cm"})
def savefig(fig,name):
    fig.savefig(os.path.join(OUT,name),bbox_inches="tight")
    fig.savefig(os.path.join(OUT,name.replace(".pdf",".png")),dpi=900,bbox_inches="tight"); plt.close(fig)

# ---------- Fig. 1: limiting roots (U,w), nested by level m
marks=["o","s","^","D","v","P"]; fills=["white","black","white","black","white","black"]
fig,ax=plt.subplots(figsize=(3.4,3.0))
for m in range(1,7):
    s=m+2; pts=[(2*s, 8*k*k/(s*s-k*k)) for k in range(m,0,-2)]
    ax.plot([p[0] for p in pts],[p[1] for p in pts],marks[m-1],mfc=fills[m-1],mec="black",ms=5,ls="none",label=fr"$\ell={m}$")
ax.set_yscale("log"); ax.set_xlabel(r"$U=\Delta^2/|\varepsilon|$ (limit $t\to0$)"); ax.set_ylabel(r"$w=g_1^2/(t|\varepsilon|)$")
ax.set_xticks([6,8,10,12,14,16]); ax.legend(frameon=False,ncol=2,fontsize=7,loc="upper left")
savefig(fig,"fig1_limit_roots.pdf")

# ---------- Fig. 2: lambda_+(t)
def lam(t):
    q=1+4*t*t; r=1-2*t; m=1+2*t
    z=q/(q-8*t+0j); return 0.5-(m/(2*r))*np.sqrt(z)
tstar=1-math.sqrt(3)/2
fig,ax=plt.subplots(figsize=(3.4,3.0))
tn=np.linspace(-0.4999,-1e-4,800); ax.plot(tn,lam(tn).real,"k-",lw=1.2,label=r"$\lambda_+$, $t<0$")
tp1=np.linspace(1e-4,tstar-2e-3,400); ax.plot(tp1,lam(tp1).real,"k--",lw=1.0,label=r"$\lambda_+<0$, $0<t<t_*$")
tp2=np.linspace(tstar+2e-3,0.4999,400); ax.plot(tp2,lam(tp2).real,"k-.",lw=1.0,label=r"Re $\lambda_+=\frac{1}{2}$, $t>t_*$")
for fr in sorted({F(k,M) for M in range(2,6) for k in range(1,M) if 2*k<M and math.gcd(k,M)==1}):
    ax.axhline(float(fr),color="0.55",lw=0.5)
    _left=float(fr)>0.3; _y=float(fr)
    ax.text(-0.495 if _left else 0.495,_y+(0.012 if _left or _y>0.22 else -0.012),f"${fr.numerator}/{fr.denominator}$",fontsize=7,va="bottom" if (_left or _y>0.22) else "top",ha="left" if _left else "right")
ax.set_ylim(-1.0,1.0); ax.set_xlim(-0.5,0.5); ax.axvline(0,color="0.7",lw=0.5)
ax.axvline(tstar,color="0.4",lw=0.6,ls=":"); ax.text(tstar+0.01,0.7,r"$t_*\approx0.134$",fontsize=6.5)
ax.plot([-0.168602],[1/3],"ko",ms=3.5,mfc="white"); ax.text(-0.46,0.02,r"$t_p$: $\lambda_+=\frac{1}{3}$ at $t=-0.1686$",fontsize=6.2)
ax.set_xlabel(r"Gaussian parameter $t$"); ax.set_ylabel(r"$\lambda_+(t)$"); ax.legend(frameon=False,fontsize=6.5,loc="lower right",bbox_to_anchor=(1.0,0.25))
savefig(fig,"fig2_lambda_plus.pdf")

# ---------- Fig. 3: n=1 existence threshold eps < A/B (closed forms of Sec. 5)
def AB(t):
    q=1+4*t*t; Q=1+6*t+4*t**2+24*t**3+16*t**4
    A=(9+36*t-24*t**2-432*t**3-96*t**4-1728*t**5-384*t**6+2304*t**7+2304*t**8)/(q*q*Q)
    B=6*(1-4*t*t)*(1+2*t)**2/Q
    return A,B
t=np.linspace(0.005,0.495,500); A,B=AB(t)
print("A/B at t=0.1:",AB(0.1)[0]/AB(0.1)[1],"(expected 1.32657)")
t0=0.308000678
fig,ax=plt.subplots(figsize=(3.4,2.8))
ax.plot(t,A/B,"k-",lw=1.2); ax.axhline(0,color="k",lw=0.5); ax.axvline(t0,color="0.5",lw=0.6,ls="--")
ax.fill_between(t,-5,A/B,color="0.88"); ax.set_ylim(-5,3.5); ax.set_xlim(0,0.5)
ax.set_xlabel(r"Gaussian parameter $t$"); ax.set_ylabel(r"$\varepsilon_{\max}/\omega=A/B$")
ax.text(0.03,-3.2,"Admissible degree-one\nsolutions",fontsize=6.5); ax.text(0.03,2.2,"No admissible degree-one\nsolutions",fontsize=6.5); ax.text(0.315,2.3,r"$t_0=0.3080$",fontsize=7)
savefig(fig,"fig3_n1_threshold.pdf")

# ---------- Fig. 4: count of admissible pairs vs eps at t=1/10 from critical biases (V1 Result 3), checked against certified CSV
crit={1:[1.32657],2:[1.23711,1.76149],3:[1.09091,1.70868,2.17690,2.28925],
      4:[0.89498,1.60869,2.14720,2.26619,2.58167,2.72505]}
def count(n,e): return sum(1 for c in crit[n] if e<c)
rows=list(csv.DictReader(open(os.path.join(HERE,"data","CERTIFIED_ROOT_COUNTS.csv"))))
mism=0;tot=0
for r in rows:
    if r["t"]=="1/10" and int(r["n"]) in crit and r["status"].startswith("certified"):
        tot+=1; mism+= (count(int(r["n"]),float(F(r["eps"])))!=int(r["count"]))
print(f"Fig.4 consistency with certified rows at t=1/10: {tot-mism}/{tot} agree")
fig,ax=plt.subplots(figsize=(3.4,2.8)); ls={1:"-",2:"--",3:"-.",4:":"}
e=np.linspace(-0.5,3.1,3000)
for n in crit:
    ax.step(e,[count(n,x) for x in e],where="post",color="k",ls=ls[n],lw=1.1,label=f"$n={n}$")
    for r in rows:
        if r["t"]=="1/10" and int(r["n"])==n and r["status"].startswith("certified") and float(F(r["eps"]))>-0.5:
            ax.plot(float(F(r["eps"])),int(r["count"]),"ko",ms=3,mfc="white")
ax.set_xlabel(r"$\varepsilon/\omega$"); ax.set_ylabel(r"number of admissible $(\Delta^2,g_1^2)$"); ax.legend(frameon=False,fontsize=7)
ax.set_xlim(-0.5,3.1); ax.set_yticks(range(0,7))
savefig(fig,"fig4_counts_vs_bias.pdf")

# ---------- Fig. 5: B_n and degrees
fig,ax=plt.subplots(figsize=(3.4,2.6)); n=np.arange(1,13)
Bn=(n+1)**2//4; ax.plot(n,Bn,"ks-",ms=4,lw=0.8,label=r"$B_n=\lfloor (n+1)^2/4\rfloor$")
ax.plot(n,(n+1)//2,"k^--",ms=4,lw=0.8,mfc="white",label=r"$d_1=\lfloor (n+1)/2\rfloor$")
ax.plot(n,n//2+1,"ko:",ms=4,lw=0.8,mfc="white",label=r"$d_0=\lfloor n/2\rfloor+1$")
ax.set_xlabel("$n$"); ax.set_ylabel("count / degree"); ax.legend(frameon=False,fontsize=7)
savefig(fig,"fig5_bezout_numbers.pdf")

# ---------- Tables (B&W, booktabs)
with open(os.path.join(TAB,"tab_bn.tex"),"w") as f:
    f.write("\\begin{tabular}{rcccc}\\toprule\n$n$ & $d_1$ & $d_0$ & $B_n=d_1d_0$ & $|\\Lambda_n|$\\\\\\midrule\n")
    def lam_count(n): return len({F(k,M) for M in range(2,n+3) for k in range(1,M) if 2*k<M and math.gcd(k,M)==1})
    for n in range(1,11): f.write(f"{n} & {(n+1)//2} & {n//2+1} & {(n+1)**2//4} & {lam_count(n)}\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n")
exc=list(csv.DictReader(open(os.path.join(HERE,"data","EXACT_EXCEPTIONAL_CASES.csv"))))
firstn={}
for r in exc:
    key=r["k_over_M"]; firstn.setdefault(key,r)
with open(os.path.join(TAB,"tab_lambda.tex"),"w") as f:
    f.write("\\begin{tabular}{cccccc}\\toprule\n$\\lambda_+$ & $M$ & First $n$ & $t$ & $g_2/\\omega$ & $\\chi=u/v$\\\\\\midrule\n")
    items=sorted(firstn.items(),key=lambda kv:float(F(kv[0])))
    for k,r in items:
        fr=F(k); nmin=min(int(x["n"]) for x in exc if x["k_over_M"]==k)
        f.write(f"$\\tfrac{{{fr.numerator}}}{{{fr.denominator}}}$ & {fr.denominator} & {nmin} & ${float(r['t_root']):.6f}$ & ${float(r['g2_over_omega']):.5f}$ & ${float(r['chi_u_over_v_for_smallest_n0']):.4f}$\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n")
print("figures and tables written")
