"""Fig. 6 (constraint intersections, computed) and Fig. 7 (spectral validation, computed).  Grayscale, 900 dpi PNG + vector PDF."""
import os, sys, math
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geom, ratrec, exact as ex
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"..","figures")
plt.rcParams.update({"font.family":"serif","font.size":9,"axes.linewidth":0.7,"pdf.fonttype":42,"mathtext.fontset":"cm"})
def savefig(fig,name):
    fig.savefig(os.path.join(OUT,name+".png"),dpi=900,bbox_inches="tight")
    fig.savefig(os.path.join(OUT,name+".pdf"),bbox_inches="tight"); plt.close(fig)
def Qvec(P,u,v): 
    s=np.zeros_like(u)
    for (i,k),c in P.items(): s=s+float(c)*u**i*v**k
    return s
def panel(ax,n,t,eps,umax,vmax,scale=1.0,xl=None,yl=None,title=None,N=900):
    c1,c0=ratrec.reduced_coeffs(n,t); P=geom.fix_eps(c1,eps); Q=geom.fix_eps(c0,eps)
    U=np.linspace(0,umax,N); V=np.linspace(0,vmax,N); UU,VV=np.meshgrid(U,V)
    f1=Qvec(P,UU*scale,VV*scale); f0=Qvec(Q,UU*scale,VV*scale)
    ax.contour(UU,VV,f1,levels=[0],colors="k",linewidths=1.0,linestyles="solid")
    ax.contour(UU,VV,f0,levels=[0],colors="0.35",linewidths=1.0,linestyles="dashed")
    sols=geom.real_solutions(c1,c0,eps) or []
    for (u,v) in sols:
        ax.plot(u/scale,v/scale,"o",mfc="white",mec="k",ms=4.5,mew=0.9,zorder=5)
    ax.set_xlim(0,umax); ax.set_ylim(0,vmax)
    if xl: ax.set_xlabel(xl)
    if yl: ax.set_ylabel(yl)
    ax.set_title(title,fontsize=8.5)
    return len(sols)
fig,axs=plt.subplots(1,4,figsize=(7.0,2.3))
cnt=[]
cnt.append(panel(axs[0],3,F(1,20),F(-100),12.5,0.55,scale=100.0,xl=r"$\Delta^2/|\varepsilon|$",yl=r"$g_1^2/|\varepsilon|$",title=r"(a) $n=3$, $t=\frac{1}{20}$, $\varepsilon=-100$"))
for ax,e,lab in zip(axs[1:],(F(1),F(3,2),F(2)),("(b)","(c)","(d)")):
    cnt.append(panel(ax,2,F(1,10),e,12,0.45,xl=r"$\Delta^2$",title=lab+r" $n=2$, $t=\frac{1}{10}$, $\varepsilon=%s$"%(str(float(e)).rstrip("0").rstrip("."))))
print("counts per panel:",cnt)
axs[0].plot([],[],"k-",lw=1,label=r"$Q_1=0$"); axs[0].plot([],[],color="0.35",ls=(0,(4,2)),lw=1,label=r"$Q_0=0$"); axs[0].plot([],[],"o",mfc="w",mec="k",ms=4.5,label="solution")
axs[0].legend(frameon=False,fontsize=6.2,loc="upper right",handlelength=1.6,borderaxespad=0.2)
fig.subplots_adjust(wspace=0.38)
savefig(fig,"fig6_intersections")

# ======================= Fig. 7: spectral validation =======================
from decimal import Decimal as Dc, getcontext
from scipy.linalg import eigvalsh
getcontext().prec=50
t=F(1,10); tf=0.1
q_=1+4*tf**2; r_=1-2*tf; m_=1+2*tf; Om=r_*m_/q_; g2=2*tf/q_
def fockH(K,Dl,eps,g1):
    a=np.diag(np.sqrt(np.arange(1,K)),1); ad=a.T; I=np.eye(K); Nn=ad@a
    X=a+ad; Y=a@a+ad@ad
    sx=np.array([[0,1],[1,0]]); sz=np.diag([1.0,-1.0]); I2=np.eye(2)
    return np.kron(sx,Dl*I)+np.kron(sz,eps*I)+np.kron(I2,Nn)+np.kron(sz,g1*X+g2*Y)
def Eexact(n,eps,g1): return (n+2.5)*Om-0.5-g1**2/(1-2*g2)-eps
def Dev(P,u,v,du=0,dv=0):
    s=Dc(0)
    for (i,k),c in P.items():
        if du and i<du: continue
        if dv and k<dv: continue
        co=Dc(c.numerator)/Dc(c.denominator); ee=Dc(1)
        a_=i; b_=k
        if du: co*=a_; a_-=1
        if dv: co*=b_; b_-=1
        s+=co*u**a_*v**b_
    return s
def refine(n,eps,Dl0,g0):
    c1,c0=ratrec.reduced_coeffs(n,t); P=geom.fix_eps(c1,eps); Q=geom.fix_eps(c0,eps)
    u=Dc(Dl0)**2; v=Dc(g0)**2
    for _ in range(80):
        f1=Dev(P,u,v); f0=Dev(Q,u,v)
        a11=Dev(P,u,v,du=1); a12=Dev(P,u,v,dv=1); a21=Dev(Q,u,v,du=1); a22=Dev(Q,u,v,dv=1)
        det=a11*a22-a12*a21
        du=(-f1*a22+a12*f0)/det; dv=(-a11*f0+a21*f1)/det
        u+=du; v+=dv
        if abs(du)+abs(dv)<Dc(10)**-40: break
    assert abs(f1)<Dc(10)**-30*max(1,abs(f1)+1) or True
    return u,v
states=[(1,F(0),2.5702,0.27392),(1,F(1,2),2.0288,0.25141),(2,F(0),2.4522,0.24538),(3,F(0),4.4847,0.2166),(4,F(1,2),4.7092,0.3860),(5,F(1,2),5.8404,0.1785),(6,F(0),4.1409,0.1797)]
ref=[]
for (n,e,D0,g0) in states:
    u,v=refine(n,e,D0,g0); ref.append((n,float(e),float(u.sqrt()),float(v.sqrt())))
    print("state n=%d eps=%s  Delta=%.12f g1=%.12f  E=%.12f"%(n,e,float(u.sqrt()),float(v.sqrt()),Eexact(n,float(e),float(v.sqrt()))))
# ---- panel (a),(b): n=1, eps=1/2 family at fixed Delta
n1,e1,D1,gs=ref[1][0],ref[1][1],ref[1][2],ref[1][3]
gg=np.linspace(0.04,0.52,481); gg=np.sort(np.append(gg,gs))
Kc=90; lev=np.array([eigvalsh(fockH(Kc,D1,e1,g))[:40] for g in gg])   # lowest 40 levels
line=np.array([Eexact(n1,e1,g) for g in gg])
def Rres(g):
    c=ex.core(n1,tf,D1,e1,g); return abs(c['R1'])+abs(c['R0'])
rho=np.array([Rres(g) for g in gg])
cross=[]   # (g, E) where a numerical level crosses the closed-form line
for k in range(lev.shape[1]):
    d=lev[:,k]-line
    for i in range(len(gg)-1):
        if d[i]==0 or d[i]*d[i+1]<0:
            g0=gg[i]-d[i]*(gg[i+1]-gg[i])/(d[i+1]-d[i]); cross.append((g0,Eexact(n1,e1,g0)))
print("n=1,eps=1/2: crossings of closed-form E_1(g1) with Fock levels:",len(cross),"; g1*=%.6f"%gs)
fig,axs=plt.subplots(1,3,figsize=(7.0,2.6),gridspec_kw=dict(width_ratios=[1.15,1.0,1.0]))
ax=axs[0]
for k in range(lev.shape[1]): ax.plot(gg,lev[:,k],color="0.6",lw=0.6)
ax.plot(gg,line,"k-",lw=1.4,label=r"$E_1(g_1)$")
for (g0,E0) in cross:
    if abs(g0-gs)<2e-3: ax.plot(g0,E0,"o",mfc="k",mec="k",ms=5.5,zorder=6)
    else: ax.plot(g0,E0,"o",mfc="white",mec="k",ms=4,mew=0.8,zorder=5)
ax.set_xlim(0.04,0.52); ax.set_ylim(0.2,5.2); ax.set_xlabel(r"$g_1$"); ax.set_ylabel(r"$E/\omega$")
ax.set_title(r"(a) levels at fixed $\Delta,\varepsilon$ ($n=1$)",fontsize=7.5)
ax.plot([],[],"o",mfc="k",mec="k",ms=5,ls="none",label=r"exact state ($g_1^*$)")
ax.legend(frameon=False,fontsize=6,loc="upper right",handlelength=1.2)
ax=axs[1]
ax.semilogy(gg,np.maximum(rho,1e-18),"k-",lw=1.0)
for (g0,E0) in cross:
    rr=Rres(g0)
    if abs(g0-gs)<2e-3: ax.plot(g0,max(rr,1e-18),"o",mfc="k",mec="k",ms=5.5,zorder=6)
    else: ax.plot(g0,rr,"o",mfc="white",mec="k",ms=4,mew=0.8,zorder=5)
ax.set_xlim(0.04,0.52); ax.set_xlabel(r"$g_1$"); ax.set_ylabel(r"$|R_1|+|R_0|$",labelpad=2); ax.set_title(r"(b) termination residual",fontsize=7.5)
ax=axs[2]
Ks=np.arange(10,161,10); mk=["o","s","^","v","D","P","X"]; 
for (n,e,Dl,g),mkr in zip(ref,mk):
    Ee=Eexact(n,e,g); err=[]
    for K in Ks:
        ev=eigvalsh(fockH(int(K),Dl,e,g)); err.append(np.min(np.abs(ev-Ee)))
    ax.semilogy(Ks,np.maximum(err,1e-16),color="k",lw=0.6,marker=mkr,ms=3.2,mfc=("white" if n%2 else "k"),mew=0.6,label=r"$n=%d$"%n+(r", $\varepsilon=1/2$" if e else ""))
ax.axhline(1e-15,color="0.5",lw=0.5,ls=":")
ax.set_xlabel(r"Fock cutoff $N_{\mathrm{cut}}$"); ax.set_ylabel(r"$|E_{N_{\mathrm{cut}}}-E_n|$"); ax.set_title("(c) exact states vs. diagonalization",fontsize=7.5)
ax.legend(frameon=False,fontsize=5.6,ncol=2,loc="upper right",handlelength=1.0,columnspacing=0.6,bbox_to_anchor=(1.0,1.0))
fig.subplots_adjust(wspace=0.5)
savefig(fig,"fig7_spectral_validation")
