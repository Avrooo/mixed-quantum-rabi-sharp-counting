"""Fig. 8: structure of exact eigenstates (Table V states, t=1/10, omega=1). All panels from the exact analytic states;
Fock diagonalization used only as an independent check. Grayscale-safe palette, 900 dpi PNG + vector PDF."""
import os,sys,io,contextlib
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import factorial,sqrt
from scipy.linalg import eigh
from scipy.special import eval_hermite,gammaln
with contextlib.redirect_stdout(io.StringIO()):
    import validation_table_run as V      # re-uses certified refinement (40-digit Newton) and fockH/Eexact
V.g2=2*V.tf/(1+4*V.tf**2)   # module left g2 sign-flipped after its unitary-image check
ex=V.ex; fockH=V.fockH; Eexact=V.Eexact; tf=V.tf; ref=V.ref
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"..","figures")
plt.rcParams.update({"font.family":"serif","font.size":9,"axes.linewidth":0.7,"pdf.fonttype":42,"mathtext.fontset":"cm"})
def savefig(fig,name):
    fig.savefig(os.path.join(OUT,name+".png"),dpi=900,bbox_inches="tight"); fig.savefig(os.path.join(OUT,name+".pdf"),bbox_inches="tight"); plt.close(fig)

def analytic(n,e,Dl,g,K):
    c=ex.core(n,tf,Dl,e,g); a=c['a']; b=c['b']; q=1+4*tf*tf; r=1-2*tf; beta=q/r*g
    ser=np.zeros(K); ser[0]=1.0
    for k in range(1,K): ser[k]=(beta*ser[k-1]+(2*tf*ser[k-2] if k>=2 else 0))/k
    mul=lambda P: sum(np.pad(pj*ser[:K-j],(j,0)) for j,pj in enumerate(P))[:K]
    fk=np.exp(0.5*gammaln(np.arange(K)+1))
    p=np.concatenate([mul(a)*fk,mul(b)*fk]); return p/np.linalg.norm(p)
def split(p,K): return p[:K],p[K:]
def hermfun(K,x):
    """Fock wavefunctions phi_k(x), k<K, via stable recurrence (omega=m=hbar=1)."""
    P=np.zeros((K,len(x))); P[0]=np.pi**-0.25*np.exp(-x*x/2)
    if K>1: P[1]=np.sqrt(2)*x*P[0]
    for k in range(2,K): P[k]=np.sqrt(2/k)*x*P[k-1]-np.sqrt((k-1)/k)*P[k-2]
    return P
K=160; sel=[1,3,6]   # indices into table: use (n,eps) pairs
pick=[(1,0.5),(3,0.0),(6,0.0)]
rows=[r for r in ref if (r[0],r[1]) in pick]
rows=[next(r for r in ref if r[0]==n and abs(r[1]-e)<1e-12) for n,e in pick]
x=np.linspace(-8,8,1601); dx=x[1]-x[0]; Phi=hermfun(K,x)
data=[]; log=[]
for (n,e,Dl,g) in rows:
    E=Eexact(n,e,g); p=analytic(n,e,Dl,g,K); up,dn=split(p,K)
    H=fockH(K,Dl,e,g); w,Vv=eigh(H); idx=int(np.argmin(abs(w-E))); pn=Vv[:,idx]; pn=pn*np.sign(pn@p)
    un,dnn=split(pn,K)
    Pk_a=up**2+dn**2; Pk_n=un**2+dnn**2
    su=Phi.T@up; sd=Phi.T@dn          # analytic position spinor
    nu=Phi.T@un; nd=Phi.T@dnn
    rho_up=su**2; rho_dn=sd**2
    # Wigner of reduced boson state, spin summed
    M=int(round(12/dx)); W=np.zeros((len(x),241)); ps=np.linspace(-6,6,241)
    for sp in (su,sd):
        for i in range(len(x)):
            k=np.arange(-min(i,len(x)-1-i,M),min(i,len(x)-1-i,M)+1)
            f=sp[i+k]*sp[i-k]
            W[i]+= (dx/np.pi)*(f[None,:]*np.cos(2*ps[:,None]*k[None,:]*dx)).sum(1)
    dp=ps[1]-ps[0]
    nfock=np.arange(K)
    rec=dict(n=n,e=e,D=Dl,g=g,E=E,norm=np.linalg.norm(p),Pa=Pk_a,Pn=Pk_n,x=x,ru=rho_up,rd=rho_dn,ps=ps,W=W,
             nbar=(nfock*Pk_a).sum(),pup=(up**2).sum(),pdn=(dn**2).sum(),
             dens_err=max(abs(rho_up-nu**2).max(),abs(rho_dn-nd**2).max()),
             P_err=abs(Pk_a-Pk_n).max(),fid=abs(pn@p),dE=abs(w[idx]-E),
             tail=(Pk_a[-20:]).sum(),dens_norm=(rho_up.sum()+rho_dn.sum())*dx,
             Wint=W.sum()*dx*dp,Wmin=W.min(),Wneg=-(W[W<0].sum())*dx*dp,
             margx=np.abs(W.sum(1)*dp-(rho_up+rho_dn)).max())
    # cutoff convergence
    conv=[]
    for Kc in (8,12,16,24,32,48,64,96):
        pk=np.concatenate([p[:K][:Kc],p[K:][:Kc]]); pk/=np.linalg.norm(pk)
        Hc=fockH(Kc,Dl,e,g); wc,Vc=eigh(Hc); i=int(np.argmin(abs(wc-E)))
        conv.append((Kc,1-abs(Vc[:,i]@pk),abs(wc[i]-E)))
    rec['conv']=conv; rec['p']=p; rec['pn']=pn; rec['K']=K; data.append(rec)
    print("n=%d eps=%.1f Delta=%.10f g1=%.10f E=%.9f | |psi|=%.15f nbar=%.4f Pup=%.4f Pdn=%.4f tail=%.1e"%(n,e,Dl,g,E,rec['norm'],rec['nbar'],rec['pup'],rec['pdn'],rec['tail']))
    print("   dens norm %.2e, max|rho_an-rho_num|=%.1e, max|P_an-P_num|=%.1e, 1-fid=%.1e, |dE|=%.1e"%(rec['dens_norm']-1,rec['dens_err'],rec['P_err'],1-rec['fid'],rec['dE']))
    print("   Wigner: int=%.6f  Wmin=%.4f  neg.vol=%.4f  max marginal err=%.1e"%(rec['Wint'],rec['Wmin'],rec['Wneg'],rec['margx']))
    print("   conv:", ["K=%d:%.0e/%.0e"%c for c in conv])
import pickle,csv; pickle.dump(data,open(os.path.join(HERE,"data","fig8.pkl"),"wb"))
os.makedirs(os.path.join(HERE,"data","fig8"),exist_ok=True)
with open(os.path.join(HERE,"data","fig8","fig8_summary.csv"),"w",newline="") as fh:
    wr=csv.writer(fh); wr.writerow(["n","eps","Delta","g1","E_exact","N_cut","norm_analytic","weight_up","weight_down","mean_photon_number","Wigner_min","Wigner_negative_volume","Wigner_integral","max_abs_dP_analytic_vs_Fock","max_abs_drho_analytic_vs_Fock","one_minus_overlap","abs_dE_Fock","truncated_weight_last20"])
    for d in data: wr.writerow([d['n'],d['e'],repr(d['D']),repr(d['g']),repr(d['E']),d['K'],repr(float(d['norm'])),repr(float(d['pup'])),repr(float(d['pdn'])),repr(float(d['nbar'])),repr(float(d['Wmin'])),repr(float(d['Wneg'])),repr(float(d['Wint'])),repr(float(d['P_err'])),repr(float(d['dens_err'])),repr(float(1-d['fid'])),repr(float(d['dE'])),repr(float(d['tail']))])
    wr.writerow([])
with open(os.path.join(HERE,"data","fig8","fig8_cutoff_convergence.csv"),"w",newline="") as fh:
    wr=csv.writer(fh); wr.writerow(["n","N_cut","one_minus_overlap","abs_dE"])
    for d in data:
        for (Kc,f_,e_) in d['conv']: wr.writerow([d['n'],Kc,repr(float(f_)),repr(float(e_))])
for d in data:
    np.savez_compressed(os.path.join(HERE,"data","fig8","fig8_state_n%d.npz"%d['n']),n=d['n'],eps=d['e'],Delta=d['D'],g1=d['g'],E=d['E'],
        fock_coeffs_analytic=d['p'],fock_coeffs_fockdiag=d['pn'],layout="[spin-up Fock 0..K-1, spin-down Fock 0..K-1]",
        P_analytic=d['Pa'],P_fockdiag=d['Pn'],x=d['x'],rho_up=d['ru'],rho_down=d['rd'],p_grid=d['ps'],W_xp=d['W'])

# ---------------- figure ----------------
C_UP="#0072B2"; C_DN="#D55E00"
fig,axs=plt.subplots(3,3,figsize=(7.0,6.3))
for j,d in enumerate(data):
    ax=axs[0,j]; ax.plot(d['x'],d['ru'],color=C_UP,lw=1.3,label=r"$|\psi_\uparrow|^2$"); ax.plot(d['x'],d['rd'],color=C_DN,lw=1.3,ls=(0,(4,2)),label=r"$|\psi_\downarrow|^2$")
    ax.set_xlim(-4.5,4.5); ax.set_xlabel(r"$x$"); 
    if j==0: ax.set_ylabel(r"$|\psi_\sigma(x)|^2$")
    ax.set_title(r"(%s) $n=%d$, $\varepsilon=%s$"%("abc"[j],d['n'],("%g"%d['e'])),fontsize=8.5)
    if j==0: ax.legend(frameon=False,fontsize=7,loc="upper left")
    ax=axs[1,j]; kk=np.arange(40); kmax=int(min(40,max(12,np.nonzero(d['Pa']>1e-6)[0].max()+2)))
    ax.bar(kk[:kmax],d['Pa'][:kmax],width=0.8,color="0.78",edgecolor="0.3",lw=0.5,label="analytic")
    ax.plot(kk[:kmax],d['Pn'][:kmax],"o",mfc="none",mec="k",ms=3.5,mew=0.8,label="Fock diag.")
    ax.set_xlabel(r"$k$"); ax.set_yscale("log"); ax.set_ylim(1e-6,1)
    if j==0: ax.set_ylabel(r"$P(k)$"); ax.legend(frameon=False,fontsize=7,loc="upper right")
    ax=axs[2,j]; vm=0.32
    im=ax.pcolormesh(d['ps'],d['x'],d['W'],cmap="RdBu_r",vmin=-vm,vmax=vm,shading="auto",rasterized=True)
    ax.set_xlim(-4.5,4.5); ax.set_ylim(-4.5,4.5); ax.set_aspect("equal"); ax.set_xlabel(r"$p$"); 
    if j==0: ax.set_ylabel(r"$x$")
    wm=d['Wmin']
    if abs(wm)>=1e-2: lab=r"$W_{\min}=%.3f$"%wm
    else:
        ex_=int(np.floor(np.log10(abs(wm)))); lab=r"$W_{\min}=-%.1f\times10^{%d}$"%(abs(wm)/10**ex_,ex_)
    ax.text(0.04,0.95,lab,transform=ax.transAxes,fontsize=7,va="top")
    ax.set_title("") if False else None
for r_,row in enumerate("abc"):
    for j in range(3): axs[r_,j].text(-0.0,1.0,"("+"abcdefghi"[3*r_+j]+")",transform=axs[r_,j].transAxes,fontsize=8,va="bottom",ha="left",fontweight="bold") if False else axs[r_,j].annotate("("+"abcdefghi"[3*r_+j]+")",xy=(0.02,1.0),xycoords="axes fraction",xytext=(0,2),textcoords="offset points",fontsize=8.5,va="bottom") if r_>0 else None
fig.subplots_adjust(hspace=0.5,wspace=0.38,right=0.9)
cax=fig.add_axes([0.92,0.08,0.015,0.22]); cb=fig.colorbar(im,cax=cax); cb.set_label(r"$W(x,p)$",fontsize=8); cb.ax.tick_params(labelsize=7)
savefig(fig,"fig8_eigenstate_structure")
