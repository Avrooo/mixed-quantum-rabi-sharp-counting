"""Conceptual workflow figure (grayscale, 900 dpi PNG + vector PDF)."""
import os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"..","figures")
plt.rcParams.update({"font.family":"serif","font.size":7,"pdf.fonttype":42,"mathtext.fontset":"cm"})
fig,ax=plt.subplots(figsize=(6.5,3.1)); ax.set_xlim(0,100); ax.set_ylim(-8,50); ax.axis("off")
def box(x,y,w,h,txt,fc="white",bold=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.25,rounding_size=1.2",fc=fc,ec="black",lw=0.8))
    ax.text(x+w/2,y+h/2,txt,ha="center",va="center",fontsize=5.6,fontweight="bold" if bold else "normal",linespacing=1.25)
def arr(x0,y0,x1,y1,ls="-"):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle="-|>",mutation_scale=7,lw=0.8,color="black",ls=ls,shrinkA=0,shrinkB=0))
def line(xs,ys): ax.plot(xs,ys,color="black",lw=0.8,solid_capstyle="butt")
W,H=18.6,14; xs=[0.5,20.5,40.5,60.5,80.5]; y1=30; y2=8
ax.text(50,47.5,"Path 1: upper bound (valid for every bias, $0<g_2<\\omega/2$)",ha="center",fontsize=6.6,style="italic")
top=["Model and\nSegal–Bargmann\nrepresentation\n(Sec. II)","Termination\nconditions\n$Q_1^{(n)}=Q_0^{(n)}=0$\n(Sec. III)","Degrees $d_1,d_0$ and\nhighest-degree parts\n$G_1,G_0$\n(Sec. IV)","Exponent argument:\nno common root of\n$G_1,G_0$ at infinity\n(Sec. V)","Bézout count:\nat most $B_n=d_1d_0$\nadmissible pairs\n(Sec. IV)"]
for x,tx in zip(xs,top): box(x,y1,W,H,tx,fc="0.9" if x==80.5 else "white",bold=(x==80.5))
for a,b in zip(xs[:-1],xs[1:]): arr(a+W+0.3,y1+H/2,b-0.3,y1+H/2)
ax.text(50,-6.8,"Path 2: attainment (nonempty region)  |  Path 3: negative coupling (does not establish sharpness)",ha="center",fontsize=6.4,style="italic")
box(xs[1],y2,W,H,"Limiting system\n$t\\to0,\\ \\varepsilon\\to-\\infty$:\nKrawtchouk roots\n(Sec. VI)")
box(xs[2],y2,W,H,"Implicit-function\ncontinuation to\nfinite $t$ and bias\n(Sec. VI)")
box(xs[4],y2,W,H,"Attainment of $B_n$\nin a nonempty\nregion\n(Sec. VI)",fc="0.9",bold=True)
box(xs[3],y2,W,H,"Negative coupling:\nexceptional coupling\nset $\\Lambda_n$; combined\nbound (Sec. VIII)")
arr(xs[1]+W/2,y1-0.3,xs[1]+W/2,y2+H+0.3)
arr(xs[1]+W+0.3,y2+H/2,xs[2]-0.3,y2+H/2)
# B -> Attainment via lower channel
line([xs[2]+W/2,xs[2]+W/2],[y2-0.3,-2.5]); line([xs[2]+W/2,xs[4]+W/2],[-2.5,-2.5]); arr(xs[4]+W/2,-2.5,xs[4]+W/2,y2-0.3)
arr(xs[4]+W/2,y1-0.3,xs[4]+W/2,y2+H+0.3)
ax.text(xs[4]+W/2+1,(y1+y2+H)/2-0.2,"complete\ncount",fontsize=5.4,va="center",ha="left")
arr(xs[3]+W/2,y1-0.3,xs[3]+W/2,y2+H+0.3)
ax.text(xs[3]+W/2+1,(y1+y2+H)/2-0.2,"fails for\n$t<0$ on $\\Lambda_n$",fontsize=5.4,va="center",ha="left")
fig.savefig(os.path.join(OUT,"fig0_workflow.png"),dpi=900,bbox_inches="tight"); fig.savefig(os.path.join(OUT,"fig0_workflow.pdf"),bbox_inches="tight")
