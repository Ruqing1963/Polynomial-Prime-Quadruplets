import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, sys
from math import gcd
def divs(n): return [d for d in range(1,n+1) if n%d==0]
def draw(ax,N,title,col,note):
    ds=divs(N); import math
    lev={d:sum(e for p,e in fac(d).items()) for d in ds}
    L=max(lev.values()); rows={}
    for d in ds: rows.setdefault(lev[d],[]).append(d)
    pos={}
    for k,v in rows.items():
        for i,d in enumerate(sorted(v)): pos[d]=((i-(len(v)-1)/2)*1.6,k*1.4)
    for a in ds:
        for b in ds:
            if b%a==0 and a!=b and lev[b]==lev[a]+1:
                ax.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color=col,lw=2,zorder=1)
    for d,(x,y) in pos.items():
        lab=r"$\{1\}$" if d==1 else (r"$\mathbb{Z}_{%d}^{*}$"%(N+1) if d==N else r"$H_{%d}$"%d)
        ax.scatter([x],[y],s=1500,color=col,edgecolor="k",lw=2,zorder=2)
        ax.text(x,y,lab,ha="center",va="center",color="white",fontsize=12,fontweight="bold",zorder=3)
    ax.set_title(title,fontsize=13,fontweight="bold"); ax.axis("off")
    ax.set_xlim(-3.5,3.5); ax.set_ylim(-1.2,L*1.4+0.8)
    ax.text(0,-1.0,note,ha="center",fontsize=11,style="italic")
def fac(n):
    f={};p=2
    while p*p<=n:
        while n%p==0: f[p]=f.get(p,0)+1; n//=p
        p+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
fig,axs=plt.subplots(1,2,figsize=(12,5))
draw(axs[0],46,"Q47: $\\mathbb{Z}_{47}^{*}\\cong C_{46}$, $46=2\\cdot 23$\n4 subgroups, $I(47)=2$","#27ae60","two incomparable proper subgroups $H_2$, $H_{23}$")
draw(axs[1],40,"Q41: $\\mathbb{Z}_{41}^{*}\\cong C_{40}$, $40=2^3\\cdot 5$\n8 subgroups, $I(41)=6$","#e74c3c","$H_d \\subset H_{d'}$ iff $d\\mid d'$")
plt.tight_layout(); plt.savefig(sys.argv[1]+"/fig2_lattice_v2.png",dpi=200)
