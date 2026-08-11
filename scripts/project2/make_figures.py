"""Manuscript figures. All data from committed CSVs; constants that appear
inline are quoted from committed reports and cited in the caption."""
import os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd
ROOT=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P=os.path.join(ROOT,"out","project2"); F=os.path.join(ROOT,"out","writeup","figs")
os.makedirs(F,exist_ok=True)
plt.rcParams.update({"font.size":9,"axes.grid":True,"grid.alpha":0.25,
                     "figure.dpi":170,"savefig.bbox":"tight"})
C={"fine":"#C44E52","medium":"#DD8452","coarse":"#4C72B0"}
O=["fine","medium","coarse"]

# ---- Fig 1: B6 impossibility ----
n=[0,1,3,5]
neckA6=[63,63,63,63]; neckA26=[63,57,15,0]; neckB=[63,63,63,63]
sA6=[0.45052,0.45195,0.45219,0.45362]; sA26=[0.45052,0.45195,0.44431,0.44216]
sB=[0.45052]*4
fig,ax=plt.subplots(1,2,figsize=(7.2,2.9))
ax[0].plot(n,neckA26,"o-",c="#4C72B0",label="A: curvature-ranked (26-conn)")
ax[0].plot(n,neckA6,"s--",c="#937860",label="A: curvature-ranked (6-conn)")
ax[0].plot(n,neckB,"^-",c="#C44E52",label="B: greedy $\Delta A\leq0$")
ax[0].set_xlabel("damage rounds $n$"); ax[0].set_ylabel("neck volume (voxels)")
ax[0].set_title("(a) neck thinning"); ax[0].legend(fontsize=7,frameon=False)
ax[1].axhline(0.45052,c="k",lw=0.8,ls=":",label="pristine")
ax[1].plot(n,sA26,"o-",c="#4C72B0"); ax[1].plot(n,sA6,"s--",c="#937860")
ax[1].plot(n,sB,"^-",c="#C44E52")
ax[1].fill_between([0,5],0.45052,0.4545,color="#C44E52",alpha=.10)
ax[1].text(2.4,0.4530,"forbidden by gate (ii)",fontsize=7,color="#C44E52")
ax[1].set_xlabel("damage rounds $n$"); ax[1].set_ylabel("$S_{spec}$")
ax[1].set_title("(b) surface area"); ax[1].set_ylim(0.4405,0.4545)
fig.savefig(os.path.join(F,"fig1_b6_impossibility.png")); plt.close(fig)

# ---- Fig 2: lattice vs real min-cut ----
lat=pd.read_csv(os.path.join(P,"audit_ni_vulnerability.csv"))
real=pd.read_csv(os.path.join(P,"c1real_results.csv")).drop_duplicates(["anode","roi"])
fig,ax=plt.subplots(figsize=(4.0,3.0))
for i,a in enumerate(O):
    l=lat[lat.analog==a].frac_mincut; r=real[real.anode==a].mincut_frac.dropna()
    ax.scatter([i-0.13]*len(l),l,c="#937860",s=26,marker="s",
               label="synthetic lattice" if i==0 else None)
    ax.scatter([i+0.13]*len(r),r,c=C[a],s=34,
               label="real ROIs" if i==0 else None)
ax.set_xticks(range(3)); ax.set_xticklabels(O); ax.set_yscale("log")
ax.set_ylabel("min-cut fraction of necks"); ax.legend(fontsize=7,frameon=False)
ax.set_title("Lattice fails at a full cross-section;\nreal networks at 1–2% of throats",fontsize=8.5)
fig.savefig(os.path.join(F,"fig2_mincut.png")); plt.close(fig)

# ---- Fig 3: the reversal ----
d=pd.read_csv(os.path.join(P,"c1real_results.csv"))
m=d[d.threshold==0.50].groupby("anode").transition.mean().loc[O]
realret={"fine":0.6795,"medium":0.8547,"coarse":0.9470}
fig,ax=plt.subplots(1,2,figsize=(6.6,2.8))
ax[0].bar(O,[realret[a] for a in O],color=[C[a] for a in O])
ax[0].set_ylabel("measured Ni retention"); ax[0].set_title("(a) reality: fine worst",fontsize=9)
ax[0].set_ylim(0,1.05)
ax[1].bar(O,[m[a] for a in O],color=[C[a] for a in O])
ax[1].set_ylabel("O6 transition (rounds)"); ax[1].set_ylim(8.5,10.2)
ax[1].set_title("(b) simulation: fine LAST",fontsize=9)
fig.savefig(os.path.join(F,"fig3_reversal.png")); plt.close(fig)

# ---- Fig 4: TPB manufacture ----
g=pd.read_csv(os.path.join(P,"c1real_rni_gate.csv"))
fig,ax=plt.subplots(figsize=(4.0,3.0))
for a in O:
    s=g[g.anode==a].sort_values("n_rounds")
    ax.plot([0]+list(s.n_rounds),[1.0]+list(s.tpb_um2/s.tpb_pristine),
            "o-",c=C[a],label=a)
ax.axhline(1.0,c="k",lw=0.8,ls=":")
ax.set_yscale("log"); ax.set_xlabel("damage rounds $n$")
ax.set_ylabel("TPB / TPB$_{pristine}$")
ax.set_title("Voxel erosion manufactures TPB\nbefore destroying it",fontsize=8.5)
ax.legend(fontsize=7,frameon=False)
fig.savefig(os.path.join(F,"fig4_tpb.png")); plt.close(fig)
print("figures written:", sorted(os.listdir(F)))
