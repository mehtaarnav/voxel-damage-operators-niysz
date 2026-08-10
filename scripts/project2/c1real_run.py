"""C1-real: Amendment D (corrected-extent graph audit) then the O6 bisection
then the Spearman analysis. Frozen design; no parameter is chosen here.
p_erode=0.35, O6 unchanged, thresholds 0.50 and 0.10, bisection [1,20]."""
import os, sys, time, itertools
import networkx as nx, numpy as np, pandas as pd, tifffile
from scipy import stats
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cmlib.damage2 import apply_o6, tpb_density_um2
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE
from cmlib.io import label_histogram, slice_paths
from cmlib.percolation import percolating_mask, percolation_summary
from cmlib.phases import assign_labels
from cmlib.pnm import extract_ni_network
from cmlib.roi import tile_rois
OUT = os.path.join(ROOT, "out", "project2"); AXIS, CONN = 2, 6
THRESH = (0.50, 0.10); DSEEDS = (300, 301, 302); NROI = 3; INF = 10**9
SIDE = {"fine": 8.0, "medium": 8.0, "coarse": 12.0}

def load(folder, mapping, r):
    ps = slice_paths(folder); sh=(r["z1"]-r["z0"], r["y1"]-r["y0"], r["x1"]-r["x0"])
    ni=np.empty(sh,bool); ysz=np.empty(sh,bool)
    for i,z in enumerate(range(r["z0"], r["z1"])):
        a=tifffile.imread(ps[z])[r["y0"]:r["y1"], r["x0"]:r["x1"]]
        ni[i]=a==mapping["Ni"]; ysz[i]=a==mapping["YSZ"]
    return ni,ysz

def spec_surf(m):
    s=0
    for ax in range(3):
        for sh in (1,-1): s+=int((m & ~np.roll(m,sh,axis=ax)).sum())
    return s/max(int(m.sum()),1)

def mincut_frac(ni,spacing):
    G,diag,ex=extract_ni_network(ni,spacing_nm=spacing,axis=AXIS,connectivity=CONN)
    if G is None or G.number_of_edges()==0: return np.nan,0,0
    H=nx.Graph()
    for i,(u,v,d) in enumerate(G.edges(data=True)): H.add_edge(u,v,idx=i,capacity=1)
    H.add_node("S"); H.add_node("T")
    for n in ex["face_lo"]:
        if n in H: H.add_edge("S",n,idx=-1,capacity=INF)
    for n in ex["face_hi"]:
        if n in H: H.add_edge(n,"T",idx=-1,capacity=INF)
    if "S" not in H or "T" not in H or not nx.has_path(H,"S","T"): return np.nan,G.number_of_nodes(),G.number_of_edges()
    val,(sside,_)=nx.minimum_cut(H,"S","T",capacity="capacity")
    idx={d["idx"] for u,v,d in H.edges(data=True) if d["idx"]>=0 and ((u in sside)!=(v in sside))}
    return len(idx)/G.number_of_edges(), G.number_of_nodes(), G.number_of_edges()

def r_ni(mask,n0): return int(percolating_mask(mask,axis=AXIS,connectivity=CONN).sum())/n0

def bisect(ni,ysz,n0,ds,thr,cache):
    def f(n):
        if n not in cache: cache[n]=apply_o6(ni,ysz,n,ds)[0]
        return r_ni(cache[n],n0)
    lo,hi=1,20
    if f(lo)<thr: return lo-0.5,"floor"
    if f(hi)>=thr: return hi+0.5,"ceiling"
    while hi-lo>1:
        m=(lo+hi)//2
        if f(m)>=thr: lo=m
        else: hi=m
    return 0.5*(lo+hi),"ok"

rows=[]; sample={s[2]:s for s in SAMPLES if s[3]=="pristine"}
for grain in ("fine","medium","coarse"):
    k=sample[grain]; nz,ny,nx_=k[6],k[5],k[4]; vz,vy,vx=k[9],k[8],k[7]
    mapping=assign_labels(label_histogram(k[1])["counts"], ZENODO_LABEL_NOTE[k[0]])
    rois=tile_rois(nz,ny,nx_,vz,vy,vx,SIDE[grain],max_rois=NROI)
    for ri,r in enumerate(rois):
        t0=time.time(); ni,ysz=load(k[1],mapping,r); n0=int(ni.sum())
        vox=float((vz*vy*vx)**(1/3))
        p0=percolation_summary(ni,axis=AXIS,connectivity=CONN,check_other_axes=False)["P_span"]
        ss=spec_surf(ni); tpb0=tpb_density_um2(ni,ysz,vox)
        mcf,nn,ne=mincut_frac(ni,(vz,vy,vx))
        print(f"[D] {grain:7s} roi{ri} P_span={p0:.4f} specSurf={ss:.4f} "
              f"mincut_frac={mcf:.4f} nodes={nn} edges={ne} TPB={tpb0:.3f} "
              f"[{time.time()-t0:.0f}s]",flush=True)
        for ds in DSEEDS:
            cache={}
            for thr in THRESH:
                mid,flag=bisect(ni,ysz,n0,ds,thr,cache)
                ntr=int(np.ceil(mid)); dm=cache.get(ntr)
                if dm is None: dm=apply_o6(ni,ysz,ntr,ds)[0]
                rows.append(dict(anode=grain,roi=ri,damage_seed=ds,threshold=thr,
                    transition=mid,flag=flag,pristine_P_span=p0,spec_surface=ss,
                    mincut_frac=mcf,n_nodes=nn,n_edges=ne,R_Ni_at_tr=r_ni(dm,n0),
                    vol_loss=1-int(dm.sum())/n0,tpb_pristine=tpb0,
                    tpb_at_tr=tpb_density_um2(dm,ysz,vox)))
                print(f"    {grain:7s} roi{ri} d{ds} thr={thr:.2f} mid={mid:.1f}({flag})",flush=True)
                pd.DataFrame(rows).to_csv(os.path.join(OUT,"c1real_results.csv"),index=False)
        del ni,ysz
df=pd.DataFrame(rows); df.to_csv(os.path.join(OUT,"c1real_results.csv"),index=False)
print("\n"+"="*70+"\nC1-REAL")
for thr in THRESH:
    s=df[df.threshold==thr]; m=s.groupby("anode").transition.mean()
    print(f"\nthreshold {thr:.2f}: "+"  ".join(f"{a}={m[a]:.2f}" for a in ("fine","medium","coarse")))
    ok=bool(m["medium"]-m["fine"]>=1.0 and m["coarse"]-m["fine"]>=1.0)
    print(f"  fine vs medium {m['medium']-m['fine']:+.2f}  fine vs coarse {m['coarse']-m['fine']:+.2f}  -> {'PASS' if ok else 'FAIL'}")
    print("  per-ROI means:"); print(s.groupby(["anode","roi"]).transition.mean().round(2).to_string())
print("\nSPEARMAN (pooled over ROIs, mean over damage seeds)")
for thr in THRESH:
    s=df[df.threshold==thr].groupby(["anode","roi"]).mean(numeric_only=True).reset_index()
    for lab,col in (("spec_surface","spec_surface"),("mincut_frac","mincut_frac")):
        v=s.dropna(subset=[col])
        if len(v)<3: print(f"  thr{thr} {lab}: insufficient"); continue
        rho,p=stats.spearmanr(v.transition,v[col])
        # partial, controlling pristine_P_span, via rank residuals
        rk=lambda x: stats.rankdata(x)
        a=rk(v.transition)-np.polyval(np.polyfit(rk(v.pristine_P_span),rk(v.transition),1),rk(v.pristine_P_span))
        b=rk(v[col])-np.polyval(np.polyfit(rk(v.pristine_P_span),rk(v[col]),1),rk(v.pristine_P_span))
        prho,_=stats.pearsonr(a,b)
        loo=[stats.spearmanr(v.drop(v.index[i]).transition,v.drop(v.index[i])[col])[0] for i in range(len(v))]
        print(f"  thr{thr:.2f} {lab:14s} rho={rho:+.3f}  partial={prho:+.3f}  "
              f"LOO[{min(loo):+.3f},{max(loo):+.3f}] signs_consistent={all(np.sign(l)==np.sign(rho) for l in loo)}")
