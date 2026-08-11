"""CSFM: cumulative-strain intergranular fracture of YSZ.
Frozen in out/project2/PREREG_CSFM.md (ffcb066) before implementation.
Frozen, never fitted: eps0=0.0065, sigma=3 vox, r=0.25, tau=0.010.
Deterministic operator -- no damage seed. Ordinal only."""
import os, sys, time
import numpy as np, pandas as pd, tifffile
from scipy import ndimage as ndi
ROOT=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,ROOT)
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE
from cmlib.io import label_histogram, slice_paths
from cmlib.percolation import percolation_summary
from cmlib.phases import assign_labels
from cmlib.pnm import extract_ni_network
from cmlib.roi import tile_rois
OUT=os.path.join(ROOT,"out","project2"); AXIS,CONN=2,6
EPS0,SIGMA,RESID,TAU=0.0065,3.0,0.25,0.010
S6=ndi.generate_binary_structure(3,1)
SIDE={"fine":8.0,"medium":8.0,"coarse":12.0}
REAL_YSZ_RET={"fine":0.958,"medium":0.865,"coarse":0.000}

def load(folder,mapping,r):
    ps=slice_paths(folder); sh=(r["z1"]-r["z0"],r["y1"]-r["y0"],r["x1"]-r["x0"])
    ni=np.empty(sh,bool); ysz=np.empty(sh,bool)
    for i,z in enumerate(range(r["z0"],r["z1"])):
        a=tifffile.imread(ps[z])[r["y0"]:r["y1"],r["x0"]:r["x1"]]
        ni[i]=a==mapping["Ni"]; ysz[i]=a==mapping["YSZ"]
    return ni,ysz

def strain_field(ni):
    """Local Ni volume fraction -> uniform volumetric expansion -> smoothed."""
    return ndi.gaussian_filter(ni.astype(np.float32),sigma=SIGMA)*EPS0

def cum(n): return EPS0*(1.0+RESID*(n-1))/EPS0   # multiplier on eps0

rows=[]
sample={s[2]:s for s in SAMPLES if s[3]=="pristine"}
for grain in ("fine","medium","coarse"):
    k=sample[grain]; nz,ny,nx_=k[6],k[5],k[4]; vz,vy,vx=k[9],k[8],k[7]
    mapping=assign_labels(label_histogram(k[1])["counts"],ZENODO_LABEL_NOTE[k[0]])
    r=tile_rois(nz,ny,nx_,vz,vy,vx,SIDE[grain],max_rois=1)[0]
    t0=time.time(); ni,ysz=load(k[1],mapping,r)
    G,diag,ex=extract_ni_network(ysz,spacing_nm=(vz,vy,vx),axis=AXIS,connectivity=CONN)
    if G is None or G.number_of_edges()==0:
        print(f"{grain}: NO YSZ NETWORK"); continue
    regions=ex["regions"]; slices=ndi.find_objects(regions.astype(np.int32))
    eps=strain_field(ni)
    p0=percolation_summary(ysz,axis=AXIS,connectivity=CONN,check_other_axes=False)
    # per-throat cumulative strain proxy = mean smoothed strain over both grains
    conns=[]; strain=[]; diam=[]
    for u,v,d in G.edges(data=True):
        conns.append((u,v)); diam.append(d["neck_nm"])
        s=0.0; c=0
        for lab in (u+1,v+1):
            sl=slices[lab-1] if lab-1<len(slices) else None
            if sl is None: continue
            m=regions[sl]==lab
            if m.any(): s+=float(eps[sl][m].mean()); c+=1
        strain.append(s/max(c,1))
    conns=np.array(conns); strain=np.array(strain); diam=np.array(diam)
    print(f"{grain:7s} ROI{SIDE[grain]:.0f}um ysz nodes={G.number_of_nodes()} throats={len(conns)} "
          f"pristine P_span={p0['P_span']:.4f} strain[min,med,max]="
          f"[{strain.min():.5f},{np.median(strain):.5f},{strain.max():.5f}] "
          f"[{time.time()-t0:.0f}s]",flush=True)

    def damaged(n):
        broke=strain*cum(n)>TAU
        m=ysz.copy()
        for t in np.flatnonzero(broke):
            a,b=int(conns[t,0])+1,int(conns[t,1])+1
            sa=slices[a-1] if a-1<len(slices) else None
            sb=slices[b-1] if b-1<len(slices) else None
            if sa is None or sb is None: continue
            box=tuple(slice(max(0,min(sa[d].start,sb[d].start)-2),
                            min(regions.shape[d],max(sa[d].stop,sb[d].stop)+2)) for d in range(3))
            sub=regions[box]; ra,rb=sub==a,sub==b
            if not ra.any() or not rb.any(): continue
            iface=(ra&ndi.binary_dilation(rb,S6))|(rb&ndi.binary_dilation(ra,S6))
            m[box]&=~ndi.binary_dilation(iface,S6)
        return m,int(broke.sum())

    for n in (1,2,4,8,12,20):
        m,nb=damaged(n)
        pr=percolation_summary(m,axis=AXIS,connectivity=CONN,check_other_axes=False)
        rows.append(dict(anode=grain,n_cycles=n,n_throats=len(conns),n_broken=nb,
            frac_broken=nb/len(conns),pristine_P_span=p0["P_span"],
            P_span=pr["P_span"],P_largest=pr["P_largest"],n_clusters=pr["n_clusters"],
            retention=pr["P_span"]/p0["P_span"] if p0["P_span"] else np.nan,
            real_retention=REAL_YSZ_RET[grain],
            broken_diam_median=float(np.median(diam[strain*cum(n)>TAU])) if nb else np.nan,
            all_diam_median=float(np.median(diam)),
            pristine_Q=1-p0["P_span"]))
        print(f"   n={n:2d} broken={nb}/{len(conns)} ({nb/len(conns):.3f}) "
              f"P_span={pr['P_span']:.4f} ret={rows[-1]['retention']:.4f}",flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(OUT,"csfm_results.csv"),index=False)
    del ni,ysz,regions,eps
print("\ndone")
