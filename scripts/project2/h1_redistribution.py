"""H1: depth-resolved Ni redistribution. Frozen in PREREG_H1_REDISTRIBUTION.md
(5dda4ac) before this file was written. Ordinal only; no fitting."""
import os, sys, time
import numpy as np, pandas as pd, tifffile
ROOT=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,ROOT)
from cmlib.ground_truth import SAMPLES, ZENODO_LABEL_NOTE
from cmlib.io import label_histogram, slice_paths
from cmlib.percolation import percolating_mask
from cmlib.phases import assign_labels
OUT=os.path.join(ROOT,"out","project2"); AXIS=2; SLAB_UM=0.5

def profile(folder, mapping, vz, vy, vx):
    """Per-slab Phi_Ni, Phi_YSZ, and connected-Ni fraction along axis 2 (x)."""
    ps=slice_paths(folder); a0=tifffile.imread(ps[0])
    nz,ny,nx=len(ps),a0.shape[0],a0.shape[1]
    step=max(1,int(round(SLAB_UM*1000.0/vx)))
    nslab=nx//step
    ni_c=np.zeros(nslab); ys_c=np.zeros(nslab); tot=np.zeros(nslab)
    ni=np.empty((nz,ny,nx),bool); ysz=np.empty((nz,ny,nx),bool)
    for i,p in enumerate(ps):
        a=tifffile.imread(p); ni[i]=a==mapping["Ni"]; ysz[i]=a==mapping["YSZ"]
    conn=percolating_mask(ni,axis=AXIS,connectivity=6)
    cn_c=np.zeros(nslab)
    for s in range(nslab):
        sl=slice(s*step,(s+1)*step)
        ni_c[s]=ni[:,:,sl].sum(); ys_c[s]=ysz[:,:,sl].sum()
        cn_c[s]=conn[:,:,sl].sum(); tot[s]=ni[:,:,sl].size
    return pd.DataFrame(dict(slab=np.arange(nslab),
        x_um=(np.arange(nslab)+0.5)*step*vx/1000.0,
        phi_ni=ni_c/tot, phi_ysz=ys_c/tot, phi_ni_conn=cn_c/tot)), step

rows=[]
for key,folder,grain,state,_nx,_ny,_nz,_vx,_vy,_vz in SAMPLES:
    nz,ny,nx=_nz,_ny,_nx; vz,vy,vx=_vz,_vy,_vx
    mapping=assign_labels(label_histogram(folder)["counts"],ZENODO_LABEL_NOTE[key])
    t=time.time(); df,step=profile(folder,mapping,vz,vy,vx)
    df.insert(0,"state",state); df.insert(0,"anode",grain); df.insert(0,"sample",key)
    rows.append(df)
    n10=max(1,len(df)//10)
    lo=df.phi_ysz.iloc[:n10].mean(); hi=df.phi_ysz.iloc[-n10:].mean()
    grad=abs(df.phi_ysz.iloc[-n10:].mean()-df.phi_ysz.iloc[:n10].mean())
    print(f"{key:12s} slabs={len(df)} step={step}vox({step*vx/1000:.2f}um) "
          f"PhiYSZ first10%={lo:.4f} last10%={hi:.4f} |grad|={grad:.4f} "
          f"{'INTERFACE PRESENT' if grad>=0.05 else 'NO SUSTAINED GRADIENT'} "
          f"PhiNi range {df.phi_ni.min():.4f}-{df.phi_ni.max():.4f} [{time.time()-t:.0f}s]",flush=True)
    pd.concat(rows).to_csv(os.path.join(OUT,"h1_depth_profiles.csv"),index=False)
print("\ndone")
