"""File readers for band_compare.py: WIEN2k .klist / .klist_band / .energy /
.fermi, Wannier90 _hr.dat (+ _wsvec.dat) / .eig / .win k-points.

hr R-vectors are PRIMITIVE-lattice integers, so a model must be evaluated at
k-points expressed in the primitive reciprocal basis of the .win it came
from (band_compare.py's "ktransform" handles conventional-cell k-lists).
"""
import numpy as np, os

Ry2eV = 13.605693

def kread(f):
    kc=[]; idv=None
    for L in open(f):
        if L.startswith('END'): break
        try:
            x=int(L[10:20]); y=int(L[20:30]); z=int(L[30:40])
            if idv is None: idv=int(L[40:50])
            kc.append((x,y,z))
        except (ValueError,IndexError): break
    return np.array(kc), idv

def klist_band_read(f):
    """WIEN2k .klist_band -> (k_frac (nk,3), labels {index: name}).

    The band path is the ONLY honest DFT reference when the model is compared
    off the mesh: dft_dots can only draw the handful of mesh k that happen to
    lie on the path, and on a coarse mesh (Rh3Si4: 2x10x6) that is a scatter of
    a few columns, not a band structure.

    TRAP the divisor guards against: klist points are integer/divisor, and a
    path whose segments do not divide the divisor evenly gets SNAPPED, making
    the steps non-uniform.  That is invisible on an index x-axis and shows up
    as fake jitter on a Cartesian one.  Callers should pick a divisor that is a
    common multiple of 2*n_seg, and the caller-side
    x-axis must still be Cartesian -- this returns fractional k only.
    """
    ks, labels = [], {}
    for L in open(f):
        if L.startswith('END') or not L.strip():
            break
        try:
            v = [int(L[10 + 10*i:20 + 10*i]) for i in range(4)]
        except (ValueError, IndexError):
            break
        nm = L[:10].strip()
        if nm:
            labels[len(ks)] = nm
        ks.append(np.array(v[:3], float) / v[3])
    if not ks:
        # An empty band path is never legitimate, and the overwhelmingly likely
        # cause is the WRONG COLUMN WIDTH: this reader is (A10,4I10), and many
        # klist_band files are the narrower (A10,4I5) that 'x kgen -band' and
        # XCrySDen write.  Returning [] silently let a caller
        # plot nothing, or index-align a model against an empty reference.
        raise SystemExit(
            f'{f}: parsed 0 k-points as (A10,4I10). If this file is the '
            f'narrower (A10,4I5) layout, call klist_band_read_auto() instead, '
            f'which detects the width. Parsing is otherwise unchanged.')
    return np.array(ks), labels


def klist_band_read_auto(f):
    """klist_band_read for BOTH WIEN2k klist_band layouts.

    XCrySDen/'x kgen -band' write (A10,4I5); the long-divisor files write
    (A10,4I10).  The I10 reader above returns an EMPTY array on an I5 file
    (int('  -70   70') fails on the first line and it breaks out), which a
    caller then plots as nothing or, worse, index-aligns against.  Detect the
    width from the first data line and return the same (k_frac, labels).
    """
    width = None
    for L in open(f):
        if L.startswith('END') or not L.strip():
            break
        for w in (10, 5):
            try:
                v = [int(L[10 + w*i:10 + w*(i+1)]) for i in range(4)]
            except (ValueError, IndexError):
                continue
            if v[3] > 0:
                width = w
                break
        break
    if width is None:
        raise SystemExit(f'{f}: first klist_band line is neither (A10,4I10) '
                         f'nor (A10,4I5)')
    if width == 10:
        return klist_band_read(f)
    ks, labels = [], {}
    for L in open(f):
        if L.startswith('END') or not L.strip():
            break
        try:
            v = [int(L[10 + 5*i:15 + 5*i]) for i in range(4)]
        except (ValueError, IndexError):
            break
        nm = L[:10].strip()
        if nm:
            labels[len(ks)] = nm
        ks.append(np.array(v[:3], float) / v[3])
    return np.array(ks), labels


def hr_load(d,seed):
    hr=open(os.path.join(d,seed+'_hr.dat')).read().split('\n')
    nw=int(hr[1]); nr=int(hr[2]); ndeg=[]; i=3
    while len(ndeg)<nr: ndeg+=[int(x) for x in hr[i].split()]; i+=1
    H={}; RL=[]
    for L in hr[i:]:
        p=L.split()
        if len(p)<7: continue
        R=(int(p[0]),int(p[1]),int(p[2])); a,b=int(p[3])-1,int(p[4])-1
        if R not in H: H[R]=np.zeros((nw,nw),complex); RL.append(R)
        H[R][a,b]=float(p[5])+1j*float(p[6])
    ws={}
    try:
        wl=open(os.path.join(d,seed+'_wsvec.dat')).read().split('\n'); j=1
        while j<len(wl):
            p=wl[j].split()
            if len(p)<5: j+=1; continue
            key=(int(p[0]),int(p[1]),int(p[2]),int(p[3])-1,int(p[4])-1); j+=1
            nd=int(wl[j]); j+=1; vs=[]
            for _ in range(nd):
                q=wl[j].split(); vs.append((int(q[0]),int(q[1]),int(q[2]))); j+=1
            ws[key]=vs
    except OSError: pass
    vecs=[];ia=[];ib=[];wt=[]
    for iR,R in enumerate(RL):
        M=H[R]
        for a in range(nw):
            for b in range(nw):
                if M[a,b]==0: continue
                Ts=ws.get(R+(a,b),[(0,0,0)]) if ws else [(0,0,0)]
                for T in Ts:
                    vecs.append((R[0]+T[0],R[1]+T[1],R[2]+T[2])); ia.append(a); ib.append(b)
                    wt.append(M[a,b]/(ndeg[iR]*len(Ts)))
    return nw,np.array(vecs,float),np.array(ia),np.array(ib),np.array(wt)

# =====================================================================
# BAND-PLOT PRIMITIVES -- one implementation of the E_F convention for
# both sides of a comparison.  Copies that derived it independently once
# plotted the model as E-E_F and left the DFT dots raw, so the dots sat
# E_F away from the bands, which a zero Fermi file hides.  Import these
# rather than re-implementing them.
# =====================================================================

def ef_ev(d, fermifile):
    """Fermi energy in eV from a WIEN2k .fermi (stored in Ry)."""
    return float(open(os.path.join(d, fermifile)).read().split()[0])*Ry2eV


def ef_for_eig(d, eigfile, fermifile, verbose=True):
    """E_F to subtract from a w2w .eig, DETERMINED not assumed.

    w2w writes its .eig ALREADY E_F-RELATIVE, so the answer is normally 0 and
    nothing on the plot should be shifted.  Getting this wrong is invisible:
    shift the model only and the dots float off it; shift both and the curves
    look perfect while the E_F = 0 line -- the one thing a reader takes from
    the figure -- sits eV away from the real Fermi level.  Verified for GaAs
    (w2w.eig), Bi2D (auto36.eig) and alpha-Sn (Wsl.eig): each matches
    (e - E_F)*Ry2eV from the corresponding WIEN2k .energy* to 0.00000 eV.

    When a .energy/.energyso IS staged this cross-checks against it (allowing
    for the w2w window starting at some DFT band > 1) and returns whatever the
    data says.  With no .energy* to check, it returns 0 and says so.
    """
    E = eig_read(d, eigfile)
    EFry = float(open(os.path.join(d, fermifile)).read().split()[0])
    cands = sorted(f for f in os.listdir(d)
                   if f.endswith(('.energyso', '.energysoup', '.energysodn',
                                  '.energy', '.energyup', '.energydn')))
    for c in cands:
        try:
            rows = _energy_read(os.path.join(d, c))
        except Exception:
            continue
        if len(rows) != E.shape[1]:
            continue
        k0 = np.array(rows[0]); tgt = E[:min(8, E.shape[0]), 0]
        n = len(tgt)
        for off in range(0, max(1, len(k0)-n)):
            if np.abs(tgt - (k0[off:off+n]-EFry)*Ry2eV).max() < 1e-3:
                if verbose:
                    print(f'   [E_F] {eigfile} is E_F-RELATIVE (checked against '
                          f'{c}, eig band 1 = DFT band {off+1}) -> shift 0')
                return 0.0
            if np.abs(tgt - k0[off:off+n]*Ry2eV).max() < 1e-3:
                if verbose:
                    print(f'   [E_F] {eigfile} is ABSOLUTE (checked against {c}) '
                          f'-> shift {EFry*Ry2eV:.4f} eV')
                return EFry*Ry2eV
    if verbose:
        print(f'   [E_F] no .energy* staged to cross-check {eigfile}; assuming '
              f'the w2w convention (E_F-relative) -> shift 0')
    return 0.0


def _energy_read(f):
    """WIEN2k .energy* -> list of per-k eigenvalue lists (Ry)."""
    ls = open(f).read().split('\n')
    E, i = [], 0
    while i < len(ls):
        L = ls[i]
        if len(L) >= 79:
            try:
                nen = int(L[73:79]); float(L[0:19]); i += 1
                E.append([float(ls[i+j].split()[1]) for j in range(nen)])
                i += nen; continue
            except (ValueError, IndexError):
                pass
        i += 1
    return E


def eig_read(d, eigfile):
    """w90-format .eig -> (nb, nk) array in eV, in the file's own reference."""
    ed = {}
    for L in open(os.path.join(d, eigfile)):
        p = L.split()
        if len(p) >= 3:
            ed[(int(p[0]), int(p[1]))] = float(p[2])
    nb = max(b for b, _ in ed); nk = max(k for _, k in ed)
    E = np.zeros((nb, nk))
    for (b, k), v in ed.items():
        E[b-1, k-1] = v
    return E


def band_path(nodes, Acell, nseg=80):
    """nodes = [(label, k_prim), ...] -> (ks, x, nodex).

    x is CARTESIAN k-distance, never an index: klist steps are non-uniform
    and an index axis manufactures fake jitter.
    """
    B = 2*np.pi*np.linalg.inv(Acell).T
    ks = []
    for i in range(len(nodes)-1):
        a, b = np.array(nodes[i][1], float), np.array(nodes[i+1][1], float)
        for t in np.linspace(0, 1, nseg, endpoint=(i == len(nodes)-2)):
            ks.append(a + t*(b-a))
    ks = np.array(ks)
    x = np.zeros(len(ks))
    for i in range(1, len(ks)):
        x[i] = x[i-1] + np.linalg.norm((ks[i]-ks[i-1]) @ B)
    nodex = [0.0]
    for i in range(len(nodes)-1):
        nodex.append(nodex[-1] + np.linalg.norm(
            (np.array(nodes[i+1][1], float)-np.array(nodes[i][1], float)) @ B))
    return ks, x, nodex


def dft_dots(ks, x, kp, Eabs, EF, ylim):
    """DFT eigenvalues at the mesh k lying on the path, ALWAYS E_F-shifted.

    Eabs is (nb,nk) in the eig's own reference and EF the matching Fermi
    energy in the SAME reference (0.0 if the eig is already E_F-relative).
    The returned dy is E - EF, so a caller MUST plot its model bands as
    E_model - EF too.  Returning anything else is what caused the bug this
    module exists to prevent.
    """
    kmap = {tuple(np.round(np.mod(k, 1), 5)): i for i, k in enumerate(kp)}
    dx, dy = [], []
    for j, k in enumerate(ks):
        i = kmap.get(tuple(np.round(np.mod(k, 1), 5)))
        if i is None or i >= Eabs.shape[1]:
            continue
        for b in range(Eabs.shape[0]):
            v = Eabs[b, i] - EF
            if ylim[0] < v < ylim[1]:
                dx.append(x[j]); dy.append(v)
    return dx, dy


def check_dots_on_bands(x, Emodel_shifted, dx, dy, tol=0.01, label='', ylim=None):
    """For nw==nb the model equals DFT at every mesh k, so every MODEL BAND
    point must coincide with a DFT eigenvalue.

    The test is directional and it matters which way round.  The model's bands
    are a SUBSET of DFT's, so band->dot is the true invariant; dot->band is
    not, because a panel showing one manifold still plots every DFT dot in the
    energy window, and dots belonging to other bands correctly miss it.  (The
    first version of this check tested dot->band and duly flagged the semicore
    panel of a 4-panel ladder, where 13 eV of "error" was simply the other
    bands' dots.)
    """
    if not dx:
        return float('nan')
    dots = {}
    for xv, yv in zip(dx, dy):
        dots.setdefault(round(xv, 9), []).append(yv)
    off = []
    for j, xv in enumerate(x):
        col = dots.get(round(xv, 9))
        if not col:
            continue
        col = np.asarray(col)
        for e in np.atleast_1d(Emodel_shifted[j]):
            # A model band OUTSIDE the dot window has no dot to match: dft_dots
            # only emits eigenvalues inside ylim.  Scoring those points measures
            # the window, not the model -- MgB2 arms B/C span to +19.9 eV under
            # ylim=[-13,8] and reported 3776 meV while reproducing the eig
            # EXACTLY at every mesh k.  The docstring's subset argument only
            # holds where dots exist, so restrict to there.
            if ylim is not None and not (ylim[0] < e < ylim[1]):
                continue
            off.append(np.abs(col-e).min())
    med = float(np.median(off)) if off else float('nan')
    flag = '' if med < tol else '   *** BANDS DO NOT MEET THE DFT DOTS ***'
    print(f'   {label}band-to-dot median offset: {med*1000:.4f} meV{flag}')
    return med


def win_kpoints(f):
    ks=[]; on=False
    for L in open(f):
        s=L.strip().lower()
        if s.startswith('begin kpoints'): on=True; continue
        if s.startswith('end kpoints'): break
        if on:
            p=L.split()
            if len(p)>=3: ks.append([float(x) for x in p[:3]])
    return np.array(ks)
