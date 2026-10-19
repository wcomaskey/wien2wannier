#!/usr/bin/env python3
"""Band comparison: a WIEN2k (or other DFT) reference against Wannier model(s) on a
k-path, with the band distances eta / eta_max (Vitale et al., npj Comput.
Mater. 6, 66 (2020)) and window RMS values the wien2wannier documentation
quotes.  Example spec: tools/band_compare_gaas.json.

Rules enforced here:
  1. x is CARTESIAN k-distance, norm((k[i]-k[i-1]) @ B) with B the reciprocal
     lattice -- never fractional distance (wrong for any anisotropic cell),
     never a point index.
  2. Wigner-Seitz correction is the DEFAULT: a _wsvec.dat beside the _hr.dat
     is used; without one the run continues but says WARNING and records
     "UNCORRECTED" in the companion.
  3. ONE E_F convention, DETERMINED from data (band_common.ef_ev/ef_for_eig).
     Both sides are plotted as E - E_F and a gate asserts they agree: the model
     is re-scored under the alternative conventions (+/- E_F) and the run fails
     closed if one of those fits clearly better, or if the model misses the
     reference by more than ef_tol in every convention.
  4. A k-grid shape/alignment mismatch between the DFT reference and a model
     (chunked band run, stale spaghetti, geninterp on another path) fails
     closed naming both files and both shapes; the model is never re-spaced
     onto the DFT axis.
Every figure ships <out>.txt: the spec as used, every input with its mtime,
E_F and its provenance per side, WS status per model, per-panel metrics.

usage: band_compare.py spec.json [--out PNG] [--ef auto|file:P|value:X]
                                 [--offset-search] [--style plain|house]
                                 [--title T] [--ylim LO HI] [--panels 0,2]
"""
import argparse
import datetime
import json
import os
import subprocess
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from band_common import (Ry2eV, hr_load, eig_read, dft_dots,         # noqa: E402
                         check_dots_on_bands, ef_ev, ef_for_eig,
                         band_path, win_kpoints, kread,
                         klist_band_read_auto, _energy_read)

BOHR = 0.529177210903
MINV = 0.5*np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)   # fcc conv -> prim
BCCT = 0.5*np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], float)
COLORS = ['#1465c0', '#c0202a', '#1e6b30', '#a61e1e', '#785ef0', '#fe6100']
LSTYLES = ['-', '-.', '-', '-', '-', '-']
PANEL_KEYS = ('dft', 'cell', 'path', 'ktransform', 'ylim', 'hlines', 'windows',
              'offset_search', 'allow_misfit', 'title_metrics', 'rms_text',
              'eta_mu', 'eta_sigma', 'eta_windows')

SPEC_DOC = r'''
SPEC (JSON).  A dict = one figure; a list of dicts = several figures.
{
 "out": "figs/x.png", "suptitle": "...",
 "root": "<base for relative paths; default: the directory of the spec file>",
 "style": "plain|house", "ylim": [-13, 8],
 "rows": ["panel", "model", [-3, 3]],       (one row per entry; "panel" = each panel's ylim [default],
                                             "model" = the models' own range, or an explicit window)
 "figsize": [6.4, 5.6], "dpi": 180, "labelmap": {"GAMMA": "$\\Gamma$"}, "nseg": 80,
 "legend": "first|each|none", "ef_tol_ev": 0.2, "eta_mu": 2.0, "eta_sigma": 0.1,
 "offset_search": false, "eta_offset": null, "allow_misfit": false, "title_metrics": false,
 "eta_windows": [[-16, 0], [0, 4]],        (window-RMS at the eta offset: index-paired RMS over
                                             states whose DFT energy is inside; Table-I convention)
 "cell": {"nnkp": "seed.nnkp"} | {"win": "seed.win"} | {"acell": [[..],[..],[..]]},
 "ktransform": null | "fcc" | "bcc" | [[..],[..],[..]] | {"lstsq": {"klist": "case.klist", "win": "seed.win"}},
 "dft": {...}, "hlines": [...], "windows": [...],       (panel defaults; a panel overrides per key)
 "panels": [
  {"title": "...",
   "dft": {"kind": "mesh", "energy": "case.energy" | "eig": "seed.eig", "fermi": "case.fermi",
           "kpoints": {"win": "seed.win"} | {"klist": "case.klist"}, "bandrange": [b1, b2]}
       or {"kind": "bandpath", "klist_band": "case.klist_band", "energy": "case.energy",
           "fermi": "case.fermi", "bandrange": [b1, b2], "drop_below": -20, "style": "lines|dots"}
       or {"kind": "spaghetti", "file": "case.spaghetti_ene", "klist_band": "case.klist_band",
           "fermi": "case.fermi"},
   "path": {"nodes": [["$\\Gamma$", [0,0,0]], ["M", [0.5,0,0]]]},          (mesh kind only)
   "models": [
     {"kind": "hr", "dir": "...", "seed": "...", "label": "...", "eig": "seed.eig", "fermi": "...",
      "ef": "auto|relative|absolute|value:<eV>", "color": "#1465c0", "ls": "-", "lw": 1.4,
      "wout": "seed.wout"},
     {"kind": "geninterp", "file": "seed_geninterp.dat", "label": "..."},
     {"kind": "npz", "file": "model.npz", "label": "..."}],
   "hlines": [2.0, {"y": 4.0, "ls": "--", "color": "#c04030", "text": "froz_emax = E_F+2"}],
   "windows": [[-13, 0], [0, 4]], "rms_text": "...", "letter": "a", "title_metrics": true}
 ]
}
"fermi" may be replaced by "ef_ry": <Ry> or "ef_ev": <eV> (logged as a spec constant).
Keys beginning with "_" are ignored.  The E_F line and all energies are E - E_F in eV.
'''

# ---------------------------------------------------------------- plumbing
INPUTS = {}          # path -> (mtime iso, size)


class Gate(SystemExit):
    """fail closed with an actionable message"""
    def __init__(self, msg):
        super().__init__('\nGATE FAILED: ' + msg + '\n')


class Missing(Exception):
    pass


def resolve(p, root):
    if p is None:
        return None
    p = os.path.expanduser(str(p))
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(root, p))


def need(p, what):
    if not p or not os.path.exists(p):
        raise Missing(f'{what}: {p} not found')
    return p


def reg(p):
    st = os.stat(p)
    INPUTS[p] = (datetime.datetime.fromtimestamp(st.st_mtime).isoformat(
        timespec='seconds'), st.st_size)
    return p


def strip_comments(o):
    if isinstance(o, dict):
        return {k: strip_comments(v) for k, v in o.items() if not k.startswith('_')}
    if isinstance(o, list):
        return [strip_comments(v) for v in o]
    return o


# ---------------------------------------------------------------- cell / k
def win_cell(f):
    ls = open(f).read().split('\n')
    for i, L in enumerate(ls):
        if L.strip().lower().startswith('begin unit_cell_cart'):
            j, scale = i + 1, 1.0
            t = ls[j].split()[0]
            if not (t[0].isdigit() or t[0] in '-.'):
                if t.lower().startswith('b'):
                    scale = BOHR
                j += 1
            return scale*np.array([[float(x) for x in ls[j + r].split()[:3]]
                                   for r in range(3)])
    raise Gate(f'no unit_cell_cart block in {f}')


def nnkp_recip(f):
    s = open(f).read().split('begin recip_lattice')[1].split(
        'end recip_lattice')[0].strip().splitlines()
    return np.array([[float(x) for x in r.split()[:3]] for r in s[:3]])


def read_cell(c, root):
    if not c:
        raise Gate('no "cell" given -- the Cartesian k-axis needs the reciprocal '
                   'lattice.  Give {"nnkp": seed.nnkp}, {"win": seed.win} or '
                   '{"acell": [[a1],[a2],[a3]]} (Angstrom rows).')
    if 'acell' in c:
        A, src = np.array(c['acell'], float), 'spec acell'
    elif 'win' in c:
        src = reg(need(resolve(c['win'], root), 'cell.win'))
        A = win_cell(src)
    elif 'nnkp' in c:
        src = reg(need(resolve(c['nnkp'], root), 'cell.nnkp'))
        A = 2*np.pi*np.linalg.inv(nnkp_recip(src)).T
    else:
        raise Gate(f'cell {c}: need one of acell / win / nnkp')
    return A, 2*np.pi*np.linalg.inv(A).T, src


def read_ktransform(t, root):
    if t in (None, 'none', 'identity'):
        return None, 'none (klist k are primitive-fractional already)'
    if t == 'fcc':
        return MINV, 'fcc: conventional klist -> primitive (k @ MINV)'
    if t == 'bcc':
        return BCCT, 'bcc: conventional klist -> primitive (k @ BCCT)'
    if isinstance(t, dict) and 'lstsq' in t:
        kl = reg(need(resolve(t['lstsq']['klist'], root), 'ktransform.lstsq.klist'))
        wf = reg(need(resolve(t['lstsq']['win'], root), 'ktransform.lstsq.win'))
        kc, idv = kread(kl)
        kw = win_kpoints(wf)
        if len(kc) != len(kw):
            raise Gate(f'ktransform.lstsq: {kl} has {len(kc)} k but {wf} has '
                       f'{len(kw)} -- not the same mesh')
        T, *_ = np.linalg.lstsq(kc/idv, kw, rcond=None)
        resid = float(np.abs((kc/idv) @ T - kw).max())
        if resid > 1e-6:
            raise Gate(f'ktransform.lstsq residual {resid:.2e} between {kl} and '
                       f'{wf}: no linear map takes this klist onto the win '
                       f'kpoints.  Either the files belong to different meshes '
                       f'or the win is not this case\'s.')
        return T, f'lstsq {os.path.basename(kl)} -> {os.path.basename(wf)} (residual {resid:.1e})'
    if isinstance(t, list):
        return np.array(t, float), 'spec matrix'
    raise Gate(f'ktransform {t!r} not understood')


def cart_x(ks, B):
    x = np.zeros(len(ks))
    for i in range(1, len(ks)):
        x[i] = x[i-1] + np.linalg.norm((ks[i]-ks[i-1]) @ B)
    return x


# ---------------------------------------------------------------- readers
def energy_blocks(f):
    """band-block .energy* parser: split at band 1."""
    blocks, cur = [], None
    for ln in open(f):
        p = ln.split()
        if len(p) != 2:
            continue
        try:
            bi, e = int(p[0]), float(p[1])
        except ValueError:
            continue
        if bi == 1:
            if cur:
                blocks.append(cur)
            cur = []
        if cur is not None:
            cur.append(e)
    if cur:
        blocks.append(cur)
    return blocks


def energy_rows(f):
    rows = _energy_read(f)
    how = 'band_common._energy_read'
    if not rows:
        rows, how = energy_blocks(f), 'band-block fallback parser'
    if not rows:
        raise Gate(f'{f}: no k-blocks parsed as a WIEN2k .energy* file')
    return rows, how


def rows_to_array(rows, scale=Ry2eV):
    nb = max(len(r) for r in rows)
    E = np.full((nb, len(rows)), np.inf)
    for k, r in enumerate(rows):
        E[:len(r), k] = np.sort(np.array(r, float))*scale
    return E


def spaghetti_read(f):
    """case.spaghetti_ene -> (kdist (nk,), E (nb, nk) eV, already E - E_F)."""
    bands, cur = [], None
    for L in open(f):
        if 'bandindex' in L:
            if cur:
                bands.append(cur)
            cur = []
            continue
        p = L.split()
        if len(p) >= 5 and cur is not None:
            try:
                cur.append((float(p[3]), float(p[4])))
            except ValueError:
                pass
    if cur:
        bands.append(cur)
    if not bands:
        raise Gate(f'{f}: no "bandindex" blocks -- not a spaghetti_ene file')
    n = {len(b) for b in bands}
    if len(n) != 1:
        raise Gate(f'{f}: bands carry different k counts {sorted(n)}')
    xd = np.array([p[0] for p in bands[0]])
    if np.any(np.diff(xd) < -1e-9):
        raise Gate(f'{f}: its k-distance column is not monotone')
    E = np.array([[p[1] for p in b] for b in bands])
    return xd, np.sort(E, axis=0)


def geninterp_read(f):
    """postw90 seed_geninterp.dat -> (kcart (nk,3) 1/A incl 2pi, E (nk,nw) eV)."""
    d = np.loadtxt(f)
    ik = d[:, 0].astype(int)
    nk = ik.max()
    nw = len(d)//nk
    if nk*nw != len(d) or np.any(ik.reshape(nk, nw) != ik.reshape(nk, nw)[:, :1]):
        raise Gate(f'{f}: {len(d)} rows are not nk x nw blocks grouped by k')
    return d[::nw, 1:4], np.sort(d[:, 4].reshape(nk, nw), axis=1)


def omega_total(wout):
    om = None
    for L in open(wout, errors='ignore'):
        if 'Omega Total' in L and '=' in L:
            try:
                om = float(L.split('=')[1].split()[0])
            except (ValueError, IndexError):
                pass
    return om


# ---------------------------------------------------------------- models
def hr_bands(model, ks):
    nw, vecs, ia, ib, wt = model
    flat = ia*nw + ib
    E = np.zeros((len(ks), nw))
    for i, k in enumerate(ks):
        ph = np.exp(2j*np.pi*(vecs @ k))*wt
        H = (np.bincount(flat, ph.real, nw*nw)
             + 1j*np.bincount(flat, ph.imag, nw*nw)).reshape(nw, nw)
        E[i] = np.linalg.eigvalsh((H + H.conj().T)/2)
    return E


def ws_vectors(mesh, Acell):
    """WS cell of the mesh superlattice with ndegen (proj_only_model's)."""
    shift = np.array([[a, b, c] for a in range(-3, 4) for b in range(-3, 4)
                      for c in range(-3, 4)], float)
    i0 = int(np.where((shift == 0).all(axis=1))[0][0])
    sup = np.array(mesh)[:, None]*Acell
    Rs, deg = [], []
    for i in range(-mesh[0], mesh[0]+1):
        for j in range(-mesh[1], mesh[1]+1):
            for l in range(-mesh[2], mesh[2]+1):
                v = np.array([i, j, l], float) @ Acell - shift @ sup
                dd = np.sqrt(np.einsum('ab,ab->a', v, v))
                dmin = dd.min()
                tol = 1e-7*max(1.0, dmin)
                if dd[i0] < dmin + tol:
                    Rs.append([i, j, l]); deg.append(int((dd < dmin + tol).sum()))
    return np.array(Rs, float), np.array(deg, float)


def npz_flat(f):
    """Hk/kp/mesh/Acell/tau grid -> the hr_load flat model, WS-distance
    corrected per orbital pair with tau (the wsvec construction)."""
    z = np.load(f)
    Hk, kp, mesh = z['Hk'], z['kp'], [int(m) for m in z['mesh']]
    Acell, tau = z['Acell'], z['tau']
    nk, nw = Hk.shape[0], Hk.shape[1]
    if len(kp) != nk:
        raise Gate(f'{f}: Hk has {nk} k but kp has {len(kp)}')
    Rs, deg = ws_vectors(mesh, Acell)
    HR = np.tensordot(np.exp(-2j*np.pi*(kp @ Rs.T)).T, Hk, axes=(1, 0))/nk
    sup = np.array(mesh)[:, None]*Acell
    shifts = np.array([[a, b, c] for a in range(-2, 3) for b in range(-2, 3)
                       for c in range(-2, 3)], float)
    S = shifts @ sup
    dtau = tau[None, :, :] - tau[:, None, :]               # [i,j] = tau_j - tau_i
    site_aware = bool(np.any(np.abs(dtau) > 1e-9))
    vecs, ia, ib, wt = [], [], [], []
    for iR, R in enumerate(Rs):
        base = (R + dtau) @ Acell                          # (nw,nw,3)
        dd = np.linalg.norm(base[:, :, None, :] + S[None, None, :, :], axis=-1)
        dmin = dd.min(-1, keepdims=True)
        mask = dd < dmin + 1e-7*np.maximum(1.0, dmin)
        nd = mask.sum(-1)
        i, j, s = np.nonzero(mask)
        vecs.append(R + shifts[s]); ia.append(i); ib.append(j)
        wt.append(HR[iR][i, j]/(deg[iR]*nd[i, j]))
    info = (f'H(R) rebuilt from {nk} mesh k ({mesh[0]}x{mesh[1]}x{mesh[2]}), '
            f'{len(Rs)} WS R-vectors, '
            + ('per-pair WS-distance with tau' if site_aware
               else 'tau all zero -> R-only WS (no site correction)'))
    return (nw, np.concatenate(vecs), np.concatenate(ia), np.concatenate(ib),
            np.concatenate(wt)), info


# ---------------------------------------------------------------- E_F
def dft_fermi(d, args, root):
    """physical E_F in eV for the panel's DFT source, with provenance."""
    if args.ef.startswith('value:'):
        v = float(args.ef[6:])
        print(f'   [E_F] *** CLI OVERRIDE --ef value:{v} eV replaces the spec\'s '
              f'fermi source ***')
        return v, f'CLI OVERRIDE --ef value:{v} (eV)', None
    if args.ef.startswith('file:'):
        f = reg(need(resolve(args.ef[5:], root), '--ef file'))
        v = ef_ev(os.path.dirname(f), os.path.basename(f))
        print(f'   [E_F] *** CLI OVERRIDE --ef file:{f} -> {v:.4f} eV ***')
        return v, f'CLI OVERRIDE --ef file:{f} -> {v:.6f} eV', f
    if args.ef != 'auto':
        raise Gate(f'--ef {args.ef!r}: use auto, file:PATH or value:FLOAT(eV)')
    if d.get('fermi'):
        f = reg(need(resolve(d['fermi'], root), 'dft.fermi'))
        v = ef_ev(os.path.dirname(f), os.path.basename(f))
        return v, f'{f} ({v/Ry2eV:.10f} Ry) -> {v:.6f} eV', f
    if 'ef_ry' in d:
        return float(d['ef_ry'])*Ry2eV, f'spec constant ef_ry={d["ef_ry"]}', None
    if 'ef_ev' in d:
        return float(d['ef_ev']), f'spec constant ef_ev={d["ef_ev"]}', None
    return None, 'none given', None


def model_shift(m, R, root, args):
    mode = str(m.get('ef', 'auto'))
    if mode.startswith('value:'):
        v = float(mode[6:])
        return v, f'spec model.ef value:{v} eV'
    if mode == 'auto' and m.get('eig'):
        eigp = reg(need(resolve(m['eig'], root), 'model.eig'))
        fp = resolve(m.get('fermi'), root) or R['fermi_path']
        if not fp:
            raise Gate(f'model.eig {eigp} needs a fermi file to be checked '
                       f'against: give model.fermi or dft.fermi')
        d = os.path.dirname(eigp)
        s = ef_for_eig(d, os.path.basename(eigp), os.path.relpath(fp, d))
        return s, (f'DETERMINED by ef_for_eig({os.path.basename(eigp)} vs '
                   f'{os.path.basename(fp)}) -> {s:.6f} eV')
    if mode in ('auto', 'relative'):
        return 0.0, ('ASSUMED E_F-relative (the w2w .eig convention); asserted '
                     'by the convention gate')
    if mode == 'absolute':
        if R['EF_phys'] is None:
            raise Gate('model.ef=absolute needs the DFT E_F: give dft.fermi or --ef')
        return R['EF_phys'], f'absolute -> DFT E_F {R["EF_phys"]:.6f} eV'
    raise Gate(f'model.ef {mode!r}: use auto|relative|absolute|value:<eV>')


# ---------------------------------------------------------------- reference
def load_reference(P, args, root, labelmap, nseg):
    d = P.get('dft')
    if not d:
        raise Gate('panel has no "dft" block')
    kind = d.get('kind', 'mesh')
    A, B, cellsrc = read_cell(P.get('cell'), root)
    T, tdesc = read_ktransform(P.get('ktransform'), root)
    R = dict(kind=kind, A=A, B=B, cellsrc=cellsrc, ktransform=tdesc, files=[])
    EF, prov, fpath = dft_fermi(d, args, root)
    R.update(EF_phys=EF, ef_prov=prov, fermi_path=fpath)

    if kind in ('bandpath', 'spaghetti'):
        kb = reg(need(resolve(d.get('klist_band'), root), 'dft.klist_band'))
        kc, labels = klist_band_read_auto(kb)
        if len(kc) == 0:
            raise Gate(f'{kb}: no k-points parsed')
        ks = kc @ T if T is not None else kc
        x = cart_x(ks, B)
        idx = sorted(labels)
        # WIEN2k labels: "GAMMA", or "K.4"-style tags from XCrySDen paths.
        nodelab = [labelmap.get(labels[i], '$\\Gamma$' if labels[i].startswith('GAMMA')
                                else labels[i].split('.')[0]) for i in idx]
        R['files'].append(kb)
    elif kind == 'mesh':
        nodes = (P.get('path') or {}).get('nodes')
        if not nodes:
            raise Gate('dft.kind=mesh needs "path": {"nodes": [[label,[k]],...]}')
        ks, x, nodex = band_path(nodes, A, nseg)
        nodelab = [n[0] for n in nodes]
        idx = None
    else:
        raise Gate(f'dft.kind {kind!r}: use mesh | bandpath | spaghetti')

    if kind == 'bandpath':
        ef = reg(need(resolve(d.get('energy'), root), 'dft.energy'))
        rows, how = energy_rows(ef)
        if len(rows) != len(ks):
            raise Gate(f'dft.bandpath: {len(rows)} k in {ef} but {len(ks)} in {kb} '
                       f'-- a chunked band run concatenated out of order, or one '
                       f'chunk incomplete.  Refusing to plot a misaligned reference.')
        if EF is None:
            raise Gate(f'{ef} holds ABSOLUTE Ry eigenvalues: give dft.fermi '
                       f'(or ef_ry / ef_ev / --ef)')
        Erel = rows_to_array(rows) - EF
        if d.get('drop_below') is not None:
            keep = [b for b in range(Erel.shape[0])
                    if np.mean(Erel[b][np.isfinite(Erel[b])]) > d['drop_below']]
            Erel = Erel[keep]
        R.update(Erel=Erel, EF_sub=EF, parser=how, files=R['files'] + [ef])
        nodex = [x[i] for i in idx]
    elif kind == 'spaghetti':
        sf = reg(need(resolve(d.get('file'), root), 'dft.file'))
        xd, Erel = spaghetti_read(sf)
        if Erel.shape[1] != len(ks):
            raise Gate(f'dft.spaghetti: {sf} has {Erel.shape[1]} k per band but '
                       f'{kb} has {len(ks)} -- stale spaghetti or another path.  '
                       f'Regenerate: x lapw1 -band; x lapwso -band; x spaghetti -so')
        x = xd                                     # the file's own Cartesian column
        nodex = [x[i] for i in idx]
        R.update(Erel=Erel, EF_sub=0.0, parser='spaghetti_ene',
                 files=R['files'] + [sf])
        R['ef_prov'] = ('spaghetti_ene is E - E_F by construction (shift 0); '
                        + ('physical E_F for the gate from ' + prov if EF is not None
                           else 'no fermi given -> convention alternatives not scanned'))
    else:                                          # mesh
        kps = d.get('kpoints') or {}
        if 'win' in kps:
            kf = reg(need(resolve(kps['win'], root), 'dft.kpoints.win'))
            kp = np.array(win_kpoints(kf))
        elif 'klist' in kps:
            kf = reg(need(resolve(kps['klist'], root), 'dft.kpoints.klist'))
            kc, idv = kread(kf)
            kp = kc/idv
            kp = kp @ T if T is not None else kp
        else:
            raise Gate('dft.kind=mesh needs "kpoints": {"win": seed.win} or '
                       '{"klist": case.klist} in the SAME order as the eigenvalue file')
        if d.get('eig'):
            eg = reg(need(resolve(d['eig'], root), 'dft.eig'))
            Eabs = eig_read(os.path.dirname(eg), os.path.basename(eg))
            if not fpath:
                raise Gate(f'dft.eig {eg} needs dft.fermi to determine its reference')
            EFsub = ef_for_eig(os.path.dirname(eg), os.path.basename(eg),
                               os.path.relpath(fpath, os.path.dirname(eg)))
            src, how = eg, f'eig_read + ef_for_eig -> shift {EFsub:.6f} eV'
        else:
            src = reg(need(resolve(d.get('energy'), root), 'dft.energy'))
            rows, how = energy_rows(src)
            Eabs = rows_to_array(rows)
            if EF is None:
                raise Gate(f'{src} holds ABSOLUTE Ry eigenvalues: give dft.fermi')
            EFsub = EF
        if len(kp) != Eabs.shape[1]:
            raise Gate(f'dft.mesh: {kf} lists {len(kp)} k-points but {src} has '
                       f'{Eabs.shape[1]} k columns -- not the same mesh/order; '
                       f'dots would be attached to the wrong k.')
        R.update(Eabs=Eabs, kp=kp, EF_sub=EFsub, parser=how,
                 files=R['files'] + [kf, src])
    R.update(ks=ks, x=x, nodex=nodex, nodelab=nodelab)
    return R


def ref_columns(R):
    """per path index: ALL DFT E-E_F at that k (or None where no reference)."""
    n = len(R['ks'])
    if R['kind'] == 'mesh':
        dx, dy = dft_dots(R['ks'], R['x'], R['kp'], R['Eabs'], R['EF_sub'],
                          (-1e9, 1e9))
        by = {}
        for xv, yv in zip(dx, dy):
            by.setdefault(round(xv, 9), []).append(yv)
        return [np.array(by[round(xv, 9)]) if round(xv, 9) in by else None
                for xv in R['x']]
    return [R['Erel'][:, j][np.isfinite(R['Erel'][:, j])] for j in range(n)]


# ---------------------------------------------------------------- gates+metrics
def misfit(E, cols, delta, lo, hi):
    offs = []
    for j, col in enumerate(cols):
        if col is None or not len(col):
            continue
        e = E[j] - delta
        e = e[(e > lo) & (e < hi)]
        if len(e):
            offs.append(np.abs(e[:, None] - col[None, :]).min(1))
    return float(np.median(np.concatenate(offs))) if offs else float('nan')


def rigid_scan(E, cols, lo, hi, span, step=0.05, maxk=40):
    """best rigid shift of the model onto the reference (diagnostic only;
    subsampled to maxk path points so a large case stays cheap)."""
    idx = [j for j, c in enumerate(cols) if c is not None and len(c)]
    if len(idx) > maxk:
        idx = [idx[i] for i in np.linspace(0, len(idx)-1, maxk).astype(int)]
    keep = set(idx)
    sub = [c if j in keep else None for j, c in enumerate(cols)]
    deltas = np.arange(-span, span + 1e-9, step)
    vals = np.array([misfit(E, sub, d, lo, hi) for d in deltas])
    if not np.isfinite(vals).any():
        return 0.0, float('nan')
    i = int(np.nanargmin(vals))
    return float(deltas[i]), float(vals[i])


def convention_gate(M, cols, R, tol, allow):
    have = [c for c in cols if c is not None and len(c)]
    if not have:
        raise Gate(f'{M["label"]}: no DFT reference point lies on the path -- '
                   f'nothing to compare.  For dft.kind=mesh choose nodes that hit '
                   f'mesh k (the path k mod 1 must equal a kpoints entry to 1e-5).')
    lo = min(c.min() for c in have) - 0.5
    hi = max(c.max() for c in have) + 0.5
    EF = R['EF_phys']
    cands = [(0.0, misfit(M['E'], cols, 0.0, lo, hi))]
    if EF is not None and abs(EF) > 1e-6:
        cands += [(EF, misfit(M['E'], cols, EF, lo, hi)),
                  (-EF, misfit(M['E'], cols, -EF, lo, hi))]
    ranked = sorted(cands, key=lambda t: (np.nan_to_num(t[1], nan=1e9)))
    m0 = cands[0][1]
    dstar, mstar = rigid_scan(M['E'], cols, lo, hi,
                              max(abs(EF) if EF is not None else 0.0, 8.0) + 1.0)
    txt = ', '.join(f'shift {s:+.3f} eV -> median nearest-band miss {v*1000:.1f} meV'
                    for s, v in ranked)
    print(f'   [E_F gate] {M["label"][:40]}: {txt}; best rigid shift '
          f'{dstar:+.2f} eV -> {mstar*1000:.1f} meV')
    reason = None
    best = ranked[0]
    if not np.isfinite(m0):
        reason = 'no model band falls inside the reference energy span'
    elif best[0] != 0.0 and best[1] < m0/3 and (m0 - best[1]) > 0.1:
        reason = (f'an ALTERNATIVE E_F convention fits far better: shifting the '
                  f'model by {best[0]:+.3f} eV drops the median miss from '
                  f'{m0*1000:.0f} to {best[1]*1000:.0f} meV.  The two sides are '
                  f'in different references.')
    elif abs(dstar) > 0.1 and mstar < m0/3 and (m0 - mstar) > 0.1:
        reason = (f'the two sides are RIGIDLY OFFSET by {dstar:+.2f} eV (median '
                  f'miss {m0*1000:.0f} -> {mstar*1000:.0f} meV), which is not '
                  f'+/- the E_F this run was given'
                  + (f' ({EF:.3f} eV)' if EF is not None else ' (none)')
                  + '.  Either that E_F is wrong (--ef / ef_ry / ef_ev / a .fermi '
                  f'from another SCF), or the reference and the model come from '
                  f'different SCF densities or structures (stale band run / '
                  f'spaghetti) -- regenerate the DFT reference from the CURRENT '
                  f'staging.')
    elif m0 > tol:
        reason = (f'the model misses the DFT reference by {m0*1000:.0f} meV median '
                  f'in EVERY convention (ef_tol {tol*1000:.0f} meV; best rigid '
                  f'shift {dstar:+.2f} eV still leaves {mstar*1000:.0f} meV).  '
                  f'Either a model defect or a reference whose bandwidths differ '
                  f'(different SCF) -- an E_F error it is not.')
    verdict = 'PASS' if reason is None else ('ALLOWED: ' + reason if allow
                                             else 'FAIL: ' + reason)
    if reason and not allow:
        alts = '\n'.join(f'      {i+1}. model shift {s:+.4f} eV  (median miss '
                         f'{v*1000:.1f} meV)' for i, (s, v) in enumerate(ranked))
        raise Gate(
            f'E_F convention gate, model "{M["label"]}":\n   {reason}\n'
            f'   DFT E_F: {R["ef_prov"]}\n   model shift: {M["shift_prov"]}\n'
            f'   ranked alternatives:\n{alts}\n'
            f'      best rigid shift {dstar:+.4f} eV  (median miss {mstar*1000:.1f} meV)\n'
            f'   Remedies, in order: (1) give model.eig so ef_for_eig determines '
            f'the model reference from data; (2) set model.ef to "relative" or '
            f'"absolute" if you KNOW the hr\'s reference; (3) if the reference is '
            f'stale, regenerate it (x lapw1 -band [-so]; x spaghetti) from the '
            f'current staging; (4) if this panel is meant to show a broken model, '
            f'set allow_misfit true (or --allow-misfit) -- the number still goes '
            f'to the companion.')
    if reason:
        print(f'   WARNING: {verdict}')
    return dict(verdict=verdict, ranked=ranked, m0=m0, rigid=(dstar, mstar))


def nearest_metrics(E, cols, ylim):
    offs = []
    for j, col in enumerate(cols):
        if col is None or not len(col):
            continue
        e = E[j]
        e = e[(e > ylim[0]) & (e < ylim[1])]
        if len(e):
            offs.append(np.abs(e[:, None] - col[None, :]).min(1))
    if not offs:
        return dict(n=0)
    o = np.concatenate(offs)*1000
    return dict(n=len(o), median=float(np.median(o)), mean=float(o.mean()),
                max=float(o.max()))


def _erfc(x):
    try:
        from scipy.special import erfc
        return erfc(x)
    except ImportError:
        import math
        return np.vectorize(math.erfc)(x)


def eta_at(W, D, off, mu, sigma):
    """Vitale eta / eta_max, model band j <-> DFT band off+j (E-E_F frame)."""
    nw = W.shape[1]
    Dw = D[:, off:off+nw]
    if Dw.shape[1] < nw:
        return None
    ok = np.isfinite(Dw)
    if not ok.any():
        return None
    Dw = np.where(ok, Dw, 0.0)
    w = np.where(np.maximum(Dw, W) > mu - 6*sigma,
                 np.sqrt(0.25*_erfc((Dw-mu)/sigma)*_erfc((W-mu)/sigma)), 1.0)
    w = np.where(ok, w, 0.0)
    if w.sum() <= 0:
        return None
    dev = np.abs(W - Dw)
    eta = float(np.sqrt((w*dev**2).sum()/w.sum()))*1000
    sel = w > 0.5
    emax = float(dev[sel].max())*1000 if sel.any() else float('nan')
    return eta, emax


def window_metrics(W, D, off, windows):
    """Index-paired RMS and max |e_W - e_DFT| over the states whose DFT energy
    lies inside each window (the convention of the paper's band tables), at the
    band offset the Vitale search settled on."""
    nw = W.shape[1]
    Dw = D[:, off:off+nw]
    out = []
    for lo, hi in windows:
        m = np.isfinite(Dw) & (Dw >= lo) & (Dw <= hi)
        if not m.any():
            out.append(dict(window=[lo, hi], n=0))
            continue
        dev = (W - Dw)[m]
        out.append(dict(window=[lo, hi], n=int(m.sum()),
                        rms=float(np.sqrt((dev**2).mean()))*1000,
                        max=float(np.abs(dev).max())*1000))
    return out


def eta_metrics(M, R, mu, sigma, search, offset=None, windows=None):
    """Pairing is by sorted band index, so a windowed model needs the DFT band
    offset: searched (--offset-search) or given.  Neither -> no eta,
    because offset 0 on a model starting at DFT band 6 reports eV of 'error'."""
    if R['kind'] == 'mesh' or (not search and offset is None):
        return None
    W = M['E']
    D = R['Erel'].T                                     # (nk, nb) sorted, inf-padded
    nb = D.shape[1]
    offs = range(0, min(260, nb - W.shape[1] + 1)) if search else [int(offset)]
    best = None
    for off in offs:
        r = eta_at(W, D, off, mu, sigma)
        if r is not None and np.isfinite(r[0]) and (best is None or r[0] < best[0]):
            best = (r[0], r[1], off)
    if best is None:
        return None
    out = dict(eta=best[0], eta_max=best[1], offset=best[2],
               searched=bool(search), mu=mu, sigma=sigma)
    if windows:
        out['windows'] = window_metrics(W, D, best[2], [tuple(w) for w in windows])
    return out


# ---------------------------------------------------------------- panels
def build_panel(P, args, root, spec, labelmap, nseg):
    R = load_reference(P, args, root, labelmap, nseg)
    print(f'   reference: {R["kind"]}, {len(R["ks"])} path k, E_F {R["ef_prov"]}')
    cols = ref_columns(R)
    ylim = tuple(P['ylim'])
    tol = float(args.ef_tol if args.ef_tol is not None else spec.get('ef_tol_ev', 0.2))
    models = []
    for im, m in enumerate(P.get('models') or []):
        kind = m.get('kind', 'hr')
        M = dict(kind=kind, label=m.get('label') or m.get('seed') or
                 os.path.basename(str(m.get('file', m.get('hr', 'model')))),
                 color=m.get('color', COLORS[im % len(COLORS)]),
                 ls=m.get('ls', LSTYLES[im % len(LSTYLES)]), lw=m.get('lw'),
                 files=[])
        t0 = time.time()
        if kind == 'hr':
            if m.get('hr'):
                hrp = resolve(m['hr'], root)
                d, seed = os.path.dirname(hrp), os.path.basename(hrp)[:-7]
            else:
                d, seed = resolve(m.get('dir', '.'), root), m['seed']
                hrp = os.path.join(d, seed + '_hr.dat')
            reg(need(hrp, 'model hr'))
            wsp = os.path.join(d, seed + '_wsvec.dat')
            if os.path.exists(wsp):
                reg(wsp)
                M['ws'] = f'WS-corrected via {wsp}'
                M['files'] += [hrp, wsp]
            else:
                print(f'WARNING: {wsp} not found -- {hrp} is interpolated '
                      f'UNCORRECTED (plain Fourier sum over the hr R-vectors; '
                      f'expect aliasing on coarse meshes).  Rerun wannier90 with '
                      f'use_ws_distance or export_wsvec.py to fix.')
                M['ws'] = f'UNCORRECTED: no {wsp}'
                M['files'] += [hrp]
            model = hr_load(d, seed)
            M['nw'] = model[0]
            E = hr_bands(model, R['ks'])
        elif kind == 'npz':
            f = reg(need(resolve(m.get('file'), root), 'model npz'))
            flat, info = npz_flat(f)
            M.update(nw=flat[0], ws=info, files=[f])
            E = hr_bands(flat, R['ks'])
        elif kind == 'geninterp':
            f = reg(need(resolve(m.get('file'), root), 'model geninterp'))
            kcart, E = geninterp_read(f)
            if E.shape[0] != len(R['ks']):
                raise Gate(f'model.geninterp: {f} holds {E.shape[0]} k but the '
                           f'path ({R["files"][0]}) has {len(R["ks"])} -- postw90 '
                           f'was run on another kpt list.')
            dev = float(np.abs(kcart - R['ks'] @ R['B']).max())
            if dev > 1e-3:
                raise Gate(f'model.geninterp: {f} k-points differ from the path '
                           f'by up to {dev:.3e} 1/A in Cartesian coordinates -- '
                           f'same count, different path or different cell '
                           f'({R["cellsrc"]}).')
            M.update(nw=E.shape[1], files=[f],
                     ws='wannier90 geninterp (WS per its own use_ws_distance)')
        else:
            raise Gate(f'model.kind {kind!r}: use hr | geninterp | npz')
        shift, sprov = model_shift(m, R, root, args)
        M.update(E=E - shift, shift=shift, shift_prov=sprov,
                 seconds=time.time() - t0)
        if m.get('wout'):
            wf = resolve(m['wout'], root)
            if os.path.exists(wf):
                reg(wf)
                M['omega'] = omega_total(wf)
                M['files'].append(wf)
        print(f'   model {M["label"][:40]}: nw={M["nw"]}, {M["ws"]}; shift {sprov}')
        M['gate'] = convention_gate(M, cols, R, tol,
                                    bool(m.get('allow_misfit', P.get('allow_misfit'))
                                         or args.allow_misfit))
        M['nearest'] = nearest_metrics(M['E'], cols, ylim)
        M['eta'] = eta_metrics(M, R, float(P.get('eta_mu', 2.0)),
                               float(P.get('eta_sigma', 0.1)),
                               bool(P.get('offset_search') or args.offset_search),
                               P.get('eta_offset'), P.get('eta_windows'))
        if R['kind'] == 'mesh':
            dx, dy = dft_dots(R['ks'], R['x'], R['kp'], R['Eabs'], R['EF_sub'], ylim)
            M['dot_median'] = check_dots_on_bands(
                R['x'], M['E'], dx, dy, ylim=ylim, tol=0.01,
                label=M['label'][:26].replace('\n', ' ') + ': ')
        nm = M['nearest']
        if nm.get('n'):
            print(f'      nearest-band miss in ylim: median {nm["median"]:.1f}  '
                  f'mean {nm["mean"]:.1f}  max {nm["max"]:.1f} meV  ({nm["n"]} pts)')
        if M['eta']:
            e = M['eta']
            print(f'      eta {e["eta"]:.1f} meV, eta_max {e["eta_max"]:.1f} meV '
                  f'(DFT band offset {e["offset"]}'
                  + (', searched' if e['searched'] else ', given') + ')')
            for wm in e.get('windows', []):
                if wm['n']:
                    print(f'      window [{wm["window"][0]:+.1f}, {wm["window"][1]:+.1f}] eV: '
                          f'RMS {wm["rms"]:.1f} meV, max {wm["max"]:.1f} meV ({wm["n"]} states)')
                else:
                    print(f'      window [{wm["window"][0]:+.1f}, {wm["window"][1]:+.1f}] eV: '
                          f'no DFT states of the paired bands inside')
        elif R['kind'] != 'mesh':
            print('      eta not computed (set offset_search or eta_offset)')
        models.append(M)
    return dict(P=P, R=R, models=models, cols=cols)


# ---------------------------------------------------------------- drawing
def row_ylim(spec_row, D):
    if spec_row == 'panel':
        return tuple(D['P']['ylim'])
    if spec_row == 'model':
        if not D['models']:
            return (-10, 10)
        return (min(M['E'].min() for M in D['models']) - 0.5,
                max(M['E'].max() for M in D['models']) + 0.5)
    return tuple(spec_row)


def draw_reference(ax, R, P, ylim, house, fs=None):
    d = P['dft']
    style = d.get('style', 'dots' if R['kind'] == 'mesh' else 'lines')
    lab = d.get('label', 'DFT')
    if R['kind'] == 'mesh':
        Eabs = R['Eabs']
        br = d.get('bandrange')
        if br:
            Eabs = Eabs[br[0]-1:br[1]]
        dx, dy = dft_dots(R['ks'], R['x'], R['kp'], Eabs, R['EF_sub'], ylim)
        c = fs.C_REF if house else '#111111'
        ax.plot(dx, dy, 'o', ms=2.7, color=c, zorder=5, label=lab)
        return len(dx)
    Erel = np.where(np.isfinite(R['Erel']), R['Erel'], np.nan)
    br = d.get('bandrange')
    n = 0
    for b in range(Erel.shape[0]):
        y = Erel[b]
        if np.all(np.isnan(y)) or np.nanmin(y) > ylim[1] + 1 or np.nanmax(y) < ylim[0] - 1:
            continue
        faint = br and not (br[0] <= b + 1 <= br[1])
        if style == 'dots':
            ax.plot(R['x'], y, 'o', ms=0.7 if faint else 1.0,
                    color='#c9c9c9' if faint else ('#000000' if house else '#111111'),
                    zorder=1 if faint else 5, label=None if (n or faint) else lab)
        elif house:
            ax.plot(R['x'], y, color=fs.C_REF, lw=fs.LW_REF, alpha=0.2 if faint else 0.45,
                    zorder=2, label=None if (n or faint) else lab)
        else:
            ax.plot(R['x'], y, color='#e6e6e6' if faint else '#ccd3da',
                    lw=1.2 if faint else 2.4, zorder=1, solid_capstyle='round',
                    label=None if (n or faint) else lab)
        if not faint:
            n += 1
    return n


def draw_models(ax, D, house, fs=None):
    for im, M in enumerate(D['models']):
        if house:
            c = fs.C_WANNIER if im == 0 else fs.PALETTE[(im - 1) % len(fs.PALETTE)]
            lw = M['lw'] or fs.LW_MAIN
            ax.plot(D['R']['x'], M['E'], color=c, lw=lw, ls=M['ls'], zorder=3)
            ax.plot([], [], color=c, lw=lw, ls=M['ls'], label=M['label'])
        else:
            lw = M['lw'] or (1.4 if im == 0 else 1.5)
            ax.plot(D['R']['x'], M['E'], color=M['color'], lw=lw, ls=M['ls'],
                    zorder=2 + im)
            ax.plot([], [], color=M['color'], lw=lw, ls=M['ls'], label=M['label'])


def draw_hlines(ax, P, x):
    for h in P.get('hlines') or []:
        h = {'y': h} if not isinstance(h, dict) else h
        ax.axhline(h['y'], color=h.get('color', '#c04030'), lw=h.get('lw', 0.8),
                   ls=h.get('ls', '--'), zorder=0)
        if h.get('text'):
            ax.text(x[min(2, len(x)-1)], h['y'] + 0.1, h['text'],
                    color=h.get('color', '#c04030'), fontsize=9, va='bottom')


def panel_title(D):
    P, t = D['P'], D['P'].get('title', '')
    if P.get('title_metrics'):
        bits = []
        for M in D['models']:
            nm = M['nearest']
            if nm.get('n'):
                bits.append(f'{M["label"].split(chr(10))[0]}: max dev '
                            f'{nm["max"]:.1f} meV, mean {nm["mean"]:.1f} meV')
        if bits:
            t += ('\n' if t else '') + ';  '.join(bits)
    return t


def render_plain(spec, data, rows, out, args):
    fw, fh = spec.get('figsize', [6.4, 5.6])
    fig, axes = plt.subplots(len(rows), len(data), figsize=(fw*len(data), fh*len(rows)),
                             sharey='row', squeeze=False)
    legend = spec.get('legend', 'first')
    for col, D in enumerate(data):
        for row, rs in enumerate(rows):
            ax = axes[row][col]
            if D.get('missing'):
                ax.text(.5, .5, D['missing'] + '\nnot available', ha='center',
                        va='center', fontsize=8, transform=ax.transAxes)
                ax.set_title(D.get('title', ''), fontsize=9)
                continue
            R, P = D['R'], D['P']
            ylim = row_ylim(rs, D)
            draw_reference(ax, R, P, ylim, False)
            draw_models(ax, D, False)
            draw_hlines(ax, P, R['x'])
            for xp in R['nodex']:
                ax.axvline(xp, color='0.85', lw=0.7, zorder=0)
            ax.axhline(0, color='0.8', ls='--', lw=0.7, zorder=0)
            ax.set_xticks(R['nodex']); ax.set_xticklabels(R['nodelab'], fontsize=8)
            ax.set_xlim(R['x'][0], R['x'][-1]); ax.set_ylim(*ylim)
            if row == 0:
                ax.set_title(panel_title(D), fontsize=9)
            if col == 0:
                ax.set_ylabel('$E - E_F$ (eV)' + ('   [model range]' if rs == 'model' else ''))
            if legend == 'each' or (legend == 'first' and col == 0 and row == 0):
                ax.legend(fontsize=8, loc='lower right')
    if spec.get('suptitle'):
        fig.suptitle(spec['suptitle'], fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=int(args.dpi or spec.get('dpi', 180)))
    plt.close(fig)
    return out


def dump_curves(path, data):
    """Write the per-panel curves band_compare computed to an npz so an
    external renderer draws the SAME
    reference path and the SAME WS-corrected interpolation (style guide
    section 9: comparisons come from the tool's own interpolation).
    Arrays: p<i>_x, p<i>_Eref (nb x nk, NaN-padded, E-E_F), p<i>_m<j>_E
    (nk x nw); 'meta' = JSON list of per-panel dicts (title, letter, ylim,
    windows, eta_windows, nodex, nodelab, models: label/nw/eta)."""
    import json
    arrays, metas = {}, []
    for i, D in enumerate(data):
        if D.get('missing'):
            metas.append(dict(missing=D['missing'], title=D.get('title', '')))
            continue
        R, P = D['R'], D['P']
        arrays[f'p{i}_x'] = np.asarray(R['x'], float)
        if R['kind'] != 'mesh':
            arrays[f'p{i}_Eref'] = np.asarray(R['Erel'], float)
        ms = []
        for j, M in enumerate(D['models']):
            arrays[f'p{i}_m{j}_E'] = np.asarray(M['E'], float)
            ms.append(dict(label=M['label'], nw=int(M.get('nw', M['E'].shape[1])),
                           eta=M.get('eta'), ws=M.get('ws')))
        metas.append(dict(title=P.get('title', ''), letter=P.get('letter'),
                          ylim=list(P['ylim']), windows=P.get('windows'),
                          eta_windows=P.get('eta_windows'), ref_kind=R['kind'],
                          nodex=[float(v) for v in R['nodex']],
                          nodelab=[str(v) for v in R['nodelab']], models=ms))
    np.savez(path, meta=json.dumps(metas, default=float), **arrays)
    print(f'wrote {path}  ({len(arrays)} arrays, {len(metas)} panels)')


def render_house(spec, data, rows, out, args):
    fsd = spec.get('figstyle_dir')
    if not fsd:
        raise Gate('--style house needs a figstyle.py module: set "figstyle_dir" in the '
                   'spec to the directory that holds it, or use --style plain.')
    sys.path.insert(0, fsd)
    try:
        import figstyle as fs
    except ImportError:
        raise Gate(f'--style house needs figstyle.py; not importable from {fsd}.  '
                   f'Fix "figstyle_dir" in the spec, or use --style plain.')
    fw, fh = spec.get('figsize', [6.4, 6.0])
    fig, axes = fs.new_fig(ncols=len(data), nrows=len(rows), panel=(fw, fh))
    for col, D in enumerate(data):
        for row, rs in enumerate(rows):
            ax = axes[row][col]
            if D.get('missing'):
                ax.text(.5, .5, D['missing'] + '\nnot available', ha='center',
                        va='center', fontsize=10, transform=ax.transAxes)
                continue
            R, P = D['R'], D['P']
            ylim = row_ylim(rs, D)
            draw_reference(ax, R, P, ylim, True, fs)
            draw_models(ax, D, True, fs)
            draw_hlines(ax, P, R['x'])
            fs.band_axes(ax, R['nodelab'], R['nodex'])
            ax.set_ylim(*ylim)
            if P.get('windows'):
                fs.window_bars(ax, [tuple(w) for w in P['windows']])
            txt = P.get('rms_text')
            if txt is None and D['models'] and D['models'][0].get('eta'):
                e = D['models'][0]['eta']
                txt = (rf'$\eta$ = {e["eta"]:.1f} meV   '
                       rf'$\eta_{{\max}}$ = {e["eta_max"]:.0f} meV')
            if txt:
                fs.rms_box(ax, txt, champion=bool(P.get('champion')))
            if row == 0:
                ax.set_title(panel_title(D), fontsize=13)
            if P.get('letter'):
                fs.panel_letter(ax, P['letter'])
            if spec.get('legend', 'none') == 'each' or (
                    spec.get('legend', 'none') == 'first' and col == 0 and row == 0):
                ax.legend(fontsize=9, loc='lower right')
    if spec.get('suptitle'):
        fig.suptitle(spec['suptitle'], fontsize=13)
    stem = os.path.splitext(os.path.basename(out))[0]
    fs.save(fig, stem, outdir=os.path.dirname(out))
    return os.path.join(os.path.dirname(out), stem + '.png')


# ---------------------------------------------------------------- companion
def git_head():
    try:
        return subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=HERE,
                              capture_output=True, text=True, timeout=5).stdout.strip()
    except Exception:
        return 'n/a'


def write_companion(out, spec, specfile, data, args, style):
    txt = os.path.splitext(out)[0] + '.txt'
    L = [f'band_compare.py companion for {out}',
         f'generated {datetime.datetime.now().isoformat(timespec="seconds")}  '
         f'interpreter {sys.executable}  git {git_head()}  style {style}',
         f'spec file: {specfile}',
         f'CLI: {" ".join(sys.argv[1:])}',
         '', '=== spec as used (after CLI overrides) ===',
         json.dumps(spec, indent=1), '', '=== inputs (path, mtime, bytes) ===']
    for p, (mt, sz) in sorted(INPUTS.items()):
        L.append(f'{p}  {mt}  {sz}')
    L += ['', '=== per panel ===']
    for i, D in enumerate(data):
        if D.get('missing'):
            L.append(f'[panel {i}] NOT RENDERED: {D["missing"]}')
            continue
        R, P = D['R'], D['P']
        L.append(f'[panel {i}] {P.get("title", "").replace(chr(10), " | ")}')
        L.append(f'  reference: kind={R["kind"]}  path k={len(R["ks"])}  '
                 f'cell={R["cellsrc"]}  ktransform={R["ktransform"]}  parser={R.get("parser")}')
        L.append(f'  DFT E_F (physical): {R["ef_prov"]}')
        L.append(f'  DFT shift applied: {R["EF_sub"]:.6f} eV')
        for M in D['models']:
            L.append(f'  model "{M["label"].replace(chr(10), " | ")}": kind={M["kind"]} '
                     f'nw={M["nw"]}  files={M["files"]}')
            L.append(f'    WS: {M["ws"]}')
            L.append(f'    model shift: {M["shift"]:.6f} eV  ({M["shift_prov"]})')
            g = M['gate']
            L.append(f'    E_F convention gate: {g["verdict"]}  '
                     + '; '.join(f'shift {s:+.4f} -> {v*1000:.1f} meV' for s, v in g['ranked'])
                     + f'; best rigid shift {g["rigid"][0]:+.2f} eV -> {g["rigid"][1]*1000:.1f} meV')
            nm = M['nearest']
            if nm.get('n'):
                L.append(f'    nearest-band miss inside ylim: median {nm["median"]:.2f} '
                         f'mean {nm["mean"]:.2f} max {nm["max"]:.2f} meV  n={nm["n"]}')
            if M.get('eta'):
                e = M['eta']
                L.append(f'    Vitale eta {e["eta"]:.2f} meV  eta_max {e["eta_max"]:.2f} meV  '
                         f'DFT band offset {e["offset"]} ({"searched" if e["searched"] else "fixed"})'
                         f'  mu=E_F+{e["mu"]} eV sigma={e["sigma"]} eV')
                for wm in e.get('windows', []):
                    if wm['n']:
                        L.append(f'    window-RMS [{wm["window"][0]:+.2f}, {wm["window"][1]:+.2f}] eV: '
                                 f'{wm["rms"]:.2f} meV  max {wm["max"]:.2f} meV  n={wm["n"]}  '
                                 f'(index-paired at offset {e["offset"]}, DFT-energy window)')
                    else:
                        L.append(f'    window-RMS [{wm["window"][0]:+.2f}, {wm["window"][1]:+.2f}] eV: '
                                 f'no paired DFT states inside')
            if 'dot_median' in M:
                L.append(f'    band-to-dot median (mesh invariant): {M["dot_median"]*1000:.4f} meV')
            if M.get('omega') is not None:
                L.append(f'    Omega_total {M["omega"]:.4f} A^2 ({M["omega"]/M["nw"]:.4f} per WF)')
            L.append(f'    evaluation time {M["seconds"]:.1f} s')
    open(txt, 'w').write('\n'.join(L) + '\n')
    return txt


# ---------------------------------------------------------------- driver
def run_figure(spec, args, specfile):
    spec = strip_comments(spec)
    specdir = os.path.dirname(os.path.abspath(specfile))
    root = args.root or spec.get('root') or specdir
    root = resolve(root, specdir)
    if args.out:
        spec['out'] = args.out
    if args.title:
        spec['suptitle'] = args.title
    if args.style:
        spec['style'] = args.style
    if args.ylim:
        spec['ylim'] = list(args.ylim)
        spec['rows'] = [list(args.ylim)]
    if 'out' not in spec:
        raise Gate('spec has no "out" (or give --out)')
    out = resolve(spec['out'], root)
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    style = spec.get('style', 'plain')
    if style not in ('plain', 'house'):
        raise Gate(f'style {style!r}: use plain | house')
    labelmap = spec.get('labelmap', {'GAMMA': '$\\Gamma$'})
    nseg = int(args.nseg or spec.get('nseg', 80))
    panels = spec.get('panels') or []
    if not panels:
        raise Gate('spec has no panels')
    if args.panels:
        sel = [int(s) for s in args.panels.split(',')]
        panels = [panels[i] for i in sel]
        spec['panels'] = panels
    INPUTS.clear()
    print(f'== {out}  ({len(panels)} panel(s), style {style}, root {root})')
    data = []
    for i, P0 in enumerate(panels):
        P = {k: spec[k] for k in PANEL_KEYS if k in spec}
        P.update(P0)
        for k in ('dft', 'cell', 'path'):
            if isinstance(spec.get(k), dict) and isinstance(P0.get(k), dict):
                P[k] = {**spec[k], **P0[k]}
        if args.ylim:
            P['ylim'] = list(args.ylim)
        if 'ylim' not in P:
            raise Gate(f'panel {i}: no ylim (spec-level or panel-level, or --ylim)')
        print(f'-- panel {i}: {P.get("title", "").splitlines()[0] if P.get("title") else ""}')
        try:
            data.append(build_panel(P, args, root, spec, labelmap, nseg))
        except Missing as e:
            if not args.allow_missing:
                raise Gate(f'panel {i} input missing -- {e}\n   Options: (1) stage '
                           f'the file (see the old script this spec replaces for '
                           f'how it was produced); (2) --panels to drop this panel; '
                           f'(3) --allow-missing to render a placeholder, recorded '
                           f'in the companion.')
            print(f'   WARNING: panel {i} skipped: {e}')
            data.append(dict(missing=str(e), title=P.get('title', '')))
    rows = spec.get('rows') or ['panel']
    rows = [r if r in ('model', 'panel') else list(r) for r in rows]
    if all(D.get('missing') for D in data):
        raise Gate('every panel is missing its inputs; nothing to draw')
    if getattr(args, 'dump_npz', None):
        dump_curves(args.dump_npz, data)
    if style == 'house':
        png = render_house(spec, data, rows, out, args)
    else:
        png = render_plain(spec, data, rows, out, args)
    txt = write_companion(png, spec, specfile, data, args, style)
    print(f'wrote {png}\nwrote {txt}')


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split('\n\n')[0] + '\n\n' +
        'Reads a JSON panel spec, evaluates each model on the DFT path and '
        'writes the figure plus its <out>.txt companion.',
        epilog=SPEC_DOC + '\nNeeds Python >= 3.9 with numpy and matplotlib.\n',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('spec', help='JSON panel spec (dict = one figure, list = several)')
    ap.add_argument('--out', help='output PNG (overrides spec "out")')
    ap.add_argument('--ef', default='auto',
                    help='auto (spec fermi sources) | file:PATH (.fermi, Ry) | '
                         'value:FLOAT (eV); an override is logged loudly')
    ap.add_argument('--offset-search', action='store_true',
                    help='search the DFT band offset that minimises Vitale eta '
                         '(for a reference whose energy zero is not trusted)')
    ap.add_argument('--style', choices=['plain', 'house'],
                    help='plain (matplotlib) or house (a figstyle.py module named '
                         'by the spec\'s "figstyle_dir")')
    ap.add_argument('--title', help='figure suptitle override')
    ap.add_argument('--ylim', nargs=2, type=float, metavar=('LO', 'HI'),
                    help='single-row energy window override (eV, E - E_F)')
    ap.add_argument('--panels', help='comma list of panel indices to render')
    ap.add_argument('--root', help='base directory for relative spec paths')
    ap.add_argument('--dump-npz', metavar='FILE',
                    help='also write every panel\'s curves (path x, node positions/'
                         'labels, reference E-E_F, each model\'s interpolated E-E_F, '
                         'eta metrics, spec fields) to FILE for external renderers')
    ap.add_argument('--allow-missing', action='store_true',
                    help='render a placeholder for a panel whose inputs are missing')
    ap.add_argument('--allow-misfit', action='store_true',
                    help='downgrade the E_F convention gate to a WARNING '
                         '(for panels that exist to show a broken model)')
    ap.add_argument('--ef-tol', type=float, metavar='EV',
                    help='E_F gate tolerance on the median nearest-band miss '
                         '(default: spec ef_tol_ev, else 0.2 eV)')
    ap.add_argument('--nseg', type=int, help='points per segment for node paths')
    ap.add_argument('--dpi', type=int)
    args = ap.parse_args()
    spec = json.load(open(args.spec))
    specs = spec if isinstance(spec, list) else [spec]
    if args.out and len(specs) > 1:
        raise Gate('--out with a multi-figure spec: select one with --panels on a '
                   'single-figure spec, or drop --out')
    for s in specs:
        run_figure(s, args, args.spec)


if __name__ == '__main__':
    main()
