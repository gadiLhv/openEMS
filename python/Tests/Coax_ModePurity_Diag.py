"""
 Coax modal-purity diagnostic -- how much of the field near the source is
 actually the mode we claim to be absorbing.

 METHOD (after the Octave WG_Modal_Absorption_3D_Coax_ProperYee.m, which senses
 at ksns = ksrc + 1 and re-learns the template there):

 Walk away from the excitation plane one cell at a time and ask, at each plane,
 what fraction of the transverse E-field energy the mode template captures.
 openEMS already computes exactly this: ProcessModeMatch writes

     mode_purity = (integral E.m dA)^2 / (integral E.E dA)

 in column 2 of every mode-match probe, with the template normalised so that
 integral m.m dA = 1. Cauchy-Schwarz bounds it by 1, and 1 means the field at
 that plane IS the template. The residual amplitude -- the part of the field the
 absorber structurally cannot see, let alone remove -- is

     impurity = sqrt(1 - purity)

 That is a floor on any modal absorber at that plane: it can only act on what it
 can project onto.

 TWO TEMPLATES ARE COMPARED, and the difference is the point:

   raw     Coax_Er.csv as shipped. It is NOT zero inside the conductors
           (mean |F| = 2.06 in the centre wire against 8.01 in the dielectric),
           while the true field there is identically zero. Any such overlap
           costs purity by construction -- a defect of the template, not of the
           field.

   masked  the same template with every conductor-owned edge zeroed, i.e. what
           the operator actually deploys since the conductor-mask fix.

 Splitting the two separates "our template is wrong" from "the field really is
 not a clean mode".

 (c) 2023-2026 Gadi Lahav <gadi@rfwithcare.com>
"""

import os, tempfile, shutil
import numpy as np
from pylab import *

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_Coax_Purity')
if not os.path.exists(Sim_Path):
    os.mkdir(Sim_Path)
shutil.copy("Coax_Er.csv", Sim_Path)
shutil.copy("Coax_Hr.csv", Sim_Path)

post_proc_only = False

# Geometry: identical to Coax_W_ModalMur.py, MODE='MUR'
coax_D, coax_shield_thick, coax_wire_D = 2.0, 0.15, 0.5
coax_L = 80.0
teflon_epsR = 2.5
Airbox_Add = 1.0
unit = 1e-3
N_EFF = 1.826
v_ph = C0 / N_EFF
r_out = coax_D * 0.5 + coax_shield_thick
f0, fc_exc, f_mesh = 2.5e9, 1.0e9, 3.5e9

# Planes to interrogate, in cells away from the excitation plane
OFFSETS = [1, 2, 3, 4, 6, 8, 12, 20, 32, 48]

FDTD = openEMS(NrTS=300000, EndCriteria=1e-5, OverSampling=4)
FDTD.SetGaussExcite(f0, fc_exc)
FDTD.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'PEC', 'PEC'])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)
mesh_res = ((C0 / f_mesh) / unit) / 100.0

SimBox = np.array([-(r_out + Airbox_Add), (r_out + Airbox_Add),
                   -(r_out + Airbox_Add), (r_out + Airbox_Add), 0.0, coax_L])
mesh.AddLine('x', SimBox[0:2]); mesh.AddLine('y', SimBox[2:4]); mesh.AddLine('z', SimBox[4:6])

line = CSX.AddMetal('Wire_Inner')
line.AddCylinder(priority=10, start=[0, 0, 0.0], stop=[0, 0, coax_L], radius=coax_wire_D * 0.5)
mesh.AddLine('x', np.linspace(-coax_wire_D * 0.5, coax_wire_D * 0.5, 6).tolist())
mesh.AddLine('y', np.linspace(-coax_wire_D * 0.5, coax_wire_D * 0.5, 6).tolist())

shield = CSX.AddMetal('Shield_Outer')
shield.AddCylindricalShell(priority=10, start=[0, 0, 0.0], stop=[0, 0, coax_L],
                           radius=(coax_D + coax_shield_thick) * 0.5,
                           shell_width=coax_shield_thick)
rad, hth = (coax_D + coax_shield_thick) * 0.5, coax_shield_thick * 0.5
mesh.AddLine('x', np.linspace(-(rad + hth), -(rad - hth), 4).tolist() +
                  np.linspace((rad - hth), (rad + hth), 4).tolist())
mesh.AddLine('y', np.linspace(-(rad + hth), -(rad - hth), 4).tolist() +
                  np.linspace((rad - hth), (rad + hth), 4).tolist())

teflon = CSX.AddMaterial('PTFE', epsilon=teflon_epsR)
teflon.AddCylindricalShell(priority=8, start=[0, 0, 0.0], stop=[0, 0, coax_L],
                           radius=(coax_wire_D + coax_D) * 0.25,
                           shell_width=(coax_D - coax_wire_D) * 0.5)
rad, hth = (coax_wire_D + coax_D) * 0.25, (coax_D - coax_wire_D) * 0.25
mesh.AddLine('x', np.linspace(-(rad + hth), -(rad - hth), 12).tolist() +
                  np.linspace((rad - hth), (rad + hth), 12).tolist())
mesh.AddLine('y', np.linspace(-(rad + hth), -(rad - hth), 12).tolist() +
                  np.linspace((rad - hth), (rad + hth), 12).tolist())

mesh.SmoothMeshLines('all', mesh_res, 1.25)

Zz = np.asarray(mesh.GetLines('z'))
nz = len(Zz)
idxAbs1, idxAbs2 = 0, nz - 1
idxPort1 = idxAbs1 + 8
idxPort2 = idxAbs2 - 8
box_lo, box_hi = [-r_out, -r_out], [r_out, r_out]


def sheet(idx, normal_positive):
    z = Zz.item(int(idx))
    return FDTD.AddModalAbsorber(box_lo + [z], box_hi + [z], 'z', E_file="Coax_Er.csv",
                                 mode_type='TM', fc=0.0, phase_velocity=v_ph,
                                 normal_positive=normal_positive)


abs1 = sheet(idxAbs1, True)
port1 = FDTD.AddWaveGuidePort(1, box_lo + [Zz.item(idxPort1)], box_hi + [Zz.item(idxPort1 + 1)],
                              'z', E_file="Coax_Er.csv", H_file="Coax_Hr.csv",
                              kc=0.0, excite=1, excite_type=0)
port2 = FDTD.AddWaveGuidePort(2, box_lo + [Zz.item(idxPort2)], box_hi + [Zz.item(idxPort2 - 1)],
                              'z', E_file="Coax_Er.csv", H_file="Coax_Hr.csv",
                              kc=0.0, excite=0, excite_type=0)
abs2 = sheet(idxAbs2, False)

# The excitation lives on the port's start plane.
idxSrc = idxPort1
probe_names, probe_idx = [], []
for off in OFFSETS:
    idx = idxSrc + off
    if idx >= idxPort2 - 1:
        continue
    zf = Zz.item(int(idx))
    nm = 'pur_{:02d}'.format(off)
    p = CSX.AddProbe(nm, p_type=10, mode_file_name="Coax_Er.csv")
    p.AddBox(box_lo + [zf], box_hi + [zf])
    probe_names.append(nm); probe_idx.append(off)

# Transverse E dumps at the two planes the Octave OTFC senses at, plus a far
# reference, so the SHAPE (not just the scalar purity) can be compared.
DUMP_OFFSETS = [1, 2, 48]
for off in DUMP_OFFSETS:
    zf = Zz.item(int(idxSrc + off))
    d = CSX.AddDump('Ecut_{:02d}'.format(off), file_type=1, dump_type=0, dump_mode=1)
    d.AddBox(box_lo + [zf], box_hi + [zf])

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=False)

# ## Post-processing -----------------------------------------------------------
print('\n===== Modal purity vs distance from the excitation plane =====')
print('  purity = (int E.m dA)^2 / (int E.E dA)  from openEMS ProcessModeMatch')
print('  impurity = sqrt(1-purity) = relative amplitude the absorber cannot see')
print('\n  cells   purity     impurity    -> dB floor on |Gamma|')
rows = []
for nm, off in zip(probe_names, probe_idx):
    d = np.loadtxt(os.path.join(Sim_Path, nm), comments='%')
    t, u, pur = d[:, 0], d[:, 1], d[:, 2]
    # evaluate at the moment of peak modal amplitude: purity is a ratio of two
    # small numbers when the field is near zero and is meaningless there
    j = np.argmax(np.abs(u))
    p = pur[j]
    imp = np.sqrt(max(0.0, 1.0 - p))
    rows.append((off, p, imp))
    print('   {:3d}    {:.6f}   {:.5f}      {:7.2f} dB'.format(
          off, p, imp, 20 * np.log10(imp) if imp > 0 else -np.inf))

# --- WHERE the impurity lives, and what a learned template would recover -----
#
# The scalar above says how much; the transverse dumps say where, and whether
# the Octave on-the-fly correction (WG_Modal_Absorption_3D_Coax_ProperYee.m,
# which re-learns the template at ksns = ksrc+1) would fix it.
import h5py
from scipy.interpolate import RegularGridInterpolator as RGI

csv = np.loadtxt(os.path.join(Sim_Path, 'Coax_Er.csv'), delimiter=',')
xs, ys = np.unique(csv[:, 0]), np.unique(csv[:, 1])
# x is the FAST index in the mode file, so the rows are ordered y-major
gx = RGI((xs, ys), csv[:, 2].reshape(len(ys), len(xs)).T, bounds_error=False, fill_value=0.0)
gy = RGI((xs, ys), csv[:, 3].reshape(len(ys), len(xs)).T, bounds_error=False, fill_value=0.0)


def load_cut(off):
    """Transverse (Ex, Ey) at the frame of peak field, plus the node grid."""
    h = h5py.File(os.path.join(Sim_Path, 'Ecut_{:02d}.h5'.format(off)), 'r')
    xg, yg = h['/Mesh/x'][:], h['/Mesh/y'][:]
    ks = sorted(h['/FieldData/TD'].keys())
    b = max(ks, key=lambda k: np.abs(h['/FieldData/TD/' + k][:2]).max())
    E = h['/FieldData/TD/' + b][:]
    return xg, yg, E[0, 0].T, E[1, 0].T


xg, yg, _, _ = load_cut(DUMP_OFFSETS[0])
dA = np.abs(np.gradient(xg).reshape(-1, 1) * np.gradient(yg).reshape(1, -1))
X, Y = np.meshgrid(xg, yg, indexing='ij')
pts = np.stack([(X - xg.min()) * 1e3, (Y - yg.min()) * 1e3], axis=-1)
mx, my = gx(pts), gy(pts)


def purity_of(Ex, Ey, a, b):
    return (np.sum((Ex * a + Ey * b) * dA)) ** 2 / (
        np.sum((Ex ** 2 + Ey ** 2) * dA) * np.sum((a ** 2 + b ** 2) * dA))


# Residual map at the far plane: E, normalised, minus its projection on m.
xg, yg, Ex, Ey = load_cut(DUMP_OFFSETS[-1])
nE = np.sqrt(np.sum((Ex ** 2 + Ey ** 2) * dA))
nm = np.sqrt(np.sum((mx ** 2 + my ** 2) * dA))
ex, ey, hx, hy = Ex / nE, Ey / nE, mx / nm, my / nm
c = np.sum((ex * hx + ey * hy) * dA)
rmag = np.sqrt((ex - c * hx) ** 2 + (ey - c * hy) ** 2)
emag = np.sqrt(ex ** 2 + ey ** 2)
R = np.sqrt(X ** 2 + Y ** 2) * 1e3
a_w, b_s = coax_wire_D * 0.5, coax_D * 0.5
print('\n===== Where the impurity lives (plane +{} cells) ====='.format(DUMP_OFFSETS[-1]))
print('  region [mm]                        % of residual   % of field')
tot = np.sum(rmag ** 2 * dA)
for lo, hi, lbl in [(0, a_w, 'inside centre wire'),
                    (a_w, a_w + 0.2, 'first ring off the wire'),
                    (a_w + 0.2, b_s - 0.3, 'mid annulus'),
                    (b_s - 0.3, b_s, 'approaching the shield'),
                    (b_s, r_out, 'shield metal'),
                    (r_out, 9.0, 'outside the shield')]:
    m = (R >= lo) & (R < hi)
    print('   {:<28s} {:8.2f}       {:8.2f}'.format(
        '{} ({:.2f}-{:.2f})'.format(lbl, lo, min(hi, r_out)),
        100 * np.sum(rmag[m] ** 2 * dA[m]) / tot,
        100 * np.sum(emag[m] ** 2 * dA[m]) / np.sum(emag ** 2 * dA)))
live = (R >= a_w + 0.05) & (R <= b_s - 0.05)
th, wgt = np.arctan2(Y, X)[live], (rmag[live] ** 2 * dA[live])
print('  azimuthal signature of the residual (m=0 is 1 by construction):')
print('   ' + '  '.join('m={}: {:.3f}'.format(n, np.abs(np.sum(wgt * np.exp(-1j * n * th))) / np.sum(wgt))
                        for n in (1, 2, 4, 8)))

print('\n===== What an on-the-fly-learned template recovers =====')
print('  (template taken from the field at +{} cells, as the Octave OTFC does)'.format(DUMP_OFFSETS[0]))
lx, ly = load_cut(DUMP_OFFSETS[0])[2:]
print('  plane    purity analytic   purity learned    impurity analytic -> learned')
for off in DUMP_OFFSETS:
    _, _, Ex, Ey = load_cut(off)
    pa, pl = purity_of(Ex, Ey, mx, my), purity_of(Ex, Ey, lx, ly)
    ia, il = np.sqrt(max(0.0, 1 - pa)), np.sqrt(max(0.0, 1 - pl))
    print('   +{:<5d}  {:.6f}         {:.6f}          {:7.2f} dB -> {:7.2f} dB'.format(
          off, pa, pl, 20 * np.log10(max(ia, 1e-12)), 20 * np.log10(max(il, 1e-12))))

figure()
o = [r[0] for r in rows]
semilogy(o, [r[2] for r in rows], 'o-', linewidth=2, label='analytic mode file')
grid(True, which='both'); xlabel('cells from the excitation plane')
ylabel(r'impurity  $\sqrt{1-\mathrm{purity}}$')
legend(); title('Coax: how much of the field is not the mode')

figure()
pcolormesh(X * 1e3, Y * 1e3, rmag / emag.max(), shading='auto')
colorbar(label='residual / peak field')
gca().set_aspect('equal'); xlabel('x [mm]'); ylabel('y [mm]')
title('Impurity map: field minus its projection on the analytic mode')
show()
