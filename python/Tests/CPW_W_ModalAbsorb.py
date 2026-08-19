"""
 Coplanar waveguide terminated by modal absorbers.

 NOT YET VALIDATED. It runs and converges, but its through-transmission
 disagrees with the reference by ~6 dB and the disagreement is real, not a
 post-processing artefact. See "WHAT MAKES CPW THE HARD ONE" below.

 Rebuilt on the geometry and mesh of CPW_With_WG_Ports.py, which is the
 validated reference and, crucially, the aperture CPW_E.csv / CPW_H.csv were
 generated for. The previous absorber script used port_w_fact = 1.0 and
 port_h_fact = 2.0, giving a 2 x 3 mm window where the mode file describes
 8.5 x 7.0 mm -- the templates were being read through the wrong aperture.

 WHAT MAKES CPW THE HARD ONE
 ---------------------------
  * It is an OPEN structure, so some radiation is expected and |S11| has a floor
    that has nothing to do with the termination.

    BUT DO NOT USE THAT TO EXCUSE THE CURRENT NUMBERS. Measured against
    CPW_With_WG_Ports.py, the validated reference on the same geometry and the
    same mode files, this script does NOT yet agree:

                              reference        this script
        |S21| at 2 GHz          -0.52 dB          -7.06 dB
        |S11|^2+|S21|^2          0.974            0.125-0.277
        raw |U2/U1| at 2 GHz    -0.26 dB          -2.76 dB

    The last row is the one that matters: it is the RAW modal voltage ratio
    between the two port probes, before any de-embedding. The reference loses
    0.26 dB between its ports and this script loses 2.76 dB, so the energy is
    genuinely going somewhere and it is not a post-processing artefact. The
    de-embedding then widens the gap further.

    So the balance sitting at 0.2 is NOT "a CPW radiates". It is unexplained.
    Suspects, in the order worth checking: the reference calls
    FDTD.AddEdges2Grid(dirs='xyz', properties=line), which this does not; the
    reference's ports sit in the airbox straddling the substrate edge at y=0
    and y=substrate_length, where the mode files were generated, whereas these
    sit well inside a substrate that now spans the whole domain; and the metal
    here runs into the transverse MUR boundaries instead of stopping short.
  * There is METAL INSIDE THE APERTURE (centre strip and both grounds), so the
    conductor mask on the deployment template is doing real work here.
  * Quasi-TEM, so no cutoff and the scalar Zw splitter applies.

 WHAT WAS LEARNED ON THE COAX AND IS APPLIED HERE
 ------------------------------------------------
  1. NO OPEN ENDS. The substrate and both metal layers span the ENTIRE y
     domain. The old script left 5 mm of air past each end of a 50 mm line, so
     the structure radiated straight past both sheets and nothing else that
     followed could be measured. That single defect was worth more than every
     absorber detail put together.

  2. PLACEMENT FOLLOWS THE METHOD. The scalar splitter measures direction from
     E = +-Zw*H, so it needs a live E at its own plane and CANNOT sit on PEC;
     it goes SCALAR_INSET cells inside a MUR face, and at one cell it overlaps
     the MUR stencil and diverges. Modal Mur is a one-way overwrite and REQUIRES
     a PEC face, with nothing behind it.

  3. ZERO-MEAN EXCITATION. A soft E source deposits charge proportional to the
     time-integral of its waveform, and inside anything closed by PEC that
     leftover is electrostatic and never leaves. On the coax it pinned a static
     field to the source plane at 16x the line background and stopped the energy
     criterion from ever tripping -- 300000 timesteps versus 31770 with the
     correction on. CPW at 2 +- 1 GHz carries less DC than that case, and the
     run prints how much was actually removed, so the cost is visible.

  4. MEASURE THE INDEX, DO NOT ASSUME IT. beta = 62.5 from a mode solver is the
     ideal; a staircased quasi-TEM line on FR4 is not obliged to agree. The
     E-only modal fit below fits the effective index per frequency and reports
     it. On the coax that gap was 15% and cost several dB.

 (c) 2023-2026 Gadi Lahav <gadi@rfwithcare.com>
"""

import os, tempfile, shutil
import numpy as np
from pylab import *

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS import utilities

# ## Which termination to test
#   'SCALAR' -- Zw splitter, inset inside a MUR face (quasi-TEM: the natural one)
#   'MUR'    -- Dispersive Modal Mur, on a PEC face, declared as the kc -> 0
#               limit of a TM mode (fc = 0)
MODE = 'SCALAR'

# Dump the whole box for inspection in ParaView. Large: budget a few hundred MB.
DUMP_FIELDS = False

Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_CPW_ModalAbsorb_' + MODE)
if not os.path.exists(Sim_Path):
    os.mkdir(Sim_Path)
shutil.copy("CPW_E.csv", Sim_Path)
shutil.copy("CPW_H.csv", Sim_Path)

post_proc_only = False
display_structure = False

# ## Geometry (drawing units = mm) -- from CPW_With_WG_Ports.py
Line_W = 1.0
CPW_gap = 0.3
cu_thick = 0.1
substrate_epsR = 4.3
substrate_width = 11.0
substrate_length = 80.0          # longer than the reference: the fit ladder needs room
substrate_thickness = 1.0
substrate_cells = 4
gap_cells = 3
trace_cells = 7

# The aperture CPW_E.csv / CPW_H.csv were generated for. Do not change these
# without regenerating the mode files.
port_w_fact = 7.5
port_h_fact = 3.5

Airbox_Add = 12.5
unit = 1e-3

f0, fc_exc = 2e9, 1e9

# Modal parameters (same values the reference uses)
beta_ref = 62.5                      # 1/m from a mode solver, at f0
v_phase = 2 * pi * f0 / beta_ref
Zw_mode = 254.0                      # modal wave impedance, Ohm
Zl = 45.0                            # line impedance for the port de-embed

# ## FDTD setup
FDTD = openEMS(NrTS=300000, EndCriteria=1e-5)
FDTD.SetGaussExcite(f0, fc_exc)
# See note 3 in the header. Harmless when the waveform already integrates to
# zero -- it reports the correction it applied, so the cost is never hidden.
FDTD.SetExciteZeroMean(True)

if MODE == 'MUR':
    # PEC on the propagation faces so the one-way sheets terminate them.
    FDTD.SetBoundaryCond(['MUR', 'MUR', 'PEC', 'PEC', 'MUR', 'MUR'])
else:
    # The splitter needs a live E at its plane, so no PEC there.
    FDTD.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)
mesh_res = ((C0 / (f0 + fc_exc)) / unit) / 50

# y is EXACTLY the structure: no air beyond either end (note 1).
SimBox = np.array([
    -substrate_width * 0.5 - Airbox_Add, substrate_width * 0.5 + Airbox_Add,
    0.0, substrate_length,
    -substrate_thickness * (port_h_fact - 1.0) - Airbox_Add,
     substrate_thickness * (1.0 + port_h_fact) + Airbox_Add])
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# centre strip -- spans the whole domain
line = CSX.AddMetal('cu_top')
line.AddBox(priority=20, start=[-Line_W / 2, 0.0, substrate_thickness],
            stop=[Line_W / 2, substrate_length, substrate_thickness + cu_thick])
mesh.AddLine('x', [-Line_W / 2, Line_W / 2])
mesh.AddLine('z', [substrate_thickness, substrate_thickness + cu_thick])
mesh.AddLine('x', linspace(-0.5 * Line_W, 0.5 * Line_W, trace_cells))

# Snap x/z lines to the port edges. The reference notes this matters a lot, and
# it is also what makes the mode file's frame land where it was generated.
p_x = Line_W * (1.0 + port_w_fact) * 0.5
p_z0 = -substrate_thickness * (port_h_fact - 1.0)
p_z1 = substrate_thickness * (1.0 + port_h_fact)
mesh.AddLine('x', [-p_x, p_x])
mesh.AddLine('z', [p_z0, p_z1])

# substrate -- spans the whole domain
sub = CSX.AddMaterial('FR4', epsilon=substrate_epsR)
sub.AddBox(priority=2, start=[-substrate_width / 2, 0.0, 0.0],
           stop=[substrate_width / 2, substrate_length, substrate_thickness])
mesh.AddLine('x', [-substrate_width / 2, substrate_width / 2])
mesh.AddLine('z', [0.0, substrate_thickness])
mesh.AddLine('z', linspace(0, substrate_thickness, substrate_cells + 1))

# ground planes -- span the whole domain
gnd = CSX.AddMetal('cu_gnd')
gnd.AddBox(priority=10, start=[-0.5 * substrate_width, 0.0, substrate_thickness],
           stop=[-0.5 * (Line_W + CPW_gap * 2), substrate_length, substrate_thickness + cu_thick])
gnd.AddBox(priority=10, start=[0.5 * (Line_W + CPW_gap * 2), 0.0, substrate_thickness],
           stop=[0.5 * substrate_width, substrate_length, substrate_thickness + cu_thick])
mesh.AddLine('x', [-0.5 * substrate_width, -0.5 * (Line_W + CPW_gap * 2),
                   0.5 * (Line_W + CPW_gap * 2), 0.5 * substrate_width])
mesh.AddLine('x', linspace(-0.5 * (Line_W + CPW_gap * 2), -Line_W * 0.5, gap_cells))
mesh.AddLine('x', linspace(Line_W * 0.5, 0.5 * (Line_W + CPW_gap * 2), gap_cells))

# UNIFORM y. The delay filter is built from a lattice dispersion relation that
# assumes a locally uniform cell along the propagation direction, and the fit
# ladder below needs integer-cell offsets.
dy = mesh_res / 2.0
mesh.AddLine('y', np.arange(0.0, substrate_length + dy / 2, dy).tolist())
mesh.SmoothMeshLines('x', mesh_res, 1.4)
mesh.SmoothMeshLines('z', mesh_res, 1.4)

Yz = np.asarray(mesh.GetLines('y'))
ny = len(Yz)

# Sheet planes, per method (note 2).
SCALAR_INSET = 4
idxAbs1 = 0 if MODE == 'MUR' else SCALAR_INSET
idxAbs2 = (ny - 1) if MODE == 'MUR' else (ny - 1 - SCALAR_INSET)
idxPort1 = idxAbs1 + 6
idxPort2 = idxAbs2 - 6

box_lo = [-p_x, p_z0]
box_hi = [p_x, p_z1]


def sheet(idx, normal_positive):
    y = Yz.item(int(idx))
    start = [box_lo[0], y, box_lo[1]]
    stop = [box_hi[0], y, box_hi[1]]
    kw = dict(E_file="CPW_E.csv", normal_positive=normal_positive,
              phase_velocity=v_phase, priority=20)
    if MODE == 'MUR':
        # Quasi-TEM is the kc -> 0 limit of a TM mode: fc = 0 leaves the
        # evanescent branch inert and the filter becomes a pure delay.
        kw.update(mode_type='TM', fc=0.0)
    else:
        kw.update(mode_type='TEM', H_file="CPW_H.csv", Zw=Zw_mode)
    return FDTD.AddModalAbsorber(start, stop, 'y', **kw)


abs1 = sheet(idxAbs1, True)

port1 = FDTD.AddWaveGuidePort(1, [box_lo[0], Yz.item(idxPort1), box_lo[1]],
                              [box_hi[0], Yz.item(idxPort1 + 1), box_hi[1]], 'y',
                              E_file="CPW_E.csv", H_file="CPW_H.csv",
                              kc=0.0, excite=1, excite_type=0)
port2 = FDTD.AddWaveGuidePort(2, [box_lo[0], Yz.item(idxPort2), box_lo[1]],
                              [box_hi[0], Yz.item(idxPort2 - 1), box_hi[1]], 'y',
                              E_file="CPW_E.csv", H_file="CPW_H.csv",
                              kc=0.0, excite=0, excite_type=0)

abs2 = sheet(idxAbs2, False)

# ## E-only modal fit ladder (note 4)
#  Fibonacci offsets: a two-plane split is blind wherever beta*d = n*pi, and no
#  in-band frequency is near a null of ALL of these at once.
FIT_OFFSETS = np.array([0, 1, 2, 3, 5, 8, 13, 21, 34])
idxFit0 = idxPort1 + 4
fit_idx = idxFit0 + FIT_OFFSETS
assert fit_idx[-1] < idxPort2 - 1, \
    'fit ladder collides with port 2 (ny=%d); lengthen substrate_length or coarsen dy' % ny
fit_files = []
for n, idx in enumerate(fit_idx):
    yf = Yz.item(int(idx))
    fn = 'mfit_ut_{:02d}'.format(n)
    p = CSX.AddProbe(fn, p_type=10, mode_file_name="CPW_E.csv")
    p.AddBox([box_lo[0], yf, box_lo[1]], [box_hi[0], yf, box_hi[1]])
    fit_files.append(fn)
y_fit_m = (Yz[fit_idx] - Yz.item(int(idxFit0))) * unit

if DUMP_FIELDS:
    Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
    Et.AddBox([SimBox[0], SimBox[2], SimBox[4]], [SimBox[1], SimBox[3], SimBox[5]])

if display_structure:
    CSX_file = os.path.join(Sim_Path, 'cpw_modal_absorb.xml')
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=False)

# ## Post-processing
f = np.linspace(1.0e9, 3.0e9, 201)


def modal_split(z_off_m, Uf_row, gam):
    M = np.stack([np.exp(-gam * z_off_m), np.exp(+gam * z_off_m)], axis=1)
    scale = np.linalg.norm(M, axis=0)
    Mn = M / scale
    sol = np.linalg.lstsq(Mn, Uf_row, rcond=None)[0]
    nrm = np.linalg.norm(Uf_row)
    res = np.linalg.norm(Mn @ sol - Uf_row) / nrm if nrm > 0 else np.inf
    A, B = sol / scale
    return A, B, np.linalg.cond(Mn), res


Uf = np.empty((len(f), len(fit_files)), dtype=complex)
for n, fn in enumerate(fit_files):
    d = np.loadtxt(os.path.join(Sim_Path, fn), comments='%')
    Uf[:, n] = utilities.DFT_time2freq(d[:, 0], d[:, 1], f)

# Fit the effective index per frequency rather than trusting beta_ref.
n_scan = np.linspace(1.10, 2.40, 261)
Gam = np.zeros(len(f), dtype=complex)
n_eff = np.zeros(len(f)); cnd = np.zeros(len(f)); res = np.zeros(len(f))
for i, fi in enumerate(f):
    k0 = 2 * np.pi * fi / C0
    best = None
    for ne in n_scan:
        A, B, c, r = modal_split(y_fit_m, Uf[i], 1j * k0 * ne)
        if (best is None) or (r < best[0]):
            best = (r, ne, A, B, c)
    res[i], n_eff[i], A, B, cnd[i] = best
    Gam[i] = B / A
Gam_dB = 20 * np.log10(np.abs(Gam))

port1.CalcPort(Sim_Path, f, ref_impedance=Zw_mode, ZL=Zl)
port2.CalcPort(Sim_Path, f, ref_impedance=Zw_mode, ZL=Zl)
s11 = np.asarray(port1.uf_ref / port1.uf_inc).ravel()
s21 = np.asarray(port2.uf_ref / port1.uf_inc).ravel()
s11_dB = 20 * np.log10(np.abs(s11))
s21_dB = 20 * np.log10(np.abs(s21))

print('\n===== CPW quasi-TEM, {} termination, no open ends ====='.format(MODE))
print('  mode-solver index at f0: {:.4f}   (beta = {:.1f} 1/m)'.format(
      beta_ref * C0 / (2 * np.pi * f0), beta_ref))
print('  f[GHz]   MODAL-FIT |Gam|   n_eff   resid    port|S11|   |S21| dB')
for fi in [1.25, 1.5, 2.0, 2.5, 2.75]:
    j = np.abs(f - fi * 1e9).argmin()
    print('   {:.2f}       {:9.2f}     {:.4f}  {:.4f}    {:8.2f}   {:8.2f}'.format(
          fi, Gam_dB[j], n_eff[j], res[j], s11_dB[j], s21_dB[j]))

jw = np.argmax(Gam_dB)
print('worst MODAL-FIT |Gamma|: {:.2f} dB at {:.2f} GHz (resid {:.4f})'.format(
      Gam_dB[jw], f[jw] / 1e9, res[jw]))
core = (f > f[0] + 0.1 * (f[-1] - f[0])) & (f < f[-1] - 0.1 * (f[-1] - f[0]))
print('   ... excluding the outer 10% of the band: {:.2f} dB'.format(np.max(Gam_dB[core])))
print('worst port |S11|       : {:.2f} dB'.format(np.max(s11_dB)))
print('modal-fit residual     : max {:.2e}  (two-exponential model)'.format(np.max(res)))
print('   A HIGH residual here is information, not noise: a CPW carries a')
print('   slot-line mode and radiates, and neither fits a single-mode split.')

bal = np.abs(s11) ** 2 + np.abs(s21) ** 2
print('energy balance |S11|^2+|S21|^2: min {:.3f}  max {:.3f}'.format(bal.min(), bal.max()))
print('   Reference CPW_With_WG_Ports.py on this geometry reads 0.974 at 2 GHz')
print('   with |S21| = -0.52 dB. A balance near 0.2 here is NOT radiation --')
print('   see the header. Treat these numbers as not yet validated.')

figure()
plot(f / 1e9, Gam_dB, 'm-', linewidth=2, label=r'$|\Gamma|$ modal fit (E-only)')
plot(f / 1e9, s11_dB, 'k--', linewidth=1.5, label='$S_{11}$ from port')
plot(f / 1e9, s21_dB, 'b-', linewidth=2, label='$S_{21}$ (through)')
grid(); legend(); ylabel('S-Parameter (dB)'); xlabel('Frequency (GHz)')
title('CPW — {} termination, structure spans the whole domain'.format(MODE))

figure()
plot(f / 1e9, n_eff, 'k-', linewidth=2, label='fitted $n_{eff}$')
axhline(beta_ref * C0 / (2 * np.pi * f0), color='r', linestyle='--', label='mode solver at $f_0$')
grid(); legend(); ylabel('effective index'); xlabel('Frequency (GHz)')
title('How far the meshed line sits from the mode-solver index')
show()
