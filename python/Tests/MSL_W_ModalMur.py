"""Microstrip line terminated by modal absorbers, 0.1 - 3 GHz.

Replaces the CPW as the quasi-TEM testbed. Propagation is +z, the substrate
normal is +y, and the cross-section is the one from MSL_With_Local_Absorbers.py.

WHAT THIS TEST IS FOR
---------------------
Two terminations on the same structure, chosen with MSL_ABS_MODE:

  'MUR'    -- Dispersive Modal Mur, sheet ON the PEC end face, declared as the
              kc -> 0 limit of a TM mode (fc = 0, phase_velocity carries the
              substrate). read_cells=2, which is worth ~6 dB at the low end.
  'SCALAR' -- the Zw splitter, SCALAR_INSET cells inside the PEC face.

BOTH get PEC propagation faces. That is the point: a MUR face behind a sheet
absorbs on its own, so a sheet in front of one is never measured in isolation.
On the CPW that flattered the splitter by ~17 dB, and the coax and rect tests
since showed the splitter does not terminate at all once the MUR face is gone.

WHAT IS INHERITED FROM THE COAX AND CPW
---------------------------------------
 1. NO OPEN ENDS. Strip, ground and substrate span the entire z domain. A
    structure that stops short of the sheet radiates past it and nothing
    downstream can be measured.
 2. np.unique() on GetLines('z'). CSRectGrid::GetLines defaults to do_sort=False
    and does not deduplicate, so it returns insertion order; every index derived
    from it was silently off on the CPW.
 3. ZERO-MEAN EXCITATION. A soft source deposits charge proportional to the time
    integral of its waveform, and inside a PEC box that residue never leaves.
 4. The E-only modal fit ladder is the ruler, not the port S-parameters: it is a
    two-exponential fit of the mode amplitude along the line and needs no
    reference impedance.

THE MODE FILE DEFINES THE GEOMETRY. MSL_E.csv / MSL_H.csv and
MSL_mode_params.json are written together by the generator, on THIS mesh, for
this cross-section. Regenerate all three or none.
"""

import os, tempfile, shutil, json
from pylab import *
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS import utilities
from openEMS.physical_constants import C0

# ## Which termination to test
MODE = os.environ.get('MSL_ABS_MODE', 'MUR')  # 'MUR' or 'SCALAR'

Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_MSL_ModalMur_' + MODE)
if not os.path.exists(Sim_Path):
    os.mkdir(Sim_Path)

post_proc_only = False
display_structure = False

# Dump the whole box for inspection in ParaView. Large: budget a few hundred MB.
DUMP_FIELDS = bool(int(os.environ.get('MSL_DUMP_FIELDS', '1')))

# ## Geometry (drawing units = mm) -- from MSL_With_Local_Absorbers.py
microstrip_W = 1.875
substrate_epsR = 4.5
substrate_width = 15.0
substrate_length = 120.0  # room for the fit ladder at the low end of the band
substrate_thickness = 1.0
substrate_cells = 8
trace_cells = 7
cu_thick = 0.1

# The port window the mode files are generated for. Do not change these without
# regenerating them.
port_w_fact = 6.0
port_h_fact = 5.5

Airbox_Add = 12.5
unit = 1e-3

f0, fc_exc = 1.55e9, 1.45e9  # 0.1 - 3.0 GHz

# ## FDTD setup
FDTD = openEMS(NrTS=int(os.environ.get('MSL_NRTS', '300000')), EndCriteria=1e-5)
FDTD.SetGaussExcite(f0, fc_exc)
FDTD.SetExciteZeroMean(True)

# x: PEC side walls.  y: PEC under the ground plane, MUR above (the open side).
# z: PEC for BOTH methods -- see the header.
FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'MUR', 'PEC', 'PEC'])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)
mesh_res = ((C0 / (f0 + fc_exc)) / unit) / 75

SimBox = np.array([-substrate_width * 0.5, substrate_width * 0.5,
                   -cu_thick, substrate_thickness * (1.0 + port_h_fact) + Airbox_Add,
                   0.0, substrate_length])  # z is EXACTLY the line
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# microstrip -- spans the whole domain
line = CSX.AddMetal('cu_top')
line.AddBox(priority=20, start=[-microstrip_W / 2, substrate_thickness, 0.0],
            stop=[microstrip_W / 2, substrate_thickness + cu_thick, substrate_length])
mesh.AddLine('x', [-microstrip_W / 2, microstrip_W / 2])
mesh.AddLine('y', [substrate_thickness, substrate_thickness + cu_thick])
mesh.AddLine('x', np.linspace(-0.5 * microstrip_W, 0.5 * microstrip_W, trace_cells))

# Snap lines to the port window edges, so the mode file's frame lands where it
# was generated.
p_x = microstrip_W * (1.0 + port_w_fact) * 0.5
p_y1 = substrate_thickness * (1.0 + port_h_fact)
mesh.AddLine('x', [-p_x, p_x])
mesh.AddLine('y', [0.0, p_y1])

# substrate -- spans the whole domain
sub = CSX.AddMaterial('FR4', epsilon=substrate_epsR)
sub.AddBox(priority=2, start=[-substrate_width / 2, 0.0, 0.0],
           stop=[substrate_width / 2, substrate_thickness, substrate_length])
mesh.AddLine('y', np.linspace(0, substrate_thickness, substrate_cells + 1))

# ground plane -- PEC, not lossy copper: conductor loss and absorption are
# indistinguishable in |Gamma|, and this test measures the absorber.
gnd = CSX.AddMetal('cu_bot')
gnd.AddBox(priority=10, start=[-substrate_width / 2, -cu_thick, 0.0],
           stop=[substrate_width / 2, 0.0, substrate_length])
mesh.AddLine('y', [-cu_thick, 0.0])

# UNIFORM z: the delay filter assumes a locally uniform cell along the
# propagation direction, and the fit ladder needs integer-cell offsets.
# linspace, not arange: arange leaves a sliver last cell (0.083 mm here) that
# drags the timestep down for nothing.
nz_cells = int(round(substrate_length / (mesh_res / 2.0)))
mesh.AddLine('z', np.linspace(0.0, substrate_length, nz_cells + 1).tolist())
dz = substrate_length / nz_cells
mesh.SmoothMeshLines('x', mesh_res, 1.4)
mesh.SmoothMeshLines('y', mesh_res, 1.4)

Zz = np.unique(np.asarray(mesh.GetLines('z')))  # np.unique, see header note 2
nz = len(Zz)
assert np.all(np.diff(Zz) > 0), 'z mesh lines are not strictly increasing'
print('z mesh: %d lines, %.4f .. %.4f mm, dz = %.4f .. %.4f mm'
      % (nz, Zz[0], Zz[-1], np.diff(Zz).min(), np.diff(Zz).max()))

# ---- MESH DONE. Everything below depends on the mode files. -------------
# ## Modal parameters -- written by the SAME generator as the mode files
with open('MSL_mode_params.json') as _fh:
    _mp = json.load(_fh)
assert abs(_mp['microstrip_W'] - microstrip_W) < 1e-9, \
    'MSL_mode_params.json describes a different geometry than this script builds'
beta_ref = _mp['beta']
v_phase = _mp['v_phase']
Zw_mode = _mp['Zw']
Zl = _mp['Zl']
print('mode file: %s, beta = %.4f 1/m (n_eff %.4f), Zw = %.2f Ohm, Zl = %.2f Ohm'
      % (_mp['mode_type'], beta_ref, C0 / v_phase, Zw_mode, Zl))

# Mode files, copied in under the fixed names every consumer below uses, so an
# alternative pair can be swapped in without touching anything else.
E_MODE_SRC = os.environ.get('MSL_E_FILE', 'MSL_E.csv')
H_MODE_SRC = os.environ.get('MSL_H_FILE', 'MSL_H.csv')
shutil.copy(E_MODE_SRC, os.path.join(Sim_Path, 'MSL_E.csv'))
shutil.copy(H_MODE_SRC, os.path.join(Sim_Path, 'MSL_H.csv'))
print('mode files: E <- %s, H <- %s' % (E_MODE_SRC, H_MODE_SRC))

# Sheet planes. MUR sits ON the PEC face; SCALAR needs a live E at its own plane
# so it goes inside, with the PEC face behind it.
SCALAR_INSET = 4
idxAbs1 = 0 if MODE == 'MUR' else SCALAR_INSET
idxAbs2 = (nz - 1) if MODE == 'MUR' else (nz - 1 - SCALAR_INSET)
idxPort1 = idxAbs1 + 8
idxPort2 = idxAbs2 - 8

box_lo = [-p_x, 0.0]  # z-normal port: (x, y) window
box_hi = [p_x, p_y1]


def sheet(idx, normal_positive):
    z = Zz.item(int(idx))
    kw = dict(E_file="MSL_E.csv", normal_positive=normal_positive)
    if MODE == 'MUR':
        # Quasi-TEM: fc = 0 is the pure-delay limit and phase_velocity carries
        # the substrate. fc and v must share one reference -- a mode solver's
        # free-space-referenced (negative) fc paired with v_phase is a 50% beta
        # error at the design frequency.
        kw.update(mode_type='TM', fc=0.0, phase_velocity=v_phase, read_cells=2)
    else:
        kw.update(mode_type='TEM', H_file="MSL_H.csv", Zw=Zw_mode)
    return FDTD.AddModalAbsorber(box_lo + [z], box_hi + [z], 'z', **kw)


# The control that keeps every A/B honest: same mesh, same ports, same
# boundaries, no sheets. Whatever it reports is what the surroundings do alone.
NO_SHEETS = bool(int(os.environ.get('MSL_NO_SHEETS', '0')))

print('  absorber1 @ z[%d]=%.4f   port1 @ z[%d]=%.4f   '
      'port2 @ z[%d]=%.4f   absorber2 @ z[%d]=%.4f'
      % (idxAbs1, Zz[idxAbs1], idxPort1, Zz[idxPort1],
         idxPort2, Zz[idxPort2], idxAbs2, Zz[idxAbs2]))

abs1 = None if NO_SHEETS else sheet(idxAbs1, True)

port1 = FDTD.AddWaveGuidePort(1, box_lo + [Zz.item(idxPort1)], box_hi + [Zz.item(idxPort1 + 1)],
                              'z', E_file="MSL_E.csv", H_file="MSL_H.csv",
                              kc=0.0, excite=1, excite_type=0)
port2 = FDTD.AddWaveGuidePort(2, box_lo + [Zz.item(idxPort2)], box_hi + [Zz.item(idxPort2 - 1)],
                              'z', E_file="MSL_E.csv", H_file="MSL_H.csv",
                              kc=0.0, excite=0, excite_type=0)

abs2 = None if NO_SHEETS else sheet(idxAbs2, False)

# ## E-only modal fit ladder, in the source-free stretch between the ports
FIT_OFFSETS = np.array([0, 1, 2, 3, 5, 8, 13, 21, 34])
idxFit0 = idxPort1 + 4
fit_idx = idxFit0 + FIT_OFFSETS
assert fit_idx[-1] < idxPort2 - 1, 'fit ladder collides with port 2; lengthen substrate_length'
fit_files = []
for n, idx in enumerate(fit_idx):
    zf = Zz.item(int(idx))
    fn = 'mfit_ut_{:02d}'.format(n)
    p = CSX.AddProbe(fn, p_type=10, mode_file_name="MSL_E.csv")
    p.AddBox(box_lo + [zf], box_hi + [zf])
    fit_files.append(fn)
z_fit_m = (Zz[fit_idx] - Zz.item(int(idxFit0))) * unit

# Dump fields
if DUMP_FIELDS:
    Et = CSX.AddDump('Et', file_type=int(os.environ.get('CPW_DUMP_TYPE', '0')), dump_type=0, dump_mode=1)
    Et.AddBox([SimBox[0], SimBox[2], SimBox[4]], [SimBox[1], SimBox[3], SimBox[5]])

if display_structure:
    CSX_file = os.path.join(Sim_Path, 'msl_modal_mur.xml')
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=False)

# ## Post-processing
f = np.linspace(f0 - fc_exc, f0 + fc_exc, 261)


def modal_split(z_off_m, Uf_row, gam):
    """LSQ split into forward/backward at the reference plane, one frequency."""
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

# Fit n_eff per frequency: a microstrip is dispersive, and a wrong beta leaks
# forward into backward.
n_scan = np.linspace(1.30, 2.20, 451)
Gam = np.zeros(len(f), dtype=complex)
n_eff = np.zeros(len(f))
cnd = np.zeros(len(f))
res = np.zeros(len(f))
for i, fi in enumerate(f):
    k0 = 2 * np.pi * fi / C0
    best = None
    for ne in n_scan:
        A, B, c, r = modal_split(z_fit_m, Uf[i], 1j * k0 * ne)
        if (best is None) or (r < best[0]):
            best = (r, ne, A, B, c)
    res[i], n_eff[i], A, B, cnd[i] = best
    Gam[i] = B / A
Gam_dB = 20 * np.log10(np.abs(Gam))

# Conventional port S-parameters, for continuity with the other tests
port1.CalcPort(Sim_Path, f, ref_impedance=np.array([Zw_mode]))
port2.CalcPort(Sim_Path, f, ref_impedance=np.array([Zw_mode]))
s11 = np.asarray(port1.uf_ref / port1.uf_inc).ravel()
s21 = np.asarray(port2.uf_ref / port1.uf_inc).ravel()
s11_dB = 20 * np.log10(np.abs(s11))
s21_dB = 20 * np.log10(np.abs(s21))

print('\n===== MSL quasi-TEM, {} termination, no open ends ====='.format(MODE))
print('  f[GHz]   MODAL-FIT |Gam|   n_eff   resid    port|S11|   |S21| dB')
for fi in [0.10, 0.25, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00]:
    j = np.abs(f - fi * 1e9).argmin()
    print('  {:6.2f}       {:9.2f}     {:.4f}  {:.4f}    {:8.2f}   {:8.2f}'.format(
          f[j] / 1e9, Gam_dB[j], n_eff[j], res[j], s11_dB[j], s21_dB[j]))
jw = int(np.argmax(Gam_dB))
print('worst MODAL-FIT |Gamma|: {:.2f} dB at {:.2f} GHz (resid {:.4f})   <-- the absorber'.format(
      Gam_dB[jw], f[jw] / 1e9, res[jw]))
core = (f > f[0] + 0.1 * (f[-1] - f[0])) & (f < f[-1] - 0.1 * (f[-1] - f[0]))
print('   ... excluding the outer 10% of the band: {:.2f} dB'.format(np.max(Gam_dB[core])))
print('worst port |S11|        in band: {:.2f} dB   <-- de-embed limited'.format(np.max(s11_dB)))
bal = np.abs(s11) ** 2 + np.abs(s21) ** 2
print('energy balance |S11|^2+|S21|^2: min {:.3f}  max {:.3f}'.format(bal.min(), bal.max()))
print('modal-fit residual: max {:.2e}   (two-exponential model)'.format(np.max(res)))

if os.environ.get('MPLBACKEND', '') != 'Agg':
    figure()
    plot(f / 1e9, Gam_dB, 'm-', linewidth=2, label=r'$|\Gamma|$ modal fit (E-only)')
    plot(f / 1e9, s11_dB, 'k--', linewidth=1.5, label='$S_{11}$ from port')
    plot(f / 1e9, s21_dB, 'b-', linewidth=2, label='$S_{21}$ (through)')
    grid(); legend(); ylim(-70, 5)
    ylabel('S-Parameter (dB)'); xlabel('Frequency (GHz)')
    title('MSL quasi-TEM — {} termination'.format(MODE))
    show()
