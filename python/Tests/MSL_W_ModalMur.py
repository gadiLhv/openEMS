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
from CSXCAD.CSProperties import ABCtype

# ## Which termination to test
MODE = os.environ.get('MSL_ABS_MODE', 'MUR')   # MUR / SCALAR / LOCALMUR / PML
assert MODE in ('MUR', 'SCALAR', 'LOCALMUR', 'PML'), 'unknown MSL_ABS_MODE ' + MODE

# 'PML' is not a sheet at all: the near z face becomes a PML of depth
# MSL_PML_DUT, so openEMS's own absorbing boundary can be measured on the same
# bench, with the same ruler, as the sheets. Only meaningful with MSL_RULER=refsub
# -- referenced against a DEEPER PML at the far end, which is what MSL_PML sets.
PML_DUT = int(os.environ.get('MSL_PML_DUT', '8'))

# ## Modal aperture sweep
# The modal sheet's aperture is not a free parameter -- it is whatever the mode
# file covers -- so enlarging it means re-solving the mode over a bigger
# rectangle. gen_msl_modes_margin.py does that and writes MSL_{E,H}_m<N>.csv
# plus MSL_mode_params_m<N>.json. MSL_MODE_MARGIN selects the set, in the same
# mesh-cell units as MSL_LOCAL_MARGIN, so the modal and local-Mur sweeps line up.
# Empty means the legacy unsuffixed files.
MODE_MARGIN = os.environ.get('MSL_MODE_MARGIN', '')
_msfx = ('_m%d' % int(MODE_MARGIN)) if MODE_MARGIN != '' else ''

# The no-sheet control runs with MODE left at its default, so the flag has to be
# in the path. Without it the control silently overwrites the arm it is meant to
# be compared against.
NO_SHEETS = bool(int(os.environ.get('MSL_NO_SHEETS', '0')))

# ## Which ruler
#
# 'ladder' (default) is the twelve-probe two-exponential fit. It needs the ladder
# to span an appreciable fraction of a guide wavelength, so it goes blind below
# a few hundred MHz on a 120 mm line -- a limit of the METHOD, not of the
# absorber.
#
# 'refsub' is reference subtraction, and has no such limit. It reconfigures the
# bench to a single-ended measurement so that exactly one sheet scatters:
#
#     MSL_REF=1                 z: PML_8 / PML_8, no sheets  -> incident alone
#     MSL_RULER=refsub          z: PEC   / PML_8, sheet 1    -> incident + its echo
#
# Gamma(f) = (U_arm - U_ref) / U_ref at a probe placed between the sheet and the
# port. Subtracted in the FREQUENCY domain, so the two runs need not have equal
# record lengths (they will not; they converge at different rates). No fitting,
# no basis functions, no reference impedance, no dependence on line length.
#
# It is NOT good down to DC, and it is not a free lunch. Its floor is the PML's
# own residual echo, and a PML a few cells deep is a poor absorber when the
# wavelength is metres: measured here, the floor is -29 dB at 3 GHz but only
# -13.7 dB at 0.5 GHz and 0.24 dB at 0.1 GHz. Build two references at different
# depths and the table below marks every row that sits within 6 dB of it.
#
# This is a diagnostic, not a production path: it costs a second run per arm and
# it measures one sheet instead of the pair.
RULER = os.environ.get('MSL_RULER', 'ladder')
assert RULER in ('ladder', 'refsub'), 'unknown MSL_RULER ' + RULER
REF_RUN = bool(int(os.environ.get('MSL_REF', '0')))
PML_N = int(os.environ.get('MSL_PML', '8'))
if REF_RUN:
    NO_SHEETS = True

_suffix = '_nosheet' if (NO_SHEETS and not REF_RUN) else ''
_suffix = _suffix + (('_am%s' % MODE_MARGIN) if MODE_MARGIN != '' else '')
if REF_RUN:
    _suffix = '_ref%d' % PML_N
elif RULER == 'refsub':
    _suffix = '_refsub'

Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_MSL_%s%s' % (MODE, _suffix))
shutil.rmtree(Sim_Path, ignore_errors=True)
os.makedirs(Sim_Path)

post_proc_only = False
display_structure = False

# Dump the whole box for inspection in ParaView. ~1 GB per arm at 300k steps and
# it dominates the runtime, so it is opt-in: MSL_DUMP_FIELDS=1.
DUMP_FIELDS = bool(int(os.environ.get('MSL_DUMP_FIELDS', '0')))

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
# z: PEC for both methods in the default bench -- see the header. The refsub
# ruler opens the far end so only the sheet under test scatters, and its
# reference opens both.
if REF_RUN:
    _zbc = ['PML_%d' % PML_N, 'PML_%d' % PML_N]
elif RULER == 'refsub':
    _zbc = ['PML_%d' % PML_DUT if MODE == 'PML' else 'PEC', 'PML_%d' % PML_N]
else:
    _zbc = ['PEC', 'PEC']
FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'MUR'] + _zbc)
print('boundaries: z = %s / %s   ruler = %s%s'
      % (_zbc[0], _zbc[1], RULER, '  [REFERENCE RUN]' if REF_RUN else ''))

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
with open('MSL_mode_params%s.json' % _msfx) as _fh:
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
E_MODE_SRC = os.environ.get('MSL_E_FILE', 'MSL_E%s.csv' % _msfx)
H_MODE_SRC = os.environ.get('MSL_H_FILE', 'MSL_H%s.csv' % _msfx)
shutil.copy(E_MODE_SRC, os.path.join(Sim_Path, 'MSL_E.csv'))
shutil.copy(H_MODE_SRC, os.path.join(Sim_Path, 'MSL_H.csv'))
print('mode files: E <- %s, H <- %s' % (E_MODE_SRC, H_MODE_SRC))

# Sheet planes. MUR sits ON the PEC face; SCALAR needs a live E at its own plane
# so it goes inside, with the PEC face behind it.
SCALAR_INSET = 4
# PML shares LOCALMUR's inset so the refsub probe lands on the same z line and
# the two are directly comparable.
INSET = {'MUR': 0, 'LOCALMUR': 1, 'SCALAR': SCALAR_INSET, 'PML': 1}[MODE]
idxAbs1 = INSET
idxAbs2 = nz - 1 - INSET
idxPort1 = idxAbs1 + 8
idxPort2 = idxAbs2 - 8

# z-normal port: the (x, y) window. It MUST be the rectangle the mode file was
# solved over, or the sheet writes a mode into cells it was never computed for.
if _msfx:
    box_lo = [_mp['win_x0'], _mp['win_y0']]
    box_hi = [_mp['win_x1'], _mp['win_y1']]
else:
    box_lo = [-p_x, 0.0]
    box_hi = [p_x, p_y1]
print('  modal aperture: x %+.4f..%+.4f  y %+.4f..%+.4f mm  (%s)'
      % (box_lo[0], box_hi[0], box_lo[1], box_hi[1],
         ('margin %s' % MODE_MARGIN) if _msfx else 'legacy window'))


def sheet(idx, normal_positive):
    z = Zz.item(int(idx))
    if MODE == 'LOCALMUR':
        # The plain local Mur sheet -- openEMS's own ABCtype.MUR_1ST, no mode
        # file, no projection. One cell in from the PEC face, which is what
        # MSL_With_Local_Absorbers.py recommends: on the face it leaks less but
        # absorbs ~3 dB worse. MUR_1ST_SA is the super-absorbing variant and is
        # noted there as unstable on this geometry.
        #
        # It needs no mode, so it is NOT confined to the port window: it spans
        # the whole cross-section, including the ground plane and the airbox
        # above the port height. MSL_LOCAL_WINDOW=1 shrinks it to the modal
        # arms' window instead, which is the like-for-like comparison.
        # MSL_LOCAL_MARGIN sets the aperture: 'full' is the whole cross-section,
        # an integer is that many mesh cells of margin around the port window on
        # every side, clipped to the SimBox. 0 is the port window exactly, which
        # is the aperture the modal sheets are locked to by their mode file.
        # This is the knob that answers "how big does it have to be".
        _m = os.environ.get('MSL_LOCAL_MARGIN', 'full')
        if bool(int(os.environ.get('MSL_LOCAL_WINDOW', '0'))):
            _m = '0'
        if _m == 'full':
            lo = [SimBox[0], SimBox[2]]
            hi = [SimBox[1], SimBox[3]]
        else:
            _m = int(_m)
            _xs = np.unique(np.asarray(mesh.GetLines('x')))
            _ys = np.unique(np.asarray(mesh.GetLines('y')))
            _ix0 = int(np.argmin(np.abs(_xs + p_x))) - _m
            _ix1 = int(np.argmin(np.abs(_xs - p_x))) + _m
            _iy0 = int(np.argmin(np.abs(_ys - 0.0))) - _m
            _iy1 = int(np.argmin(np.abs(_ys - p_y1))) + _m
            _ix0 = max(0, _ix0); _ix1 = min(len(_xs) - 1, _ix1)
            _iy0 = max(0, _iy0); _iy1 = min(len(_ys) - 1, _iy1)
            lo = [float(_xs[_ix0]), float(_ys[_iy0])]
            hi = [float(_xs[_ix1]), float(_ys[_iy1])]
        print('  local Mur aperture: x %+.4f..%+.4f  y %+.4f..%+.4f mm   '
              '(window is x %+.4f..%+.4f  y %+.4f..%+.4f; box is x %+.4f..%+.4f  y %+.4f..%+.4f)'
              % (lo[0], hi[0], lo[1], hi[1], box_lo[0], box_hi[0], box_lo[1], box_hi[1],
                 SimBox[0], SimBox[1], SimBox[2], SimBox[3]))
        a = CSX.AddAbsorbingBC('abs_%d' % idx, NormalSignPositive=normal_positive,
                               AbsorbingBoundaryType=ABCtype.MUR_1ST)
        a.AddBox(lo + [z], hi + [z], priority=20)
        a.SetPhaseVelocity(v_phase)
        return a
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


# NO_SHEETS (read at the top, because Sim_Path depends on it) is the control
# that keeps every A/B honest: same mesh, same ports, same boundaries, no
# sheets. Whatever it reports is what the surroundings do alone.
print('  absorber1 @ z[%d]=%.4f   port1 @ z[%d]=%.4f   '
      'port2 @ z[%d]=%.4f   absorber2 @ z[%d]=%.4f'
      % (idxAbs1, Zz[idxAbs1], idxPort1, Zz[idxPort1],
         idxPort2, Zz[idxPort2], idxAbs2, Zz[idxAbs2]))

abs1 = None if (NO_SHEETS or MODE == 'PML') else sheet(idxAbs1, True)

port1 = FDTD.AddWaveGuidePort(1, box_lo + [Zz.item(idxPort1)], box_hi + [Zz.item(idxPort1 + 1)],
                              'z', E_file="MSL_E.csv", H_file="MSL_H.csv",
                              kc=0.0, excite=1, excite_type=0)
port2 = FDTD.AddWaveGuidePort(2, box_lo + [Zz.item(idxPort2)], box_hi + [Zz.item(idxPort2 - 1)],
                              'z', E_file="MSL_E.csv", H_file="MSL_H.csv",
                              kc=0.0, excite=0, excite_type=0)

# refsub measures ONE sheet: a second one would put its echo in the same record
# and we would be back to untangling a standing wave.
abs2 = None if (NO_SHEETS or RULER == 'refsub') else sheet(idxAbs2, False)

# ## The reference-subtraction probe
# Between the sheet and the port, 4 cells from each: everything travelling one
# way past it is incident, everything travelling the other way is the sheet's
# echo. Present in every run, including the reference, at the same z.
idxRef = idxAbs1 + 4
REFSUB_FILE = 'refsub_ut'
_p = CSX.AddProbe(REFSUB_FILE, p_type=10, mode_file_name="MSL_E.csv")
_p.AddBox(box_lo + [Zz.item(idxRef)], box_hi + [Zz.item(idxRef)])
print('  refsub probe @ z[%d] = %.4f mm' % (idxRef, Zz[idxRef]))

# ## E-only modal fit ladder, in the source-free stretch between the ports
# The fit separates forward from backward by their phase difference across the
# ladder, so the ladder has to be an appreciable fraction of a guide wavelength
# at the LOWEST frequency of interest. The old span of 34 cells (22.7 mm) is
# 0.014 lambda_g at 0.1 GHz, where n_eff simply railed at the edge of n_scan.
# 144 cells is 96 mm, which makes the ruler usable from a few hundred MHz --
# still not 0.1 GHz, but that is a property of a 120 mm line, not of a bug.
FIT_OFFSETS = np.array([0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144])
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


# ## Divergence guard
# A run that blew up still writes probe files, and the fit below will turn 1e19
# into a confident-looking table -- that is exactly how a diverged MSL arm got
# reported as converged once already. Check the raw traces first: a healthy
# record decays, a diverged one ends at its peak.
def _stability_report():
    bad = []
    # Only probes NEAR the source. A far ladder rung carries almost no signal
    # once the far end is a PML, so its tail/peak ratio is noise and it
    # false-positived on half an aperture sweep. Non-finite values are checked
    # everywhere, though: those are unambiguous.
    for fn in fit_files:
        fp = os.path.join(Sim_Path, fn)
        if os.path.exists(fp):
            vv = np.loadtxt(fp, comments='%')[:, 1]
            if not np.all(np.isfinite(vv)):
                bad.append((fn, 'non-finite samples'))
    for fn in ['port_ut_1', REFSUB_FILE, fit_files[0]]:
        fp = os.path.join(Sim_Path, fn)
        if not os.path.exists(fp):
            continue
        v = np.loadtxt(fp, comments='%')[:, 1]
        if not np.all(np.isfinite(v)):
            bad.append((fn, 'non-finite samples'))
            continue
        # Compare the tail against the record's PEAK, not against its head: on
        # a far ladder rung the head is pre-arrival, so head-relative growth
        # just measures where the wavefront was, and it false-positives on any
        # run that converged early enough to leave a short record.
        n = len(v)
        a = np.abs(v)
        tail = a[-n // 3:].max()
        if a.max() > 0 and tail / a.max() > 0.9:
            bad.append((fn, 'not decaying: tail/peak = %.3g' % (tail / a.max())))
    return bad


_bad = _stability_report()
if _bad:
    print('\n' + '=' * 72)
    print('NOT DECAYING: the field did not die down. If these are non-finite the')
    print('run blew up and everything below is meaningless; if they are merely')
    print('flat the box is still ringing and the fit is unreliable. Probes:')
    for fn, why in _bad[:6]:
        print('   %-14s %s' % (fn, why))
    print('=' * 72)

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

# ## Reference-subtraction ruler
if REF_RUN:
    print('\nThis was the REFERENCE run. Its record is the incident wave; the')
    print('table above is meaningless for it. Now run the arm:')
    print('   MSL_ABS_MODE=%s MSL_RULER=refsub python3 %s'
          % (MODE, os.path.basename(__file__)))
elif RULER == 'refsub':
    ref_path = os.path.join(tempfile.gettempdir(), 'Test_MSL_%s_ref%d' % (MODE, PML_N))
    ref_file = os.path.join(ref_path, REFSUB_FILE)
    if not os.path.exists(ref_file):
        print('\nrefsub: no reference record at %s' % ref_file)
        print('        run this first:')
        print('   MSL_ABS_MODE=%s MSL_REF=1 MSL_PML=%d python3 %s'
              % (MODE, PML_N, os.path.basename(__file__)))
    else:
        da = np.loadtxt(os.path.join(Sim_Path, REFSUB_FILE), comments='%')
        dr = np.loadtxt(ref_file, comments='%')
        # Subtract in the frequency domain: linearity makes this identical to
        # subtracting the records, but it does not care that the two runs
        # stopped at different timesteps.
        Ua = utilities.DFT_time2freq(da[:, 0], da[:, 1], f)
        Ur = utilities.DFT_time2freq(dr[:, 0], dr[:, 1], f)
        Gref = (Ua - Ur) / Ur
        Gref_dB = 20 * np.log10(np.abs(Gref))
        _dut = ('PML_%d' % PML_DUT) if MODE == 'PML' else ('%s sheet' % MODE)
        print('\n===== reference subtraction, {}, PML_{} far end ====='
              .format(_dut, PML_N))
        print('  arm record {} samples, reference {} samples'
              .format(len(da), len(dr)))
        # The floor: difference two references built with different PML depths.
        # Whatever they disagree by is echo the PML failed to swallow, and no
        # reflection smaller than that can be resolved.
        floor_dB = None
        for alt in (32, 24, 16, 12, 8, 6):
            if alt == PML_N:
                continue
            af = os.path.join(tempfile.gettempdir(),
                              'Test_MSL_%s_ref%d' % (MODE, alt), REFSUB_FILE)
            if os.path.exists(af):
                d2 = np.loadtxt(af, comments='%')
                U2 = utilities.DFT_time2freq(d2[:, 0], d2[:, 1], f)
                floor_dB = 20 * np.log10(np.abs((U2 - Ur) / Ur))
                print('  floor from the PML_%d vs PML_%d references' % (PML_N, alt))
                break

        print('  f[GHz]   |Gamma| dB     floor dB    verdict')
        for fi in [0.10, 0.25, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00]:
            j = np.abs(f - fi * 1e9).argmin()
            if floor_dB is None:
                print('  {:6.2f}     {:8.2f}         n/a      floor not measured'
                      .format(f[j] / 1e9, Gref_dB[j]))
            else:
                margin = Gref_dB[j] - floor_dB[j]
                verdict = ('FLOOR-LIMITED' if margin < 6.0
                           else 'ok (%.0f dB above floor)' % margin)
                print('  {:6.2f}     {:8.2f}     {:8.2f}    {}'
                      .format(f[j] / 1e9, Gref_dB[j], floor_dB[j], verdict))
        print('  worst in band: {:.2f} dB at {:.2f} GHz'
              .format(np.max(Gref_dB), f[np.argmax(Gref_dB)] / 1e9))
        if floor_dB is None:
            print('  NO FLOOR MEASURED -- these numbers are upper bounds of unknown')
            print('  tightness. Build a second reference: MSL_REF=1 MSL_PML=16')
        print('  The MODAL-FIT column is NOT comparable here: in this single-ended')
        print('  configuration the ladder sits downstream of the source with the')
        print('  PML beyond it, so its backward wave comes from the PML, not from')
        print('  the sheet under test.')

if os.environ.get('MPLBACKEND', '') != 'Agg':
    figure()
    plot(f / 1e9, Gam_dB, 'm-', linewidth=2, label=r'$|\Gamma|$ modal fit (E-only)')
    plot(f / 1e9, s11_dB, 'k--', linewidth=1.5, label='$S_{11}$ from port')
    plot(f / 1e9, s21_dB, 'b-', linewidth=2, label='$S_{21}$ (through)')
    grid(); legend(); ylim(-70, 5)
    ylabel('S-Parameter (dB)'); xlabel('Frequency (GHz)')
    title('MSL quasi-TEM — {} termination'.format(MODE))
    show()
