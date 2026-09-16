"""Learn the CPW mode from openEMS itself, instead of from the mode solver.

A long CPW, centre-fed and PML-terminated at both ends, is driven with a CW
source. The steady-state field is sampled a few cells away from the source,
directly on the simulation mesh, as time-domain plane dumps turned into a phasor
by a DFT over the last whole periods (see the dump block). That field IS the
discrete mode this mesh supports, so the resulting E/H files carry no FEM
continuum, corner or interpolation error at all.

Cross-section and x/z mesh are copied verbatim from CPW_W_ModalAbsorb.py,
including f0/fc_exc, which set mesh_res and therefore the smoothed x/z lines.
Change one and the learned file no longer lines up with the test.

Usage:  python3 CPW_Learn_Mode.py          run, then build the files
        CPW_LEARN_POST_ONLY=1 python3 ...   rebuild the files from the last run
"""
import os, tempfile
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

POST_ONLY = bool(int(os.environ.get('CPW_LEARN_POST_ONLY', '0')))
Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_CPW_LearnMode')
if not os.path.exists(Sim_Path):
    os.mkdir(Sim_Path)

# ## Geometry -- identical to CPW_W_ModalAbsorb.py
Line_W = 1.15
CPW_gap = 0.3
cu_thick = 0.1
substrate_epsR = 4.3
substrate_width = 11.0
substrate_thickness = 1.0
substrate_cells = 4
gap_cells = 3
trace_cells = 7
port_w_fact = 7.5
port_h_fact = 3.5
Airbox_Add = 12.5
unit = 1e-3
f0, fc_exc = 1.55e9, 1.45e9          # ONLY for mesh_res -- keeps the x/z mesh identical

# ## Learning run
F_CW = 2e9                           # CW frequency; same f0 the mode solver used
N_HALF = 75                          # y cells each side of the source (~75 mm)
N_PERIODS = 40                       # simulated CW periods
SAMPLE_OFFSETS = (5, 10, -5)         # planes, in cells from the source; +5 is the one used

mesh_res = ((C0 / (f0 + fc_exc)) / unit) / 50
dy = mesh_res / 2.0
L_half = N_HALF * dy

CSX = ContinuousStructure()
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

SimBox = np.array([
    -substrate_width * 0.5 - Airbox_Add, substrate_width * 0.5 + Airbox_Add,
    -L_half, L_half,
    -substrate_thickness * (port_h_fact - 1.0) - Airbox_Add,
     substrate_thickness * (1.0 + port_h_fact) + Airbox_Add])
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

y0, y1 = SimBox[2], SimBox[3]        # every layer runs through both PMLs: no open ends

line = CSX.AddMetal('cu_top')
line.AddBox(priority=20, start=[-Line_W / 2, y0, substrate_thickness],
            stop=[Line_W / 2, y1, substrate_thickness + cu_thick])
mesh.AddLine('x', [-Line_W / 2, Line_W / 2])
mesh.AddLine('z', [substrate_thickness, substrate_thickness + cu_thick])
mesh.AddLine('x', np.linspace(-0.5 * Line_W, 0.5 * Line_W, trace_cells))

p_x = Line_W * (1.0 + port_w_fact) * 0.5
p_z0 = -substrate_thickness * (port_h_fact - 1.0)
p_z1 = substrate_thickness * (1.0 + port_h_fact)
mesh.AddLine('x', [-p_x, p_x])
mesh.AddLine('z', [p_z0, p_z1])

sub = CSX.AddMaterial('FR4', epsilon=substrate_epsR)
sub.AddBox(priority=2, start=[-substrate_width / 2, y0, 0.0],
           stop=[substrate_width / 2, y1, substrate_thickness])
mesh.AddLine('x', [-substrate_width / 2, substrate_width / 2])
mesh.AddLine('z', [0.0, substrate_thickness])
mesh.AddLine('z', np.linspace(0, substrate_thickness, substrate_cells + 1))

gnd = CSX.AddMetal('cu_gnd')
gnd.AddBox(priority=10, start=[-0.5 * substrate_width, y0, substrate_thickness],
           stop=[-0.5 * (Line_W + CPW_gap * 2), y1, substrate_thickness + cu_thick])
gnd.AddBox(priority=10, start=[0.5 * (Line_W + CPW_gap * 2), y0, substrate_thickness],
           stop=[0.5 * substrate_width, y1, substrate_thickness + cu_thick])
mesh.AddLine('x', [-0.5 * substrate_width, -0.5 * (Line_W + CPW_gap * 2),
                   0.5 * (Line_W + CPW_gap * 2), 0.5 * substrate_width])
mesh.AddLine('x', np.linspace(-0.5 * (Line_W + CPW_gap * 2), -Line_W * 0.5, gap_cells))
mesh.AddLine('x', np.linspace(Line_W * 0.5, 0.5 * (Line_W + CPW_gap * 2), gap_cells))

mesh.AddLine('y', (np.arange(-N_HALF, N_HALF + 1) * dy).tolist())
mesh.SmoothMeshLines('x', mesh_res, 1.4)
mesh.SmoothMeshLines('z', mesh_res, 1.4)

Yz = np.unique(np.asarray(mesh.GetLines('y')))
iy_src = int(np.argmin(np.abs(Yz)))
assert abs(Yz[iy_src]) < 1e-9, 'no mesh line at the source plane'
assert np.allclose(np.diff(Yz), dy), 'y mesh is not uniform'

# ## Source: even-mode gap voltage, straight down the middle.
#  E_x in the two slots with opposite signs -- the CPW (even) mode, not the
#  slotline (odd) one. Signs follow CPW_E.csv: Ex > 0 in the x < 0 slot.
GI = 0.5 * (Line_W + CPW_gap * 2)
LW2 = 0.5 * Line_W
for name, sgn, xa, xb in (('exc_slot_lo', +1.0, -GI, -LW2), ('exc_slot_hi', -1.0, LW2, GI)):
    e = CSX.AddExcitation(name, exc_type=0, exc_val=[sgn, 0, 0])
    e.AddBox([xa, 0.0, substrate_thickness], [xb, 0.0, substrate_thickness + cu_thick])

# ## TIME-domain plane dumps, window only, raw Yee values and node-interpolated.
#  NOT frequency-domain dumps: ProcessFieldsFD samples every Nyquist timestep
#  (processing.cpp:164), and for a CW source openEMS trims dt so one period is
#  exactly two Nyquist intervals. The "phasor" is then two samples per period at
#  fixed phases -- a real snapshot of sin(wt - beta*y), with no phase progression
#  along y and a Poynting ratio of exactly 1. Measured: n_eff = 0.006 and a 3x
#  power gain over 5 mm. The DFT is done here instead, over whole periods.
OVERSAMPLING = 8                     # TD dump every Nyquist/8 steps -> 16 per period
LAST_PERIODS = 10                    # steady-state window for the DFT
for off in SAMPLE_OFFSETS:
    yp = Yz[iy_src + off]
    for fld, dtyp in (('E', 0), ('H', 1)):
        for tag, dmode in (('raw', 0), ('node', 1)):
            d = CSX.AddDump('%s_%s_%+d' % (fld, tag, off), dump_type=dtyp, dump_mode=dmode,
                            file_type=1)
            d.AddBox([-p_x, yp, p_z0], [p_x, yp, p_z1])

if not POST_ONLY:
    FDTD = openEMS(NrTS=int(1e8), EndCriteria=0, MaxTime=N_PERIODS / F_CW,
                   OverSampling=OVERSAMPLING)
    FDTD.SetSinusExcite(F_CW)
    FDTD.SetBoundaryCond(['MUR', 'MUR', 'PML_8', 'PML_8', 'MUR', 'MUR'])
    FDTD.SetCSX(CSX)
    print('learning run: %d y lines, source at y[%d], dy = %.5f mm, %d CW periods at %.2f GHz'
          % (len(Yz), iy_src, dy, N_PERIODS, F_CW / 1e9))
    FDTD.Run(Sim_Path, verbose=3, cleanup=False)


# ############################################################################
# ## Build the mode files from the steady state
# ############################################################################
import h5py

W = 2.0 * np.pi * F_CW


def phasor(name):
    """DFT at F_CW over the last LAST_PERIODS whole periods of a TD plane dump."""
    f = h5py.File(os.path.join(Sim_Path, name + '.h5'), 'r')
    td = f['/FieldData/TD']
    keys = sorted(td.keys(), key=int)
    t = np.array([td[k].attrs['time'] for k in keys]).ravel()
    sel = t >= t[-1] - LAST_PERIODS / F_CW + 0.5 / F_CW / (2 * OVERSAMPLING)
    acc = None
    for k, tk in zip(np.array(keys)[sel], t[sel]):
        v = td[k][:].astype(float)[:, :, 0, :] * np.exp(-1j * W * tk)
        acc = v if acc is None else acc + v
    P = 2.0 * acc / np.count_nonzero(sel)
    x = f['/Mesh/x'][:] / unit
    z = f['/Mesh/z'][:] / unit
    y = f['/Mesh/y'][0] / unit
    return P, x, y, z, np.count_nonzero(sel)


def cosim(a, b):
    return abs(np.vdot(a.ravel(), b.ravel())) / (np.linalg.norm(a) * np.linalg.norm(b))


E5, x, y5, z, nsamp = phasor('E_node_+5')
H5 = phasor('H_node_+5')[0]
E10, _, y10, _, _ = phasor('E_node_+10')
Em5 = phasor('E_node_-5')[0]
print('\nDFT over %d samples (%d periods x %d/period)' % (nsamp, LAST_PERIODS, nsamp // LAST_PERIODS))

# ---- is it a clean forward travelling mode? -------------------------------
dA = np.outer(np.abs(np.gradient(z)), np.abs(np.gradient(x))) * unit * unit
Sy = 0.5 * np.sum((E5[2] * np.conj(H5[0]) - E5[0] * np.conj(H5[2])) * dA)
dphi = np.angle(np.vdot(E5.ravel(), E10.ravel()))
beta = -dphi / ((y10 - y5) * unit)
print('beta +5 -> +10  : %.4f rad/m  -> n_eff = %.5f   (solver 1.50347, modal fit 1.5400)'
      % (beta, beta * C0 / W))
print('|E(+10)|/|E(+5)|: %.5f   (1 for a lossless forward wave)' % (np.linalg.norm(E10) / np.linalg.norm(E5)))
print('shape +5 vs +10 : %.6f   +5 vs -5: %.6f' % (cosim(E5, E10), cosim(E5, Em5)))
print('Re(S)/|S| at +5 : %+.5f' % (Sy.real / abs(Sy)))

# ---- real mode shape: strip the common phase ------------------------------
def realise(P):
    phi = 0.5 * np.angle(np.sum(P * P))
    R = P * np.exp(-1j * phi)
    return R.real, np.linalg.norm(R.imag) / np.linalg.norm(R)

eE, impE = realise(E5)
eH, impH = realise(H5)
print('residual imaginary part after phase removal: E %.2e, H %.2e' % (impE, impH))

# H sign: the pair must carry power in +y, i.e. Ez*Hx - Ex*Hz > 0
if np.sum((eE[2] * eH[0] - eE[0] * eH[2]) * dA) < 0:
    eH = -eH
# E sign: match the existing file's convention (Ex > 0 in the x < 0 slot)
ref = np.loadtxt('CPW_E.csv', delimiter=',')
zl, xl = np.unique(ref[:, 0]), np.unique(ref[:, 1])
assert len(zl) == len(z) and len(xl) == len(x) and \
    np.allclose(zl, z - z[0], atol=1e-6) and np.allclose(xl, x - x[0], atol=1e-6), \
    'learned plane does not sit on the CPW_E.csv grid'
refEz = ref[:, 2].reshape(len(xl), len(zl)).T
refEx = ref[:, 3].reshape(len(xl), len(zl)).T
if np.sum(eE[0] * refEx + eE[2] * refEz) < 0:
    eE, eH = -eE, -eH
print('shape vs CPW_E.csv (edge-integrated solver file): %.5f'
      % cosim(np.stack([eE[2], eE[0]]), np.stack([refEz, refEx])))

# scale both by one factor so E matches the existing file's peak (cosmetic only)
s_fac = np.abs(np.stack([refEz, refEx])).max() / np.abs(eE).max()
eE, eH = eE * s_fac, eH * s_fac


def write(fn, F):
    # y-normal port: col0 = z-local, col1 = x-local (col0 fastest), col2 = F_z, col3 = F_x
    Zg, Xg = np.meshgrid(z - z[0], x - x[0], indexing='ij')
    out = np.column_stack([Zg.ravel(order='F'), Xg.ravel(order='F'),
                           F[2].ravel(order='F'), F[0].ravel(order='F')])
    np.savetxt(fn, out, delimiter=',')
    print('wrote %s (%d rows)' % (fn, len(out)))


write('CPW_E_sim.csv', eE)
write('CPW_H_sim.csv', eH)
