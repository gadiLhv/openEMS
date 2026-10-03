"""Short coax between two IPML waveguide ports: the original vs. a refined port-window mesh.

The colleague's coax_ipml_repro.py, configured from the block below instead of
command-line flags, with one switch for the transverse mesh:

  MESH = 'original'   the array's window lines (0.75-1.06 mm cells across the
                      r = 4 mm opening): an out-of-band mode (~40 GHz; the
                      decimated port probe shows it aliased at ~8.1 GHz) grows
                      without bound
  MESH = 'refined'    uniform REFINED_RES cells across the window: stable

Each mesh runs into its own folder, coax_mesh_compare/<MESH>/, kept after the
run, with time-domain E dumps (VTK) on two planes for a side-by-side look in
ParaView: through the coax axis (xz), and across the line at mid-length (xy).

  coax       pin r = 1.12 mm, opening r = 4.0 mm through a PEC slab, eps_r 2.2
  ports      8 x 8 mm windows (circumscribe the opening), propagation along z
  absorbers  PML_8 sheet the size of the window on each port plane, and a PEC
             block of the window one mesh cell deep behind it
  source     Gaussian 2-3 GHz (f0 2.5, fc 0.5 GHz), zero-mean
"""
import glob, os, shutil
import numpy as np
from CSXCAD import ContinuousStructure, AppCSXCAD_BIN
from CSXCAD.CSProperties import ABCtype
from openEMS import openEMS

# ============================== CONFIG =======================================
MESH = 'refined'  # 'original' (the array's window lines) or 'refined'
REFINED_RES = 0.5  # mm, uniform window cells for MESH = 'refined' (0.5 and 0.25 were stable)

LENGTH = 6.0  # mm, line length
EPS = 2.2  # coax dielectric eps_r
MAX_NS = 60.0  # run length; the energy end-criterion is off
TSF = None  # timestep factor, None = openEMS default
SHEETS = 'both'  # 'both', '1', '2' or 'none'
FULL_SHEET = False  # sheet across the whole domain cross-section
BACKING = True  # PEC block behind each port
ZERO_MEAN = True  # zero time-integral excitation

DUMP_FIELDS = False  # time-domain E on the xz (axis) and xy (mid-line) planes
SHOW_CSXCAD = True  # write the XML and open AppCSXCAD before the run
RUN_FDTD = True
# =============================================================================

A, B = 1.12, 4.0  # pin radius, opening radius (mm)
W = 4.0  # half-width of the port window (mm)
L = LENGTH
F0, FC = 2.5e9, 0.5e9  # Blit: centre and half-span of 2-3 GHz
EXC_PEAK_S = 9.0 / (2 * np.pi * FC)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'coax_mesh_compare', MESH)
assert MESH in ('original', 'refined')


def graded(start, step, ratio, maxstep, length):
    """Lines from `start` outward by `length`, growing from `step` at `ratio`."""
    out, x, h = [], start, step
    while abs(x - start) < length:
        h = min(h * ratio, maxstep)
        x += h
        out.append(x)
    return out


# ---- mesh (mm): the port window, graded out to padding + PML_8
if MESH == 'original':
    win = [-4.0, -3.2487, -2.1844, -1.12, -0.3733, 0.3733, 1.12, 2.1844, 3.2487, 4.0]
else:
    win = list(np.linspace(-W, W, int(round(2 * W / REFINED_RES)) + 1))
side = graded(W, win[-1] - win[-2], 1.3, 3.0, 16.0)
XY = sorted(set(np.round([-x for x in side] + win + side, 4)))
h = 0.4083  # the array's cells in front of a port
n_line = max(3, int(round(L / h)))
line_z = list(np.linspace(0.0, L, n_line + 1))  # uniform, so >= 3 equal cells in front
back = 1.2583  # the cell behind each port, as in the array
pad = graded(back, back, 1.3, 3.0, 12.0)
Z = sorted(set(np.round([-back] + [-z for z in pad] + line_z + [L + back] + [L + z for z in pad], 4)))

FDTD = openEMS(EndCriteria=1e-30, MaxTime=MAX_NS * 1e-9)
if TSF:
    FDTD.SetTimeStepFactor(TSF)
FDTD.SetGaussExcite(F0, FC)
if ZERO_MEAN:
    FDTD.SetExciteZeroMean(True)
FDTD.SetBoundaryCond(["PML_8"] * 6)
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
g = CSX.GetGrid()
g.SetDeltaUnit(1e-3)
g.SetLines("x", XY); g.SetLines("y", XY); g.SetLines("z", Z)

# ---- geometry: dielectric everywhere along the line, a PEC slab with a round
# opening (the opening is the dielectric cylinder at higher priority), the pin.
xmax = max(XY)
diel = CSX.AddMaterial("diel", epsilon=EPS)
diel.AddBox([-xmax, -xmax, 0], [xmax, xmax, L], priority=1)
pec = CSX.AddMetal("PEC")
pin = CSX.AddMetal("PIN")

pec.AddBox([-xmax, -xmax, 0], [xmax, xmax, L], priority=10)  # the slab
diel.AddCylinder([0, 0, 0], [0, 0, L], B, priority=20)  # its opening
pin.AddCylinder([0, 0, -0.0], [0, 0, L], A, priority=30)  # the pin

# ---- analytic TEM mode, written on a regular grid over the window
#      (u, v relative to the window corner; u varies fastest)
os.makedirs(OUT, exist_ok=True)
uu = np.linspace(0, 2 * W, 161)
U, V = np.meshgrid(uu, uu, indexing="xy")
X, Y = U - W, V - W
R = np.hypot(X, Y)
inside = (R > A) & (R < B)
Er = np.where(inside, 1.0 / (np.where(R > 0, R, 1) * np.log(B / A)), 0.0)
Ex, Ey = Er * X / np.where(R > 0, R, 1), Er * Y / np.where(R > 0, R, 1)
eta = 376.73 / np.sqrt(EPS)
Hx, Hy = -Ey / eta, Ex / eta  # H = z x E / eta
e_file, h_file = os.path.join(OUT, "mode_E.csv"), os.path.join(OUT, "mode_H.csv")
np.savetxt(e_file, np.column_stack([U.ravel(), V.ravel(), Ex.ravel(), Ey.ravel()]), delimiter=",")
np.savetxt(h_file, np.column_stack([U.ravel(), V.ravel(), Hx.ravel(), Hy.ravel()]), delimiter=",")

zi = lambda v: int(np.argmin(np.abs(np.asarray(Z) - v)))
ports = []
for nr, zp, d in ((1, 0.0, +1), (2, L, -1)):
    i0 = zi(zp)
    if SHEETS in ("both", str(nr)):
        s = CSX.AddAbsorbingBC("port_%d_ipml" % nr, NormalSignPositive=(d == 1),
                               AbsorbingBoundaryType=ABCtype.PML_8)
        sw = xmax if FULL_SHEET else W
        s.AddBox([-sw, -sw, zp], [sw, sw, zp], priority=60)
    if BACKING:
        pec.AddBox([-W, -W, Z[i0 - d]], [W, W, zp], priority=50)  # one cell, outward
    ports.append(FDTD.AddWaveGuidePort(
        nr, [-W, -W, Z[i0 + d]], [W, W, Z[i0 + 2 * d]], "z",
        kc=0.0, excite=1.0 if nr == 1 else 0.0, excite_type=0,
        E_file=e_file, H_file=h_file, mode_type="TEM"))

if DUMP_FIELDS:
    # by-eye comparison: the coax axis (xz) and a cross-section at mid-line (xy),
    # both reaching one window-width past the opening
    xz = CSX.AddDump('Et_xz_axis', dump_type=0, file_type=0)
    xz.AddBox([-2 * W, 0, Z[0]], [2 * W, 0, Z[-1]])
    xy = CSX.AddDump('Et_xy_midline', dump_type=0, file_type=0)
    xy.AddBox([-2 * W, -2 * W, L / 2], [2 * W, 2 * W, L / 2])

print("MESH = %s: %d x %d x %d; window lines %s; line %.3f mm in %d cells"
      % (MESH, len(XY), len(XY), len(Z), [round(float(v), 4) for v in win], L, len(line_z) - 1))

sim = os.path.join(OUT, "sim")
if os.path.exists(sim):
    shutil.rmtree(sim)
os.makedirs(sim)
CSX_file = os.path.join(sim, 'coax_ipml.xml')
CSX.Write2XML(CSX_file)
if SHOW_CSXCAD:
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))
if not RUN_FDTD:
    raise SystemExit
FDTD.Run(sim, cleanup=False)

# ---- verdict (the port probe is decimated: its "dominant frequency" can be an alias)
t_end, rep = 0.0, []
for f in sorted(glob.glob(os.path.join(sim, "port_ut*"))):
    dd = np.loadtxt(f, comments="%")
    t, u = dd[:, 0], dd[:, 1]
    t_end = max(t_end, t[-1])
    drive = np.nanmax(np.abs(u[t <= EXC_PEAK_S])) if (t <= EXC_PEAK_S).any() else np.nan
    q = t >= 0.75 * t[-1]
    late = np.nanmax(np.abs(u[q]))
    seg = u[q] - u[q].mean()
    sp = np.abs(np.fft.rfft(seg * np.hanning(len(seg))))
    fq = np.fft.rfftfreq(len(seg), t[1] - t[0])
    rep.append((os.path.basename(f), drive, late, fq[int(np.argmax(sp))] / 1e9,
                (not np.isfinite(u).all()) or late > drive))
cut = t_end < 0.98 * MAX_NS * 1e-9
grew = cut or any(r[4] for r in rep)
print("\nRESULT mesh=%s eps=%.2f length=%.3f sheets=%s backing=%s zero_mean=%s tsf=%s  ran %.2f of %.2f ns  -> %s"
      % (MESH, EPS, L, SHEETS, BACKING, ZERO_MEAN, TSF, t_end * 1e9, MAX_NS,
         "GROWS / DIVERGED" if grew else "stable"))
for name, drive, late, fdom, g_ in rep:
    print("    %-10s max|u| to the peak %.3g, last quarter %.3g, last-quarter dominant %.2f GHz%s"
          % (name, drive, late, fdom, "  <- grew" if g_ else ""))
print("fields: %s/Et_*.vtr" % sim)
