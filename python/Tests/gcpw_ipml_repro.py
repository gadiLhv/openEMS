"""MSL -> GCPW line terminated by invisible-PML sheets: does it diverge?

Plain openEMS + CSXCAD + numpy. Boxes only. The geometry, mesh lines, boundary
conditions and excitation are exactly what Blit hands openEMS for project
"MSL_Turned_GCPW (Gadi copy)" -- the mesh lines below are Blit's, verbatim.
Units: mm.

  python3 gcpw_ipml_repro.py                   # lumped source, IPML sheets on both ports
  python3 gcpw_ipml_repro.py --no-ipml         # control: same, sheets removed
  python3 gcpw_ipml_repro.py --drive wg        # Blit's waveguide ports (port_*_E/H.csv)
  python3 gcpw_ipml_repro.py --tsf 0.7         # scale openEMS's timestep
  python3 gcpw_ipml_repro.py --sheets 1        # sheet on port 1 (z=0) only; 2, both, none
  python3 gcpw_ipml_repro.py --refine-ports    # 1 mm cells in front of each sheet (Blit: 3.15)
  python3 gcpw_ipml_repro.py --max-ns 3        # stop at 3 ns (a stable run otherwise goes ~200 ns)

It reports where the run stopped against the excitation's peak, and the largest
port voltage: a diverged run stops before the peak with |u| in the 1e6..1e12
range; a healthy one runs well past it.
"""
import argparse, glob, os, shutil, tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from CSXCAD.CSProperties import ABCtype
from openEMS import openEMS

MESH_X = [-103.885, -99.4432, -95.0018, -90.5605, -86.1191, -81.6777, -77.2364, -72.795, -68.3536, -62.4215, -56.4893, -50.5572, -44.625, -38.6929, -33.3601, -29.2159, -26.3077, -24.2983, -22.9526, -21.9684, -20.9842, -20, -19.0135, -18.027, -17.0406, -15.6845, -13.6546, -10.8001, -7.97274, -5.91565, -4.55953, -3.57305, -2.58658, -1.6001, -0.80005, 0, 0.80005, 1.6001, 2.58658, 3.57305, 4.55953, 5.91565, 7.97274, 10.8001, 13.6546, 15.6845, 17.0406, 18.027, 19.0135, 20, 20.9842, 21.9684, 22.9526, 24.2983, 26.3077, 29.2159, 33.3601, 38.6929, 44.625, 50.5572, 56.4893, 62.4215, 68.3536, 72.795, 77.2364, 81.6777, 86.1191, 90.5605, 95.0018, 99.4432, 103.885]
MESH_Y = [-83.9846, -79.5432, -75.1018, -70.6605, -66.2191, -61.7777, -57.3364, -52.895, -48.4536, -42.5215, -36.5893, -30.6572, -24.725, -18.7929, -13.4601, -9.31588, -6.40772, -4.39835, -3.05258, -2.06839, -1.08419, -0.1, 0.45, 1, 1.98419, 2.96839, 3.95258, 5.29835, 7.30772, 9, 10.2159, 14.3601, 19.6929, 25.625, 31.5572, 37.4893, 43.4215, 49.3536, 53.795, 58.2364, 62.6777, 67.1191, 71.5605, 76.0018, 80.4432, 84.8846]
MESH_Z = [-83.8846, -79.4432, -75.0018, -70.5605, -66.1191, -61.6777, -57.2364, -52.795, -48.3536, -41.8958, -35.438, -28.9802, -22.5224, -16.3106, -10.4764, -4.65362, 0, 3.15418, 6.30836, 9.46254, 12.6167, 15.7709, 18.9251, 22.0793, 25.2334, 28.3876, 31.5418, 34.696, 37.8502, 41.0043, 44.1585, 47.3127, 50.4669, 53.4733, 55.6035, 57.0087, 58.0058, 59.0029, 60, 60.9971, 61.9942, 62.9913, 64.3965, 66.5267, 69.5331, 72.6873, 75.8415, 78.9957, 82.1498, 85.304, 88.4582, 91.6124, 94.7666, 97.9207, 101.075, 104.229, 107.383, 110.537, 113.692, 116.846, 120, 124.654, 130.476, 136.311, 142.522, 148.98, 155.438, 161.896, 168.354, 172.795, 177.236, 181.678, 186.119, 190.56, 195.002, 199.443, 203.885]

F_MIN, F_MAX = 0.1e9, 3.0e9
F0, FC = (F_MIN + F_MAX) / 2, (F_MAX - F_MIN) / 2  # Blit: centre and half-span
EXC_PEAK_S = 9.0 / (2 * np.pi * FC)  # openEMS Gaussian peak

ap = argparse.ArgumentParser()
ap.add_argument("--drive", choices=("lumped", "wg"), default="lumped")
ap.add_argument("--no-ipml", action="store_true")
ap.add_argument("--tsf", type=float, default=None)
ap.add_argument("--sheets", choices=("both", "1", "2", "none"), default="both")
ap.add_argument("--refine-ports", action="store_true")
ap.add_argument("--max-ns", type=float, default=None)
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))

max_time = 100 * 9.0 / (np.pi * FC)  # Blit: 100 pulse widths
if a.max_ns:
    max_time = a.max_ns * 1e-9
if a.no_ipml:
    a.sheets = "none"
if a.refine_ports:
    # Replace the 3.15 mm cells next to each port plane with 1 mm ones, four
    # cells deep, so the sheet's virtual cells inherit 1 mm instead.
    MESH_Z = sorted(set([z for z in MESH_Z if not (0 < z < 4) and not (116 < z < 120)]
                        +[1.0, 2.0, 3.0, 4.0, 116.0, 117.0, 118.0, 119.0]))
# With --max-ns the energy end-criterion is switched OFF, so a healthy run
# always reaches the cap. openEMS leaves its loop when `change > endCrit` fails,
# change = currE/maxE, and diverged fields make that NaN -- so stopping before
# the cap is then a deterministic sign of divergence, wherever it starts.
FDTD = openEMS(EndCriteria=1e-30 if a.max_ns else 1e-5, MaxTime=max_time)
if a.tsf:
    FDTD.SetTimeStepFactor(a.tsf)
FDTD.SetGaussExcite(F0, FC)
FDTD.SetExciteZeroMean(True)  # zero time-integral: the pulse leaves no static charge at the source
FDTD.SetBoundaryCond(["PML_8"] * 6)
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
g = CSX.GetGrid()
g.SetDeltaUnit(1e-3)
g.SetLines("x", MESH_X); g.SetLines("y", MESH_Y); g.SetLines("z", MESH_Z)

# ---- geometry: every Blit STL of this model is an exact box
fr4 = CSX.AddMaterial("FR4", epsilon=4.37)
fr4.AddBox([-20, 0.0, 0], [20, 1.0, 120], priority=1)  # Substrate
pec = CSX.AddMetal("PEC")
pec.AddBox([-20, -0.1, 0], [20, 0.0, 120], priority=10)  # Gnd (bottom)
pec.AddBox([-1, 1.0, 0], [1, 1.1, 120], priority=10)  # Trace
pec.AddBox([-20, 1.0, 0], [-1.6, 1.1, 60], priority=10)  # Top_GND, left
pec.AddBox([1.6, 1.0, 0], [20, 1.1, 60], priority=10)  # Top_GND, right

pec = CSX.AddMetal("VIA")
# Ground connecting blocks
pec.AddBox([-20, 0.0, 4.5], [-4.5, 1.0, 55.5], priority=11)  # Top_GND, left
pec.AddBox([4.5, 0.0, 4.5], [20, 1.0, 55.5], priority=11)  # Top_GND, right


def line_index(lines, v):
    return int(np.argmin(np.abs(np.asarray(lines) - v)))


# The port window: x -8..8, y -0.1..9 snapped to the mesh. The IPML sheets and
# the waveguide ports share it.
xw = (MESH_X[line_index(MESH_X, -8)], MESH_X[line_index(MESH_X, 8)])
yw = (MESH_Y[line_index(MESH_Y, -0.1)], MESH_Y[line_index(MESH_Y, 9)])
N_BOUNDARY_PML = 8  # PML_8 on all six faces, see SetBoundaryCond above


def clear_of_boundary_pml(start, stop):
    # At least one mesh cell between the box and the boundary PML, on every face.
    for lines, lo, hi, axis in zip((MESH_X, MESH_Y, MESH_Z), start, stop, "xyz"):
        i_lo, i_hi = line_index(lines, min(lo, hi)), line_index(lines, max(lo, hi))
        if i_lo < N_BOUNDARY_PML + 1 or i_hi > len(lines) - 1 - (N_BOUNDARY_PML + 1):
            raise ValueError("box %s..%s touches the boundary PML along %s" % (start, stop, axis))


def ipml_sheet(name, z, domain_positive):
    # The sheet covers the port window on the port plane. Behind it (away from
    # the domain) is a PEC block of the same window, one mesh cell deep.
    s = CSX.AddAbsorbingBC(name, NormalSignPositive=domain_positive,
                           AbsorbingBoundaryType=ABCtype.PML_8)
    s.AddBox([xw[0], yw[0], z], [xw[1], yw[1], z], priority=60)
    i0 = line_index(MESH_Z, z)
    zb = MESH_Z[i0 - 1] if domain_positive else MESH_Z[i0 + 1]
    clear_of_boundary_pml([xw[0], yw[0], zb], [xw[1], yw[1], z])
    pec.AddBox([xw[0], yw[0], zb], [xw[1], yw[1], z], priority=50)


# Port 1 at z = 0 (GCPW end, guide on +z), port 2 at z = 120 (MSL end, guide on -z).
ports = []
if a.sheets in ("both", "1"):
    ipml_sheet("port_1_ipml", 0.0, True)
if a.sheets in ("both", "2"):
    ipml_sheet("port_2_ipml", 120.0, False)

if a.drive == "lumped":
    # A 50-ohm lumped source trace -> bottom ground, in the GCPW section.
    zs = MESH_Z[line_index(MESH_Z, 30.0)]
    ports.append(FDTD.AddLumpedPort(1, 50, [-1, 0.0, zs], [1, 1.0, zs], "y", excite=1.0))
else:
    # Blit's waveguide ports on the port window: excitation one cell in front
    # of the sheet, measurement two cells in.
    for nr, zp, d in ((1, 0.0, +1), (2, 120.0, -1)):
        i0 = line_index(MESH_Z, zp)
        ports.append(FDTD.AddWaveGuidePort(
            nr, [xw[0], yw[0], MESH_Z[i0 + d]], [xw[1], yw[1], MESH_Z[i0 + 2 * d]], "z",
            kc=0.0, excite=1.0 if nr == 1 else 0.0, excite_type=0,
            E_file=os.path.join(HERE, "port_%d_E.csv" % nr),
            H_file=os.path.join(HERE, "port_%d_H.csv" % nr), mode_type="TEM"))

sim = tempfile.mkdtemp(prefix="gcpw_ipml_")

# Define dump box...
SimBox = [MESH_X[0], MESH_X[-1], MESH_Y[0], MESH_Y[-1], MESH_Z[0], MESH_Z[-1]]
Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
start = [float(SimBox[0]), float(SimBox[2]), float(SimBox[4])];
stop = [float(SimBox[1]), float(SimBox[3]), float(SimBox[5])];
Et.AddBox(start, stop);

if 1:  # debugging only
    CSX_file = os.path.join(sim, 'gcpw_ipml_repro.xml')
    if not os.path.exists(sim):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

FDTD.Run(sim, cleanup=True)

# ---- verdict, from the port voltage probes.
# Diverged if the voltage goes non-finite, or grows well past what the source
# drove (latest-quarter max vs max up to the excitation peak). Stopping before
# the cap alone is not enough to judge either way: a NaN in the energy makes
# openEMS leave its loop exactly as if it had converged.
t_end, verdicts = 0.0, []
for f in sorted(glob.glob(os.path.join(sim, "port_ut*"))):
    d = np.loadtxt(f, comments="%")
    t, u = d[:, 0], np.abs(d[:, 1])
    t_end = max(t_end, t[-1])
    drive = np.nanmax(u[t <= EXC_PEAK_S]) if (t <= EXC_PEAK_S).any() else np.nan
    late = np.nanmax(u[t >= 0.75 * t[-1]])
    grew = (not np.isfinite(u).all()) or (np.isfinite(drive) and drive > 0 and late > 100 * drive)
    verdicts.append((os.path.basename(f), drive, late, grew))
cut_short = (t_end < 0.98 * max_time) if a.max_ns else (t_end < EXC_PEAK_S)
diverged = cut_short or any(v[3] for v in verdicts)
print("\nRESULT drive=%s sheets=%s refine=%s tsf=%s  stopped at %.3f ns of %.3f  "
      "(excitation peak %.3f ns)  -> %s" % (a.drive, a.sheets, a.refine_ports, a.tsf,
      t_end * 1e9, max_time * 1e9, EXC_PEAK_S * 1e9, "DIVERGED" if diverged else "stable"))
for name, drive, late, grew in verdicts:
    print("    %-10s max|u| up to the peak %.3g V, in the last quarter %.3g V%s"
          % (name, drive, late, "  <- grew" if grew else ""))
shutil.rmtree(sim, ignore_errors=True)
