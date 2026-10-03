"""In-band reflection of the boundary PML terminating a coax: |Gamma(f)| by reference subtraction.

The coax of minimal_coax_pml.py (pin r 1.12, outer r 4.0, Dk 2.2, PEC box)
runs into the domain's PML_8 at z-min. A soft source launches a pulse toward
it; an E probe sits PROBE_DIST in front of the PML. The reference run moves
the z-min end REF_EXTRA further away (PEC), so nothing comes back from that end
within the run. The difference of the two probe signals is the PML's
reflection: |Gamma| = |FFT(run - ref)| / |FFT(ref)|. The run is a few ns long,
too short for any late-time growth to matter. The PML settings come from the
environment (UPML_CFS, UPML_CFS_KMAX, ...), as for minimal_coax_pml.py.
"""
import os, shutil, tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

# ============================== CONFIG =======================================
XY_RES = 1.0           # mm, transverse cells
DZ = 0.245             # mm, cells along the coax (and in the PML)
N_PML = 8
EPSR = 2.2
R_IN, R_OUT = 1.12, 4.0
HALF = 6.0             # mm, half-width of the PEC box (integer lines: y = 0 is a line, as in minimal_coax_pml.py)
PROBE_DIST = 100.0     # mm, PML face -> probe
SRC_DIST = 30.0        # mm, probe -> source
TOP = 450.0            # mm, source -> z-max PEC wall (its echo comes after T_RUN)
REF_EXTRA = 300.0      # mm, the reference run's z-min extension
F0, FC = 2.5e9, 2.0e9  # 0.5 .. 4.5 GHz
T_RUN = 3.5e-9
FREQS = (0.5e9, 1.0e9, 2.0e9, 2.5e9, 3.0e9, 4.0e9)
# =============================================================================
E = os.environ.get
XY_RES = float(E('XY_RES', XY_RES))


def run(ref):
    FDTD = openEMS(EndCriteria=1e-30, MaxTime=T_RUN)
    FDTD.SetGaussExcite(F0, FC); FDTD.SetExciteZeroMean(True); FDTD.SetOverSampling(20)
    FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'PEC' if ref else 'PML_%d' % N_PML, 'PEC'])
    CSX = ContinuousStructure(); FDTD.SetCSX(CSX)
    g = CSX.GetGrid(); g.SetDeltaUnit(1e-3)
    n = int(round(2 * HALF / XY_RES))
    g.SetLines('x', np.linspace(-HALF, HALF, n + 1)); g.SetLines('y', np.linspace(-HALF, HALF, n + 1))
    z_lo = -(REF_EXTRA if ref else N_PML * DZ)
    z = np.round(np.arange(z_lo, PROBE_DIST + SRC_DIST + TOP + DZ / 2, DZ), 6)
    g.SetLines('z', z)
    metal = CSX.AddMetal('PEC'); metal.AddBox([-HALF, -HALF, z[0]], [HALF, HALF, z[-1]], priority=1)
    diel = CSX.AddMaterial('diel', epsilon=EPSR); diel.AddCylinder([0, 0, z[0]], [0, 0, z[-1]], R_OUT, priority=5)
    metal.AddCylinder([0, 0, z[0]], [0, 0, z[-1]], R_IN, priority=10)
    zi = lambda v: z[int(np.argmin(abs(z - v)))]
    zp, zs = zi(PROBE_DIST), zi(PROBE_DIST + SRC_DIST)
    xs = [v for v in np.linspace(-HALF, HALF, n + 1) if 1.2 < v < 3.9]
    CSX.AddExcitation('src', exc_type=0, exc_val=[1, 0, 0]).AddBox([xs[0], 0, zs], [xs[-1], 0, zs])
    CSX.AddProbe('E_p', p_type=2).AddBox([2.5, 0, zp], [2.5, 0, zp])
    sim = tempfile.mkdtemp(prefix='pmlrefl_')
    FDTD.Run(sim, cleanup=True, verbose=0)
    a = np.loadtxt(os.path.join(sim, 'E_p'), comments='%')
    shutil.rmtree(sim, ignore_errors=True)
    return a[:, 0], a[:, 1]


t, v = run(False)
tr, vr = run(True)
m = min(len(t), len(tr)); t, v, vr = t[:m], v[:m], vr[:m]
spec = lambda s, f: abs(np.sum(s * np.exp(-2j * np.pi * f * t)))
cfg = ' '.join('%s=%s' % (k, os.environ[k]) for k in sorted(os.environ) if k.startswith('UPML_CFS')) or 'plain UPML'
print('REFL xy %.2f  %-52s ' % (XY_RES, cfg) + '  '.join('%.1fGHz %6.1f dB' % (f / 1e9, 20 * np.log10(spec(v - vr, f) / spec(vr, f))) for f in FREQS))
