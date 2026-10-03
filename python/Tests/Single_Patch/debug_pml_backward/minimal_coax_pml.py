"""Minimal repro: a short coax cavity terminated by openEMS's own boundary UPML.

No IPML, no STEP, no ports. An ideal coax (pin r 1.12, outer r 4.0, Dk 2.2) in a
12 x 12 mm PEC box, the domain's PML_8 at z-min, a PEC short at z-max. A soft
E source in the annulus, an E probe beside it. The report: the late-time growth
rate of |E| (dB/ns); positive = unstable.
"""
import os, sys, shutil, tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

# ============================== CONFIG =======================================
COAX_LEN = 1.225       # mm, PML face -> PEC short
XY_RES = 0.5           # mm, transverse cells (1.0 is unstable: backward-wave modes)
DZ = 0.245             # mm, cells along the coax (and in the PML)
N_PML = 8
EPSR = 2.2
R_IN, R_OUT = 1.12, 4.0
HALF = 6.0             # mm, half-width of the PEC box
NRTS = 20000
TSM = 3                # timestep method: 3 Rennings (default), 1 CFL
# =============================================================================
E = os.environ.get     # sweep overrides (harness use only)
COAX_LEN = float(E('COAX_LEN', COAX_LEN)); XY_RES = float(E('XY_RES', XY_RES)); DZ = float(E('DZ', DZ))
N_PML = int(E('N_PML', N_PML)); NRTS = int(E('NRTS', NRTS)); TSM = int(E('TSM', TSM)); EPSR = float(E('EPSR', EPSR))

FDTD = openEMS(EndCriteria=1e-30, NrTS=NRTS)
FDTD.SetTimeStepMethod(TSM)
FDTD.SetGaussExcite(2.5e9, 1.0e9)
FDTD.SetExciteZeroMean(True)
FDTD.SetOverSampling(50)
FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'PML_%d' % N_PML, 'PEC'])
CSX = ContinuousStructure(); FDTD.SetCSX(CSX)
g = CSX.GetGrid(); g.SetDeltaUnit(1e-3)
n = int(round(2 * HALF / XY_RES))
g.SetLines('x', np.linspace(-HALF, HALF, n + 1)); g.SetLines('y', np.linspace(-HALF, HALF, n + 1))
nz = int(round(COAX_LEN / DZ))
z = np.concatenate([-DZ * np.arange(N_PML, 0, -1), np.linspace(0, COAX_LEN, max(nz, 1) + 1)])
g.SetLines('z', z)

metal = CSX.AddMetal('PEC')
metal.AddBox([-HALF, -HALF, z[0]], [HALF, HALF, z[-1]], priority=1)
diel = CSX.AddMaterial('diel', epsilon=EPSR)
SHAPE = E('SHAPE', 'coax')   # coax | tube (no pin) | rect (dielectric-filled square guide, no curves)
if SHAPE == 'rect':
    diel.AddBox([-R_OUT, -R_OUT, z[0]], [R_OUT, R_OUT, z[-1]], priority=5)
else:
    diel.AddCylinder([0, 0, z[0]], [0, 0, z[-1]], R_OUT, priority=5)
if SHAPE == 'coax':
    metal.AddCylinder([0, 0, z[0]], [0, 0, z[-1]], R_IN, priority=10)
zs = z[N_PML + max(nz // 2, 0)]
ex = CSX.AddExcitation('src', exc_type=0, exc_val=[1, 0, 0]); ex.AddBox([2.0, 0, zs], [3.0, 0, zs])
CSX.AddProbe('E_p', p_type=2).AddBox([2.5, 1.0, zs], [2.5, 1.0, zs])

sim = tempfile.mkdtemp(prefix='mincoax_')
FDTD.Run(sim, cleanup=True, verbose=0)
a = np.loadtxt(os.path.join(sim, 'E_p'), comments='%'); t = a[:, 0]; e = np.linalg.norm(a[:, 1:4], axis=1)
w = (t > 2.5e-9) & (e > 0) & np.isfinite(e)
tt, ee = t[w], e[w]; nb = max(len(tt) // 10, 1)
et = [tt[i:i + nb].mean() for i in range(0, len(tt) - nb + 1, nb)]
em = [ee[i:i + nb].max() for i in range(0, len(tt) - nb + 1, nb)]
g_ = np.polyfit(et, np.log(em), 1)[0]; rate = 8.686 * g_ * 1e-9
_comp = np.argmax(np.abs(a[w, 1:4]).max(0)); y = a[w, 1 + _comp] * np.exp(-g_ * tt); y = y - y.mean()
_F = np.fft.rfftfreq(len(y) * 4, np.median(np.diff(tt))); f_peak = _F[np.argmax(np.abs(np.fft.rfft(y * np.hanning(len(y)), len(y) * 4)))]
print('RESULT %-5s len' % SHAPE + ' %.3f mm  xy %.3f  dz %.3f  N_PML %d  TSM %d  epsr %.2f  -> %+.1f dB/ns at %.1f GHz  (|E| end %.3g, %.2f ns)'
      % (COAX_LEN, XY_RES, DZ, N_PML, TSM, EPSR, rate, f_peak / 1e9, e[-1], t[-1] * 1e9))
shutil.rmtree(sim, ignore_errors=True)
