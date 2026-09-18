"""
 Rectangular waveguide terminated by the invisible PML

 The invisible PML (ABCtype.PML_8/16/32) is an N-cell UPML that lives only in
 the engine extension, extruded behind a sheet that lies on a PEC face. This
 test runs the same physical waveguide [0, length] three ways:

   IPML_MODE=IPML   sheets on the PEC domain boundaries z=0 and z=length
   IPML_MODE=BLOCK  sheets on the faces of PEC blocks inside the domain
                    (the domain runs on 2 cells behind each block face)
   IPML_MODE=REAL   no sheets; the mesh runs on N cells past each end and the
                    z boundaries are openEMS' own PML_N

 The invisible PML uses the equations of the real one, so IPML and BLOCK must
 reproduce REAL to floating-point round-off. Run all three and compare the
 data lines of the port_* files in the Sim_Path folders; they are identical.

 Knobs (IPML_PLOT=0 skips the S-parameter figure, IPML_ZMEAN=0 the zero-mean excitation):
   IPML_MODE  IPML | BLOCK | REAL  (IPML)
   IPML_N     8 | 16 | 32          (8)
   IPML_GAP   in BLOCK mode, the number of air cells between the PEC block face
              and the sheet (0). The sheet owns every field that couples its
              two sides, so a gap must change nothing.

 (c) 2026 Gadi Lahav <gadi@rfwithcare.com>
"""

import os, tempfile, shutil
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

from CSXCAD.CSProperties import ABCtype

MODE = os.environ.get('IPML_MODE', 'IPML').upper()
N = int(os.environ.get('IPML_N', '8'))
GAP = int(os.environ.get('IPML_GAP', '0'))
assert MODE in ('IPML', 'BLOCK', 'REAL')
assert N in (8, 16, 32)
PML_TYPE = {8: ABCtype.PML_8, 16: ABCtype.PML_16, 32: ABCtype.PML_32}[N]

Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WG_IPML_%s_%d%s' % (MODE, N, '_gap%d' % GAP if GAP else ''))
shutil.rmtree(Sim_Path, ignore_errors=True)
os.makedirs(Sim_Path)

unit = 1e-6  # drawing unit in um

# WR42
a = 10700
b = 4300
length = 50000

f_start = 20e9
f_0 = 24e9
f_stop = 26e9
lambda0 = C0 / f_0 / unit

TE_mode = 'TE10'

# a uniform z mesh, so every variant shares the physical region cell by cell
dz = 250.0
nz = int(round(length / dz))
length = nz * dz
mesh_res = lambda0 / 50

FDTD = openEMS(NrTS=20000, EndCriteria=1e-6)
FDTD.SetGaussExcite(0.5 * (f_start + f_stop), 0.5 * (f_stop - f_start))
# zero time-integral excitation, so the pulse leaves no static charge at the
# port (IPML_ZMEAN=0 turns it off)
FDTD.SetExciteZeroMean(bool(int(os.environ.get('IPML_ZMEAN', '1'))))

zbc = 'PML_%d' % N if MODE == 'REAL' else 'PEC'
FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', zbc, zbc])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.AddLine('x', [0, a])
mesh.AddLine('y', [0, b])
mesh.SmoothMeshLines('x', mesh_res, ratio=1.4)
mesh.SmoothMeshLines('y', mesh_res, ratio=1.4)

ext = {'IPML': 0, 'BLOCK': 2 + GAP, 'REAL': N}[MODE]
mesh.AddLine('z', dz * np.arange(-ext, nz + ext + 1))

ports = []
ports.append(FDTD.AddRectWaveGuidePort(0, [0, 0, 3 * dz], [a, b, 4 * dz], 'z', a * unit, b * unit, TE_mode, 1))
ports.append(FDTD.AddRectWaveGuidePort(1, [0, 0, length - 3 * dz], [a, b, length - 4 * dz], 'z', a * unit, b * unit, TE_mode))

if MODE == 'BLOCK':
    pec = CSX.AddMetal('PEC_blocks')
    pec.AddBox(priority=5, start=[0, 0, -ext * dz], stop=[a, b, -GAP * dz])
    pec.AddBox(priority=5, start=[0, 0, length + GAP * dz], stop=[a, b, length + ext * dz])

if MODE in ('IPML', 'BLOCK'):
    abs1 = CSX.AddAbsorbingBC('abs1', NormalSignPositive=True, AbsorbingBoundaryType=PML_TYPE)
    abs1.AddBox([0, 0, 0], [a, b, 0], priority=6)
    abs2 = CSX.AddAbsorbingBC('abs2', NormalSignPositive=False, AbsorbingBoundaryType=PML_TYPE)
    abs2.AddBox([0, 0, length], [a, b, length], priority=6)

print('mode = %s, N = %d, Sim_Path = %s' % (MODE, N, Sim_Path))
FDTD.Run(Sim_Path, verbose=3, cleanup=False)

freq = np.linspace(f_start, f_stop, 201)
for port in ports:
    port.CalcPort(Sim_Path, freq)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc
print('\n  f/GHz   |S11| dB   |S21| dB')
for f in (20e9, 21.5e9, 23e9, 24.5e9, 26e9):
    k = np.argmin(abs(freq - f))
    print('  %5.2f   %8.2f   %8.3f' % (freq[k] / 1e9, 20 * np.log10(abs(s11[k])), 20 * np.log10(abs(s21[k]))))
print('  worst |S11| over the band: %.2f dB' % (20 * np.log10(abs(s11).max())))

# ## Plot the S-parameters (IPML_PLOT=0 skips the figure, for batch runs)
if int(os.environ.get('IPML_PLOT', '1')):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(freq / 1e9, 20 * np.log10(np.abs(s11)), 'k-', linewidth=2, label='$S_{11}$')
    plt.plot(freq / 1e9, 20 * np.log10(np.abs(s21)), 'b-', linewidth=2, label='$S_{21}$')
    plt.grid()
    plt.legend()
    plt.ylabel('S-Parameter (dB)')
    plt.xlabel('Frequency (GHz)')
    plt.title('Rect. waveguide, %s, N = %d' % (MODE, N))
    plt.show()
