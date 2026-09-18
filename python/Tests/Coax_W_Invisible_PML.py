"""
 Coaxial line terminated by the invisible PML

 The cross-section, the mode files (Coax_Er.csv / Coax_Hr.csv) and the
 excitation are those of Coax_W_WG_Ports.py. What changes:

 * Every boundary is PEC. The shield is the outer conductor, so the side
   walls only close the box; nothing outside the shield needs absorbing.
   (A side PML in this box -- the SimBox ends at the shield's outside, so
   8 cells reach in to r ~ 0.66 mm -- damps the TEM mode in the teflon.)
 * The line runs to the z ends of the domain, and a full cross-section
   invisible PML sheet sits on each end plane. The z mesh is uniform, so
   the virtual PML cells have the size of the line's cells.

   IPML_MODE=IPML   sheets on the PEC z domain boundaries
   IPML_MODE=BLOCK  sheets on the faces of PEC blocks inside the domain
   IPML_MODE=REAL   no sheets; mesh, wire, shield and teflon run on N cells
                    past each end into openEMS' own PML_N z boundaries
 With every wall PEC, IPML and BLOCK reproduce REAL to round-off.

 Knobs:
   IPML_MODE   IPML | BLOCK | REAL   (IPML)
   IPML_N      8 | 16 | 32           (8)
   IPML_ZMEAN  zero-mean excitation  (1)
   IPML_DUMP   full E field dump     (0)
   IPML_PLOT   S-parameter figure    (1)

 The mode file is on its own 67 x 67 grid, not this mesh's lines, so openEMS
 interpolates it.

 (c) 2026 Gadi Lahav <gadi@rfwithcare.com>
"""

import os, tempfile, shutil
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

from CSXCAD.CSProperties import ABCtype

MODE = os.environ.get('IPML_MODE', 'IPML').upper()
N = int(os.environ.get('IPML_N', '8'))
assert MODE in ('IPML', 'BLOCK', 'REAL')
assert N in (8, 16, 32)
PML_TYPE = {8: ABCtype.PML_8, 16: ABCtype.PML_16, 32: ABCtype.PML_32}[N]

Sim_Path = os.path.join(tempfile.gettempdir(), 'Coax_IPML_%s_%d' % (MODE, N))
shutil.rmtree(Sim_Path, ignore_errors=True)
os.makedirs(Sim_Path)
shutil.copy("Coax_Er.csv", Sim_Path)
shutil.copy("Coax_Hr.csv", Sim_Path)

# ## Cross-section -- as in Coax_W_WG_Ports.py
coax_D = 2
coax_shield_thick = 0.15
coax_wire_D = 0.5
coax_L = 25

teflon_epsR = 2.5

unit_res = 1e-3

f0 = 1.55e9
fc = 1.45e9

# from the mode solve of Coax_Er.csv / Coax_Hr.csv
kz = 82.84554871
Zw = 238.26517157
Zl = 52.43928218

FDTD = openEMS(NrTS=300000, EndCriteria=1e-6)
FDTD.SetGaussExcite(f0, fc)
# zero time-integral excitation, so the pulse leaves no static charge at the
# port (IPML_ZMEAN=0 turns it off)
FDTD.SetExciteZeroMean(bool(int(os.environ.get('IPML_ZMEAN', '1'))))
zbc = 'PML_%d' % N if MODE == 'REAL' else 'PEC'
FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', zbc, zbc])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit_res)
mesh_res = ((C0 / (f0 + fc)) / unit_res) / 100

# uniform z mesh; the variants differ only past the ends of the line
dz = mesh_res / 2
nz = int(round(coax_L / dz))
dz = coax_L / nz
ext = {'IPML': 0, 'BLOCK': 2, 'REAL': N}[MODE]
z0, z1 = -ext * dz, coax_L + ext * dz
mesh.AddLine('z', dz * np.arange(-ext, nz + ext + 1))

r_box = coax_D * 0.5 + coax_shield_thick
mesh.AddLine('x', [-r_box, r_box])
mesh.AddLine('y', [-r_box, r_box])

# center wire
wire = CSX.AddMetal('Wire_Inner')
wire.AddCylinder(priority=10, start=[0, 0, z0], stop=[0, 0, z1], radius=coax_wire_D * 0.5)
mesh.AddLine('x', np.linspace(-coax_wire_D * 0.5, coax_wire_D * 0.5, 6))
mesh.AddLine('y', np.linspace(-coax_wire_D * 0.5, coax_wire_D * 0.5, 6))

# outer shield
shield = CSX.AddMetal('Shield_Outer')
shield.AddCylindricalShell(priority=10, start=[0, 0, z0], stop=[0, 0, z1],
                           radius=(coax_D + coax_shield_thick) * 0.5, shell_width=coax_shield_thick)
rad = (coax_D + coax_shield_thick) * 0.5
hthick = coax_shield_thick * 0.5
for ax in ('x', 'y'):
    mesh.AddLine(ax, np.linspace(-(rad + hthick), -(rad - hthick), 4).tolist() + 
                     np.linspace(rad - hthick, rad + hthick, 4).tolist())

# teflon fill
teflon = CSX.AddMaterial('PTFE', epsilon=teflon_epsR)
teflon.AddCylindricalShell(priority=8, start=[0, 0, z0], stop=[0, 0, z1],
                           radius=(coax_wire_D + coax_D) * 0.25, shell_width=(coax_D - coax_wire_D) * 0.5)
rad = (coax_wire_D + coax_D) * 0.25
hthick = (coax_D - coax_wire_D) * 0.25
for ax in ('x', 'y'):
    mesh.AddLine(ax, np.linspace(-(rad + hthick), -(rad - hthick), 12).tolist() + 
                     np.linspace(rad - hthick, rad + hthick, 12).tolist())

mesh.SmoothMeshLines('x', mesh_res, 1.25)
mesh.SmoothMeshLines('y', mesh_res, 1.25)

# ports 4 cells in from each end, as in the other invisible PML tests
zP1, zP2 = 4 * dz, coax_L - 4 * dz
port1 = FDTD.AddWaveGuidePort(1, [-r_box, -r_box, zP1], [r_box, r_box, zP1 + dz], 'z',
                              E_file="Coax_Er.csv", H_file="Coax_Hr.csv", kc=0.0, excite=1, excite_type=0)
port2 = FDTD.AddWaveGuidePort(2, [-r_box, -r_box, zP2], [r_box, r_box, zP2 - dz], 'z',
                              E_file="Coax_Er.csv", H_file="Coax_Hr.csv", kc=0.0, excite=0, excite_type=0)

if MODE == 'BLOCK':
    pec = CSX.AddMetal('PEC_blocks')
    pec.AddBox(priority=50, start=[-r_box, -r_box, z0], stop=[r_box, r_box, 0.0])
    pec.AddBox(priority=50, start=[-r_box, -r_box, coax_L], stop=[r_box, r_box, z1])

if MODE in ('IPML', 'BLOCK'):
    abs1 = CSX.AddAbsorbingBC('abs1', NormalSignPositive=True, AbsorbingBoundaryType=PML_TYPE)
    abs1.AddBox([-r_box, -r_box, 0.0], [r_box, r_box, 0.0], priority=60)
    abs2 = CSX.AddAbsorbingBC('abs2', NormalSignPositive=False, AbsorbingBoundaryType=PML_TYPE)
    abs2.AddBox([-r_box, -r_box, coax_L], [r_box, r_box, coax_L], priority=60)

if int(os.environ.get('IPML_DUMP', '0')):
    Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
    Et.AddBox([-r_box, -r_box, z0], [r_box, r_box, z1])

print('mode = %s, N = %d, dz = %.4f mm, Sim_Path = %s' % (MODE, N, dz, Sim_Path))
FDTD.Run(Sim_Path, verbose=3, cleanup=False)

# ## Post-processing -- as in Coax_W_WG_Ports.py, with the solved Zl
freq = np.linspace(max(1e9, f0 - fc), f0 + fc, 401)
port1.CalcPort(Sim_Path, freq, ref_impedance=Zw, ZL=Zl)
port2.CalcPort(Sim_Path, freq, ref_impedance=Zw, ZL=Zl)
s11 = port1.uf_ref / port1.uf_inc
s21 = port2.uf_ref / port1.uf_inc
print('\n  f/GHz   |S11| dB   |S21| dB')
for f in (1e9, 1.5e9, 2e9, 2.5e9, 3e9):
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
    plt.title('Coax, %s, N = %d' % (MODE, N))
    plt.show()
