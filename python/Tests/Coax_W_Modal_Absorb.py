"""
 Coaxial waveguide with a waveguide port on one end and a modal absorber on the other.

 Demonstrates the modal absorbing BC by replacing the second waveguide port with
 a ModalAbsorber that terminates the coaxial TEM mode using the same CSV mode files
 as the waveguide port.

 Expected result: S11 should be very low (good absorption at the far end).

 (c) 2023-2025 Gadi Lahav <gadi@rfwithcare.com>
"""

# ## Import Libraries
import os, tempfile, shutil
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

# ## General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_Coax_Modal_Absorb')

print("Copying files 'Coax_Er.csv' and 'Coax_Hr.csv' to {}".format(Sim_Path))
if not os.path.exists(Sim_Path):
    os.mkdir(Sim_Path)
shutil.copy("Coax_Er.csv", Sim_Path)
shutil.copy("Coax_Hr.csv", Sim_Path)

post_proc_only = False
display_structure = False

# substrate setup
coax_D             = 2
coax_shield_thick  = 0.15
coax_wire_D        = 0.5
coax_L             = 25

teflon_epsR = 2.5

mesh_res = 0.5

Airbox_Add = 0
unit_res = 1e-3

# simulation box
SimBox = np.array([
            -(coax_D * 0.5 + coax_shield_thick + Airbox_Add),
             (coax_D * 0.5 + coax_shield_thick + Airbox_Add),
            -(coax_D * 0.5 + coax_shield_thick + Airbox_Add),
             (coax_D * 0.5 + coax_shield_thick + Airbox_Add),
            -Airbox_Add,
             coax_L + Airbox_Add])

# excitation
f0 = 2.5e9
fc = 1e9

# ## FDTD setup
FDTD = openEMS(NrTS=300000, EndCriteria=1e-4)
FDTD.SetGaussExcite(f0, fc)
FDTD.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)

mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit_res)
mesh_res = ((C0 / (f0 + fc)) / unit_res) / 100

# initialize mesh with simulation box
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# center wire
line = CSX.AddMetal('Wire_Inner')
start = [0.0, 0.0, 0.0]
stop  = [0.0, 0.0, coax_L]
line.AddCylinder(priority=10, start=start, stop=stop, radius=coax_wire_D * 0.5)
mesh.AddLine('x', np.linspace(start[0] - coax_wire_D * 0.5, stop[0] + coax_wire_D * 0.5, 6).tolist())
mesh.AddLine('y', np.linspace(start[1] - coax_wire_D * 0.5, stop[1] + coax_wire_D * 0.5, 6).tolist())
mesh.AddLine('z', [start[2], stop[2]])

# outer shield
shield = CSX.AddMetal('Shield_Outer')
start = [0.0, 0.0, 0.0]
stop  = [0.0, 0.0, coax_L]
shield.AddCylindricalShell(priority=10, start=start, stop=stop,
                            radius=(coax_D + coax_shield_thick) * 0.5,
                            shell_width=coax_shield_thick)
rad    = (coax_D + coax_shield_thick) * 0.5
hthick = coax_shield_thick * 0.5
mesh.AddLine('x',
             np.linspace(start[0] - (rad + hthick), start[0] - (rad - hthick), 4).tolist() +
             np.linspace(stop[0]  + (rad - hthick), stop[0]  + (rad + hthick), 4).tolist())
mesh.AddLine('y',
             np.linspace(start[1] - (rad + hthick), start[1] - (rad - hthick), 4).tolist() +
             np.linspace(stop[1]  + (rad - hthick), stop[1]  + (rad + hthick), 4).tolist())
mesh.AddLine('z', [start[2], stop[2]])

# teflon fill
teflon = CSX.AddMaterial('PTFE', epsilon=teflon_epsR)
start = [0.0, 0.0, 0.0]
stop  = [0.0, 0.0, coax_L]
teflon.AddCylindricalShell(priority=8, start=start, stop=stop,
                            radius=(coax_wire_D + coax_D) * 0.25,
                            shell_width=(coax_D - coax_wire_D) * 0.5)
rad    = (coax_wire_D + coax_D) * 0.25
hthick = (coax_D - coax_wire_D) * 0.25
mesh.AddLine('x',
             np.linspace(start[0] - (rad + hthick), start[0] - (rad - hthick), 12).tolist() +
             np.linspace(stop[0]  + (rad - hthick), stop[0]  + (rad + hthick), 12).tolist())
mesh.AddLine('y',
             np.linspace(start[1] - (rad + hthick), start[1] - (rad - hthick), 12).tolist() +
             np.linspace(stop[1]  + (rad - hthick), stop[1]  + (rad + hthick), 12).tolist())
mesh.AddLine('z', [start[2], stop[2]])

# dense mesh near ports
mesh.AddLine('z', np.array([0.25, 0.5, 0.8, 1]) * mesh_res)
mesh.AddLine('z', coax_L - np.array([0.25, 0.5, 0.8, 1]) * mesh_res)

mesh.SmoothMeshLines('all', mesh_res, 1.25)

# find port / absorber mesh positions
Zz = mesh.GetLines('z')
idxPort1 = (np.where(Zz == 0.0)[0] + 1).item(0)
idxPort2 = (np.where(Zz == coax_L)[0] - 1).item(0)
idxAbs2  = idxPort2 + 1   # one cell outside the second port plane → absorber location

# --- Port 1: waveguide port with excitation (same as Coax_W_WG_Ports.py) ---
start = [-coax_D * 0.5 - coax_shield_thick, -coax_D * 0.5 - coax_shield_thick, Zz.item(idxPort1 + 0)]
stop  = [ coax_D * 0.5 + coax_shield_thick,  coax_D * 0.5 + coax_shield_thick, Zz.item(idxPort1 + 1)]
port1 = FDTD.AddWaveGuidePort(1, start, stop, 'z',
                               E_file="Coax_Er.csv", H_file="Coax_Hr.csv",
                               kc=0.0, excite=1, excite_type=0)

# --- Absorber: modal absorber replacing the second waveguide port ---
# The sheet is placed at the mesh line just outside the waveguide end.
# normal_positive=False because the incoming wave travels in the +z direction
# and the absorber faces it from the far (high-z) end.
abs_z = Zz.item(idxAbs2)
abs_start = [-coax_D * 0.5 - coax_shield_thick, -coax_D * 0.5 - coax_shield_thick, abs_z]
abs_stop  = [ coax_D * 0.5 + coax_shield_thick,  coax_D * 0.5 + coax_shield_thick, abs_z]
modal_abs = FDTD.AddModalAbsorber(abs_start, abs_stop, 'z',
                                   E_file="Coax_Er.csv",
                                   H_file="Coax_Hr.csv",
                                   normal_positive=False,
                                   Zw=238.26517157)

# ## Run the simulation
if display_structure:
    CSX_file = os.path.join(Sim_Path, 'coax_modal_absorb.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)

    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=0, cleanup=False)

# ## Post-processing
# Only port 1 is available; measure reflection to assess absorber quality.
f = np.linspace(max(1e9, f0 - fc), f0 + fc, 401)

Zw  = np.array([238.26517157])   # waveguide impedance (same as Coax_W_WG_Ports reference)
port1.CalcPort(Sim_Path, f, ref_impedance=Zw, ZL=50)

s11     = port1.uf_ref / port1.uf_inc
s11_dB  = 20.0 * np.log10(np.abs(s11))

figure()
plot(f / 1e9, s11_dB, 'k-', linewidth=2, label='$S_{11}$ (reflection, modal absorber at far end)')
grid()
legend()
ylabel('S-Parameter (dB)')
xlabel('Frequency (GHz)')
title('Coaxial line — modal absorber vs. WG port termination')

show()
