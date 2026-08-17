"""
 Coplanar waveguide terminated by modal absorbers, with S11 and S21.

 NOT YET VALIDATED -- read this before trusting any number it prints.

 CPW is the hardest of the three test cases in this family, and not because of
 the absorber:

  * It is an OPEN structure. The modal absorber terminates the guided mode;
    anything that radiates off the line still has to be caught by the airbox
    boundary. So |S11| cannot go arbitrarily low -- radiation and the accuracy
    of the mode file set a floor that has nothing to do with the termination.

  * There is METAL INSIDE THE PORT APERTURE (the centre strip and both ground
    planes). That matters: a mode file is not obliged to be zero inside
    conductors, and the coax CSV in this directory is not. The operator now
    masks the template onto non-conductor edges, which is what stopped the coax
    from diverging; CPW exercises the same path with a more complicated
    conductor layout.

  * Quasi-TEM. There is no cutoff, so MODE = 'TEM' below takes the scalar Zw
    absorber. In an earlier run on this geometry the scalar absorber STALLED --
    total energy flattened at about -3.5 dB and stopped decaying, where the
    MUR_1ST sheets of CPW_With_WG_Ports.py converge to -38.6 dB in a third of
    the timesteps. That was measured BEFORE the conductor-mask fix, so it is
    worth re-checking, but do not assume this script passes.

 Set MODE = 'TE' to force the dispersive Modal Mur path instead (with fc = 0,
 since a quasi-TEM line has no cutoff -- kc = 0 makes the evanescent branch
 inert and leaves a pure dispersive delay). On the PTFE coax that path measured
 7-22 dB better than the scalar one, so it is worth trying here too.

 PLACEMENT: the propagation faces (y) are PEC and the absorber sheets sit
 directly on them. Modal Mur requires that -- it is an overwrite, so anything
 behind the sheet is decoupled and a gap becomes a sealed cavity. The other
 four faces stay MUR so the structure can still radiate.

 (c) 2023-2026 Gadi Lahav <gadi@rfwithcare.com>
"""

# ## Import Libraries
import os, tempfile, shutil
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

# 'TEM' -> scalar Zw absorber.  'TE' -> dispersive Modal Mur (fc = 0).
MODE = 'TEM'

# ## General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_CPW_Modal_Absorb')
if not os.path.exists(Sim_Path):
    os.mkdir(Sim_Path)
shutil.copy("CPW_E.csv", Sim_Path)
shutil.copy("CPW_H.csv", Sim_Path)

post_proc_only = False
display_structure = False

# CPW geometry (drawing units = mm)
Line_W = 1.0
CPW_gap = 0.3
cu_thick = 0.05

substrate_epsR = 4.3
substrate_width = 11
substrate_length = 50
substrate_thickness = 1
substrate_cells = 4

port_w_fact = 1.0
port_h_fact = 2.0
portThick_mm = 0.0

Airbox_Add = 5
unit = 1e-3

SimBox = np.array([
    -substrate_width * 0.5 - Airbox_Add,
     substrate_width * 0.5 + Airbox_Add,
    -Airbox_Add,
     substrate_length + Airbox_Add,
    -substrate_thickness * (port_h_fact - 1.0) - Airbox_Add,
     substrate_thickness * (1.0 + port_h_fact) + Airbox_Add])

# excitation
f0 = 2e9
fc = 1e9

# Modal parameters of the CPW mode (same values CPW_With_WG_Ports.py uses).
beta = 62.5                       # 1/m, from a mode solver
v_phase = 2 * pi * f0 / beta
Zw_mode = 254                     # modal wave impedance, Ohm
Zl = 45                           # line impedance for the port de-embed

FDTD = openEMS(NrTS=60000, EndCriteria=1e-5)
FDTD.SetGaussExcite(f0, fc)

# PEC on the propagation faces (y), so the absorber sheets sit directly on a
# PEC block as Modal Mur requires. The rest stays MUR: a CPW radiates, and that
# energy is not the modal absorber's to remove.
FDTD.SetBoundaryCond(['MUR', 'MUR', 'PEC', 'PEC', 'MUR', 'MUR'])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)
mesh_res = (C0 / (f0 + fc) / unit) / 40

mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', [0.0, substrate_length])
mesh.AddLine('z', SimBox[4:6])

# substrate
sub = CSX.AddMaterial('FR4', epsilon=substrate_epsR)
sub.AddBox(priority=1,
           start=[-substrate_width / 2, 0.0, 0.0],
           stop=[substrate_width / 2, substrate_length, substrate_thickness])
mesh.AddLine('x', [-substrate_width / 2, substrate_width / 2])
mesh.AddLine('z', linspace(0, substrate_thickness, substrate_cells + 1))

# centre strip
line = CSX.AddMetal('CPW_line')
line.AddBox(priority=10,
            start=[-Line_W / 2, 0.0, substrate_thickness],
            stop=[Line_W / 2, substrate_length, substrate_thickness + cu_thick])

# ground planes
gnd = CSX.AddMetal('CPW_gnd')
gnd.AddBox(priority=10,
           start=[-0.5 * substrate_width, 0.0, substrate_thickness],
           stop=[-0.5 * (Line_W + CPW_gap * 2), substrate_length, substrate_thickness + cu_thick])
gnd.AddBox(priority=10,
           start=[0.5 * (Line_W + CPW_gap * 2), 0.0, substrate_thickness],
           stop=[0.5 * substrate_width, substrate_length, substrate_thickness + cu_thick])

gap_cells = 4
mesh.AddLine('x', linspace(-0.5 * (Line_W + CPW_gap * 2), -Line_W * 0.5, gap_cells))
mesh.AddLine('x', linspace(Line_W * 0.5, 0.5 * (Line_W + CPW_gap * 2), gap_cells))
mesh.AddLine('x', [-Line_W / 2, Line_W / 2])
mesh.AddLine('z', [substrate_thickness, substrate_thickness + cu_thick])

# Uniform y mesh: the delay filter is built from a lattice dispersion relation
# that assumes a locally uniform cell along the propagation direction.
dy = mesh_res
mesh.AddLine('y', np.arange(0.0, substrate_length + dy / 2, dy).tolist())
mesh.SmoothMeshLines('x', mesh_res, 1.4)
mesh.SmoothMeshLines('z', mesh_res, 1.4)

Yz = mesh.GetLines('y')

# Absorbers ON the PEC end faces (first and last y planes); ports well inside.
idxAbs1 = 0
idxAbs2 = len(Yz) - 1
idxPort1 = 10
idxPort2 = len(Yz) - 1 - 10

# Port / absorber cross-section: the same box for both, so the mode file's
# coordinate frame (anchored at the box's start corner) is shared.
box_lo = [-Line_W * (1.0 + port_w_fact) * 0.5, -substrate_thickness * (port_h_fact - 1.0)]
box_hi = [ Line_W * (1.0 + port_w_fact) * 0.5,  substrate_thickness * (1.0 + port_h_fact)]


def add_absorber(y, normal_positive):
    """Absorber sheet at plane y, method chosen by MODE."""
    start = [box_lo[0], y, box_lo[1]]
    stop  = [box_hi[0], y, box_hi[1]]
    if MODE == 'TEM':
        return FDTD.AddModalAbsorber(start, stop, 'y',
                                     E_file="CPW_E.csv", H_file="CPW_H.csv",
                                     mode_type='TEM',
                                     Zw=Zw_mode,
                                     phase_velocity=v_phase,
                                     normal_positive=normal_positive,
                                     priority=20)
    # Quasi-TEM has no cutoff, so fc = 0: kc = 0 leaves the evanescent branch
    # inert and the filter is a pure dispersive delay.
    return FDTD.AddModalAbsorber(start, stop, 'y',
                                 E_file="CPW_E.csv",
                                 mode_type='TE', fc=0.0,
                                 phase_velocity=v_phase,
                                 normal_positive=normal_positive,
                                 priority=20)


# --- Absorber 1: on the low-y PEC face ---------------------------------------
abs1 = add_absorber(Yz.item(idxAbs1), True)

# --- Port 1: excited ---------------------------------------------------------
start = [box_lo[0], Yz.item(idxPort1 + 0), box_lo[1]]
stop  = [box_hi[0], Yz.item(idxPort1 + 1), box_hi[1]]
port1 = FDTD.AddWaveGuidePort(1, start, stop, 'y',
                              E_file="CPW_E.csv", H_file="CPW_H.csv",
                              kc=0.0, excite=1, excite_type=0)

# --- Port 2: passive, for S21 ------------------------------------------------
start = [box_lo[0], Yz.item(idxPort2 - 0), box_lo[1]]
stop  = [box_hi[0], Yz.item(idxPort2 - 1), box_hi[1]]
port2 = FDTD.AddWaveGuidePort(2, start, stop, 'y',
                              E_file="CPW_E.csv", H_file="CPW_H.csv",
                              kc=0.0, excite=0, excite_type=0)

# --- Absorber 2: on the high-y PEC face --------------------------------------
abs2 = add_absorber(Yz.item(idxAbs2), False)

### Field export -- disabled. Uncomment to dump E(t) over the whole box.
# Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
# Et.AddBox([SimBox[0], SimBox[2], SimBox[4]], [SimBox[1], SimBox[3], SimBox[5]])

if display_structure:
    CSX_file = os.path.join(Sim_Path, 'cpw_modal_absorb.xml')
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=False)

# ## Post-processing
f = np.linspace(max(0.5e9, f0 - fc), f0 + fc, 401)

port1.CalcPort(Sim_Path, f, ref_impedance=Zw_mode, ZL=Zl)
port2.CalcPort(Sim_Path, f, ref_impedance=Zw_mode, ZL=Zl)

s11 = np.asarray(port1.uf_ref / port1.uf_inc).ravel()
s21 = np.asarray(port2.uf_ref / port1.uf_inc).ravel()
s11_dB = 20.0 * np.log10(np.abs(s11))
s21_dB = 20.0 * np.log10(np.abs(s21))

print('\n===== CPW, modal absorbers (MODE = {}) ====='.format(MODE))
print('  f[GHz]    |S11| dB    |S21| dB')
for fi in [1.0, 1.5, 2.0, 2.5, 3.0]:
    if fi * 1e9 < f[0] or fi * 1e9 > f[-1]:
        continue
    j = np.abs(f - fi * 1e9).argmin()
    print('   {:.2f}    {:8.2f}    {:8.2f}'.format(fi, s11_dB[j], s21_dB[j]))
print('worst |S11| in band: {:.2f} dB'.format(np.max(s11_dB)))

# Energy balance. In a CLOSED lossless two-port |S11|^2+|S21|^2 <= 1; a CPW is
# open, so genuine radiation legitimately pushes the sum BELOW 1. A sum ABOVE 1
# is the diagnostic worth watching: it means the port de-embedding is not
# separating incident from reflected consistently, and |S11| is then reporting
# the de-embedding rather than the absorber.
bal = np.abs(s11) ** 2 + np.abs(s21) ** 2
print('energy balance |S11|^2+|S21|^2: min {:.3f}  max {:.3f}{}'.format(
      bal.min(), bal.max(),
      '   <-- ABOVE 1: |S11| is de-embed limited, not the absorber' if bal.max() > 1.02 else ''))

figure()
plot(f / 1e9, s11_dB, 'k-', linewidth=2, label='$S_{11}$ (reflection)')
plot(f / 1e9, s21_dB, 'b-', linewidth=2, label='$S_{21}$ (through)')
grid()
legend()
ylabel('S-Parameter (dB)')
xlabel('Frequency (GHz)')
title('CPW — modal absorbers on both PEC faces (MODE = {})'.format(MODE))
show()
