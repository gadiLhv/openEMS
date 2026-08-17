"""
 Coaxial waveguide with a waveguide port on one end and a modal absorber on the other.

 Demonstrates the modal absorbing BC by replacing the second waveguide port with
 a ModalAbsorber that terminates the coaxial TEM mode using the same CSV mode files
 as the waveguide port.

 Expected result: S11 should be very low (good absorption at the far end).

 (c) 2023-2025 Gadi Lahav <gadi@rfwithcare.com>
"""

# ## Import Libraries
import os, tempfile, shutil, glob, re
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

Airbox_Add = 1
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
FDTD = openEMS(NrTS=300000, EndCriteria=1e-4, OverSampling=4)
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
idxPort1 = (np.where(Zz == 0.0)[0] + 15).item(0)
idxPort2 = (np.where(Zz == coax_L)[0] - 3).item(0)
idxAbs1  = idxPort1 - 10   # one cell outside the second port plane → absorber location
idxAbs2  = idxPort2 + 1   # one cell outside the second port plane → absorber location

# --- Port 1: waveguide port with excitation (same as Coax_W_WG_Ports.py) ---
start = [-coax_D * 0.5 - coax_shield_thick, -coax_D * 0.5 - coax_shield_thick, Zz.item(idxPort1 + 0)]
stop  = [ coax_D * 0.5 + coax_shield_thick,  coax_D * 0.5 + coax_shield_thick, Zz.item(idxPort1 + 1)]
port1 = FDTD.AddWaveGuidePort(1, start, stop, 'z',
                               E_file="Coax_Er.csv", H_file="Coax_Hr.csv",
                               kc=0.0, excite=1, excite_type=0)

abs_z = Zz.item(idxAbs1)
abs_start = [-coax_D * 0.5 - coax_shield_thick, -coax_D * 0.5 - coax_shield_thick, abs_z]
abs_stop  = [ coax_D * 0.5 + coax_shield_thick,  coax_D * 0.5 + coax_shield_thick, abs_z]
modal_abs_1 = FDTD.AddModalAbsorber(abs_start, abs_stop, 'z',
                                    E_file="Coax_Er.csv",
                                    H_file="Coax_Hr.csv",
                                    mode_type='TEM',   # no cutoff -> scalar Zw
                                    normal_positive=True,
                                    Zw=238.26517157)


# --- Port 2: passive waveguide port, for the S21 measurement ---
# Reversed z ordering so the port faces the incoming +z wave, exactly as in
# RectWG_W_ModalMur.py. It sits one cell inside absorber 2.
start = [-coax_D * 0.5 - coax_shield_thick, -coax_D * 0.5 - coax_shield_thick, Zz.item(idxPort2 - 0)]
stop  = [ coax_D * 0.5 + coax_shield_thick,  coax_D * 0.5 + coax_shield_thick, Zz.item(idxPort2 - 1)]
port2 = FDTD.AddWaveGuidePort(2, start, stop, 'z',
                               E_file="Coax_Er.csv", H_file="Coax_Hr.csv",
                               kc=0.0, excite=0, excite_type=0)

# --- Absorber 2: terminates the far (high-z) end, behind port 2 ---
# normal_positive=False because the incoming wave travels in the +z direction
# and the absorber faces it from the far (high-z) end.
abs_z = Zz.item(idxAbs2)
abs_start = [-coax_D * 0.5 - coax_shield_thick, -coax_D * 0.5 - coax_shield_thick, abs_z]
abs_stop  = [ coax_D * 0.5 + coax_shield_thick,  coax_D * 0.5 + coax_shield_thick, abs_z]
modal_abs_2 = FDTD.AddModalAbsorber(abs_start, abs_stop, 'z',
                                    E_file="Coax_Er.csv",
                                    H_file="Coax_Hr.csv",
                                    mode_type='TEM',   # no cutoff -> scalar Zw
                                    normal_positive=False,
                                    Zw=238.26517157)

### Field export -- disabled. Uncomment to dump E(t) over the whole box.
# Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
# start = [SimBox[0], SimBox[2], SimBox[4]]
# stop  = [SimBox[1], SimBox[3], SimBox[5]]
# Et.AddBox(start, stop)

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
# Port 1 gives the reflection (absorber quality), port 2 the through wave.
f = np.linspace(max(1e9, f0 - fc), f0 + fc, 401)

Zw  = np.array([238.26517157])   # modal impedance (same as Coax_W_WG_Ports)
port1.CalcPort(Sim_Path, f, ref_impedance=Zw, ZL=50)
port2.CalcPort(Sim_Path, f, ref_impedance=Zw, ZL=50)

s11 = port1.uf_ref / port1.uf_inc
s21 = port2.uf_ref / port1.uf_inc
s11_dB = 20.0 * np.log10(np.abs(s11)).ravel()
s21_dB = 20.0 * np.log10(np.abs(s21)).ravel()

print('\n===== Coax TEM, modal absorbers (scalar Zw, mode_type=TEM) =====')
print('  f[GHz]    |S11| dB    |S21| dB')
for fi in [1.5, 2.0, 2.5, 3.0, 3.5]:
    if fi * 1e9 < f[0] or fi * 1e9 > f[-1]:
        continue
    j = np.abs(f - fi * 1e9).argmin()
    print('   {:.2f}    {:8.2f}    {:8.2f}'.format(fi, s11_dB[j], s21_dB[j]))
print('worst |S11| in band: {:.2f} dB'.format(np.max(s11_dB)))

# Energy balance sanity check. In a lossless line |S11|^2 + |S21|^2 must be <= 1.
# A value meaningfully above 1 means the port de-embedding is not separating
# incident from reflected consistently, and the |S11| column above is then
# reporting the de-embedding rather than the absorber. This is a property of
# WaveguidePort.CalcPort, not of the termination.
bal = np.abs(s11).ravel() ** 2 + np.abs(s21).ravel() ** 2
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
title('Coaxial line (TEM) — modal absorbers at both ends')

# ## Modal-absorber mode-match graphical debugging
# Each modal absorber dumps two time-domain mode-match files into Sim_Path:
#   modal_absorber_<i>_E  -> column 'voltage' = E-field mode amplitude (E)
#   modal_absorber_<i>_H  -> column 'current' = H-field mode amplitude (H)
# (file columns:  % t/s   <voltage|current>   mode_purity)
# The travelling-wave amplitudes are  a+- = E +- Zw*H, i.e. the forward/backward
# mode amplitudes. For a working absorber the incident wave should carry almost
# all of the energy in one of a+/a-, and the counter-propagating component
# (the reflection the absorber failed to swallow) should be strongly suppressed.

# Wave impedance used when each absorber was created (see AddModalAbsorber above).
# Map absorber index -> Zw if they ever differ; here both share the same value.
Zw_abs_default = 238.26517157
Zw_abs = {}   # e.g. {0: 238.26517157} to override per absorber

def _read_modematch_td(fname):
    """Return (t, value) from a ProcessModeMatch time-domain dump (data column 1)."""
    data = np.loadtxt(fname, comments='%', ndmin=2)
    if data.size == 0:
        return None, None
    return data[:, 0], data[:, 1]

_idx_re = re.compile(r'modal_absorber_(\d+)_E$')
e_files = sorted(glob.glob(os.path.join(Sim_Path, 'modal_absorber_*_E')),
                 key=lambda p: int(_idx_re.search(p).group(1)))

if not e_files:
    print("No modal-absorber mode-match files found in {} "
          "(run the simulation first).".format(Sim_Path))

for e_file in e_files:
    idx    = int(_idx_re.search(e_file).group(1))
    h_file = os.path.join(Sim_Path, 'modal_absorber_{}_H'.format(idx))
    if not os.path.exists(h_file):
        print("Absorber {}: missing H file {}, skipping.".format(idx, h_file))
        continue

    tE, E = _read_modematch_td(e_file)
    tH, H = _read_modematch_td(h_file)
    if tE is None or tH is None:
        print("Absorber {}: empty mode-match file, skipping.".format(idx))
        continue

    # H lives on the dual (half-timestep-shifted) Yee grid -> resample onto tE
    # so E and H can be combined sample-by-sample.
    H_on_E = np.interp(tE, tH, H)

    Zw_i    = Zw_abs.get(idx, Zw_abs_default)
    a_plus  = E + Zw_i * H_on_E
    a_minus = E - Zw_i * H_on_E

    figure()
    subplot(2, 1, 1)
    plot(tE / 1e-9, E,            'b-',  label='E  (voltage mode match)')
    plot(tE / 1e-9, Zw_i * H_on_E, 'r--', label=r'$Z_w\,H$  (current mode match)')
    grid(); legend(); ylabel('amplitude')
    title('Modal absorber #{}  —  mode-match constituents ($Z_w$ = {:.3f} $\\Omega$)'.format(idx, Zw_i))

    subplot(2, 1, 2)
    plot(tE / 1e-9, a_plus,  'k-', linewidth=2, label=r'$a_+ = E + Z_w H$')
    plot(tE / 1e-9, a_minus, 'g-', linewidth=2, label=r'$a_- = E - Z_w H$')
    grid(); legend(); ylabel('amplitude'); xlabel('time (ns)')

show()
