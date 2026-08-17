"""
 Rectangular waveguide (TE10) terminated by Dispersive Modal Mur absorbers.

 The C++ twin of the Octave testbed WG_ModalMur_BoundaryTermination.m: same
 cross-section (fc = 1.2 GHz), same 8 mm cells, absorbers sitting directly on
 the PEC end faces of the guide.

 THE METHOD
 ----------
 mode_type='TE' selects it: a TE mode has a real cutoff, so no scalar wave
 impedance can work near it. (A TEM line has no cutoff, Zw is constant, and
 mode_type='TEM' takes the cheaper scalar absorber instead -- see the coax and
 CPW tests.)

 The sheet's modal component is overwritten every timestep with the delayed
 modal amplitude one cell inside,

     a_sheet(w) = a_inside(w) * exp(-j*beta(w)*dz)

 with beta from the exact discrete lattice dispersion relation. |exp(-j*beta*dz)|
 never exceeds 1, so the termination is passive by construction -- no wave
 impedance, no drain cap, no stability guard, and evanescent content below
 cutoff decays exactly rather than being misread by a direction test that
 cannot work there.

 PLACEMENT
 ---------
 The absorber sheets are at the FIRST and LAST z planes, i.e. directly on the
 PEC domain boundary. That is deliberate. Because the condition is an
 overwrite, whatever lies behind the sheet is driven but can never drive back,
 so a gap between the sheet and the PEC becomes a sealed cavity that fills with
 energy it can never release. Measured in Octave: the guide's own decay is
 identical to within 1e-13 dB whether that gap is 0 or 10 cells -- the trapped
 energy never touches the guide -- but it does stall any global energy
 convergence test, and openEMS uses exactly such a test. With the sheet on the
 PEC face there is no such region at all.

 EXCITATION
 ----------
 The default band (1.30 - 2.90 GHz) sits entirely ABOVE the 1.2 GHz cutoff.
 That matters, and not because of the absorber: a waveguide cannot carry energy
 below its cutoff in the first place, so any below-cutoff content the source
 injects just sits there as a near-static field with no outgoing wave for ANY
 absorber to remove. Right at cutoff the group velocity is zero, so that energy
 never reaches the boundary either. Measured in Octave, moving the source off
 DC was worth ~20 dB across the band and turned "never converges" into a clean
 convergence -- with the absorber completely unchanged. Set DC_HEAVY = True to
 run the stress case and watch the run refuse to settle; that is physics, not a
 defect, and commercial solvers do the same thing.

 (c) 2023-2026 Gadi Lahav <gadi@rfwithcare.com>
"""

# ## Import Libraries
import os, tempfile
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

# ## General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_RectWG_ModalMur')
if not os.path.exists(Sim_Path):
    os.mkdir(Sim_Path)

post_proc_only = False
display_structure = False

# Excitation band. False = off-DC (recommended, see header); True = the DC-heavy
# stress case that deliberately slams a decade below cutoff.
DC_HEAVY = False

# waveguide setup (Octave testbed dimensions, drawing units = mm)
fc_wg = 1.2e9  # TE10 cutoff frequency
unit = 1e-3
wg_a = (C0 / (2.0 * fc_wg)) / unit  # 124.913 mm
wg_b = wg_a / 2.0
wg_L = 1120.0

# mesh (uniform: 24 x 12 cells transverse, dz = 8 mm, as in the Octave testbed)
Nx, Ny = 24, 12
dz = 8.0

if DC_HEAVY:
    f0, fc = 1.55e9, 1.45e9     # reaches ~0.1 GHz: a decade below cutoff
else:
    f0, fc = 2.10e9, 0.80e9     # band 1.30 - 2.90 GHz, entirely above cutoff

# ## Cutoff frequency handed to the absorber
#  The absorber builds beta from the DISCRETE lattice dispersion relation, so
#  the cutoff it should be given is the discrete transverse eigenvalue of this
#  mesh, not the continuum pi/a. For a uniform mesh with PEC side walls the
#  discrete eigenvector is sin(pi*x/a) sampled at the nodes -- exactly the
#  template below -- with eigenvalue kc = (2/dx)*sin(pi*dx/(2a)). It differs
#  from pi/a by only ~0.07% here, but beta = sqrt(k0^2 - kc^2) is extremely
#  sensitive near cutoff, which is precisely where the absorber has to work
#  hardest. Costs one line; pass fc_wg instead if the mode is not analytic.
dx_m = (wg_a * unit) / Nx
a_m = wg_a * unit
kc_num = (2.0 / dx_m) * np.sin(np.pi * dx_m / (2.0 * a_m))
fc_abs = C0 * kc_num / (2.0 * np.pi)
print('TE10 cutoff: continuum {:.4f} GHz, discrete {:.4f} GHz (handed to the absorber)'.format(
      fc_wg / 1e9, fc_abs / 1e9))

# ## Generate the TE10 mode files (local coords in drawing units, columns x,y,Fx,Fy)
#    E:  Ey ~  sin(pi*x/a)          H:  Hx ~ -sin(pi*x/a)   (pair for +z travel)
#  Modal Mur needs only the E file; the H file is here for the waveguide ports.
E_mode_file = os.path.join(Sim_Path, 'RectWG_Er.csv')
H_mode_file = os.path.join(Sim_Path, 'RectWG_Hr.csv')
xm, ym = np.meshgrid(np.linspace(0.0, wg_a, 101), np.linspace(0.0, wg_b, 51))
prof = np.sin(np.pi * xm / wg_a).reshape(-1, 1)
zero = np.zeros_like(prof)
np.savetxt(E_mode_file, np.hstack([xm.reshape(-1, 1), ym.reshape(-1, 1), zero, prof]), delimiter=',')
np.savetxt(H_mode_file, np.hstack([xm.reshape(-1, 1), ym.reshape(-1, 1), -prof, zero]), delimiter=',')

# ## FDTD setup
FDTD = openEMS(NrTS=30000, EndCriteria=1e-5)
FDTD.SetGaussExcite(f0, fc)
FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'PEC', 'PEC', 'PEC'])  # honest PEC guide

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)

mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.AddLine('x', np.linspace(0.0, wg_a, Nx + 1).tolist())
mesh.AddLine('y', np.linspace(0.0, wg_b, Ny + 1).tolist())
mesh.AddLine('z', np.arange(0.0, wg_L + dz / 2, dz).tolist())

Zz = mesh.GetLines('z')

# Absorbers ON the PEC end faces -- the first and last z planes. No pocket.
idxAbs1 = 0
idxAbs2 = len(Zz) - 1

# Ports well clear of both the sheets and each other. A port sitting in a
# sheet's near field measures the sheet, not a travelling wave, and the
# absorber's read plane must not see raw source injection either.
idxPort1 = 12
idxPort2 = len(Zz) - 1 - 12

kc_TE10 = np.pi / (wg_a * unit)  # 1/m, for the dispersive port math

# --- Absorber 1: terminates the low-z PEC face -------------------------------
#  normal_positive=True: the guide lies at HIGHER index than the sheet, so the
#  read plane is one cell in the +z direction.
abs_z = Zz.item(idxAbs1)
modal_mur_1 = FDTD.AddModalAbsorber([0.0, 0.0, abs_z], [wg_a, wg_b, abs_z], 'z',
                                    E_file=E_mode_file,
                                    mode_type='TE',      # -> dispersive Modal Mur
                                    fc=fc_abs,
                                    normal_positive=True)

# --- Port 1: waveguide port with excitation (mode files, TE -> excite_type 0) ---
start = [0.0, 0.0, Zz.item(idxPort1 + 0)]
stop = [wg_a, wg_b, Zz.item(idxPort1 + 1)]
port1 = FDTD.AddWaveGuidePort(1, start, stop, 'z',
                              E_file=E_mode_file, H_file=H_mode_file,
                              kc=kc_TE10, excite=1, excite_type=0)

# --- Port 2: passive waveguide port for the S21 measurement ------------------
# Reversed z order (same convention as the coax test) so the port faces the
# incoming +z wave.
start = [0.0, 0.0, Zz.item(idxPort2 - 0)]
stop = [wg_a, wg_b, Zz.item(idxPort2 - 1)]
port2 = FDTD.AddWaveGuidePort(2, start, stop, 'z',
                              E_file=E_mode_file, H_file=H_mode_file,
                              kc=kc_TE10, excite=0, excite_type=0)

# --- Absorber 2: terminates the high-z PEC face ------------------------------
#  normal_positive=False: the guide lies at LOWER index, read plane is -z.
abs_z = Zz.item(idxAbs2)
modal_mur_2 = FDTD.AddModalAbsorber([0.0, 0.0, abs_z], [wg_a, wg_b, abs_z], 'z',
                                    E_file=E_mode_file,
                                    mode_type='TE',      # -> dispersive Modal Mur
                                    fc=fc_abs,
                                    normal_positive=False)

# Define dump box...
# Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
# Et.AddBox([0.0, 0.0, 0.0], [wg_a, wg_b, wg_L])

if display_structure:
    CSX_file = os.path.join(Sim_Path, 'rectwg_modal_mur.xml')
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=False)

# ## Post-processing
#
# TWO measurements are reported, and they disagree. That disagreement is the
# point, so read both.
#
#  (1) TIME-GATED reflection -- the honest absorber metric.
#      The incident pulse passes port 1, and the reflection off the far
#      absorber comes back a full round trip later. Split the port's raw
#      VOLTAGE record at that boundary and take the spectrum ratio. No
#      impedance, no incident/reflected separation, no E/H colocation: only
#      one probe and the clock. Validated against a bare-PEC far end, which
#      this method reports at -0.16 dB, i.e. the perfect reflector it is.
#
#  (2) |S11| from WaveguidePort.CalcPort -- the convenient metric, and the
#      one that is limited by its own de-embedding rather than by the
#      absorber. On this testbed it reads a monotonic -31 -> -19 dB while the
#      gated measurement shows the dispersive bowl bottoming at -50 dB. The
#      port de-embedding separates incident from reflected using an E/H ratio
#      half a cell apart, and the residual forward leak floors the result.
#
# Believe (1) about the absorber; use (2) for continuity with other tests.
f = np.linspace(0.1e9, 3e9, 401)


def gated_reflection(sim_path, port_file, f_eval, t_gate):
    """|Gamma|(f) from one voltage record, split at t_gate into incident/reflected."""
    d = np.loadtxt(os.path.join(sim_path, port_file), comments='%')
    t, u = d[:, 0], d[:, 1]
    dt_s = t[1] - t[0]
    W = np.exp(-2j * np.pi * np.outer(f_eval, t)) * dt_s
    inc = W @ (u * (t < t_gate))
    ref = W @ (u * (t >= t_gate))
    return np.abs(ref / inc)


port1.CalcPort(Sim_Path, f)
port2.CalcPort(Sim_Path, f)

s11 = port1.uf_ref / port1.uf_inc
s21 = port2.uf_ref / port1.uf_inc
s11_dB = 20.0 * np.log10(np.abs(s11))
s21_dB = 20.0 * np.log10(np.abs(s21))

# Octave reference (WG_ModalMur_BoundaryTermination.m, pocket 0, off-DC source).
# Those numbers come from that script's own de-embedding at +-12 cells, so they
# are an order-of-magnitude orientation, not a bit-for-bit target.
f_tab = np.array([1.3, 1.5, 1.7, 1.9, 2.11, 2.30, 2.50, 2.70, 2.90]) * 1e9
oct_offdc = [-35.4, -46.8, -54.9, -63.2, -66.3, -70.4, -78.6, -75.0, -64.9]
oct_dc = [-37.1, -44.1, -63.4, -51.4, -52.6, -51.6, -49.9, -46.6, -42.3]
ref = oct_dc if DC_HEAVY else oct_offdc

# Earliest possible arrival of the far-end reflection at port 1: the round trip
# at the speed of light (the group velocity is always slower, so nothing can
# beat this). Offset by the incident peak so the gate sits in the quiet stretch
# between the two.
d_round = 2.0 * (Zz.item(idxAbs2) - Zz.item(idxPort1)) * unit
u1 = np.loadtxt(os.path.join(Sim_Path, 'port_ut_1'), comments='%')
t_peak = u1[np.argmax(np.abs(u1[:, 1])), 0]
t_gate = t_peak + 0.95 * d_round / C0
gam_dB = 20.0 * np.log10(gated_reflection(Sim_Path, 'port_ut_1', f, t_gate))
print('\ntime gate at {:.2f} ns (incident peak {:.2f} ns + round trip {:.2f} ns)'.format(
      t_gate * 1e9, t_peak * 1e9, d_round / C0 * 1e9))

print('\n===== Dispersive Modal Mur, {} source ====='.format(
      'DC-heavy' if DC_HEAVY else 'off-DC'))
print('  f[GHz]   octave-ref   GATED |Gam|   port |S11|   |S21| dB')
for i, fi in enumerate(f_tab):
    j = np.abs(f - fi).argmin()
    print('   {:.2f}     {:8.2f}     {:9.2f}     {:8.2f}   {:8.2f}'.format(
          fi / 1e9, ref[i], gam_dB[j], s11_dB[j], s21_dB[j]))

# Restrict the verdict to where the gating is trustworthy. Right above cutoff
# the group velocity collapses, so the reflection has not fully returned within
# the record and BOTH gated numbers are optimistic there -- the bare-PEC control
# reads -2.7 dB at 1.3 GHz instead of 0 dB, which is how that shows up.
band = (f >= 1.4e9) & (f <= 2.9e9)
print('\nworst GATED |Gamma| over 1.4-2.9 GHz: {:.2f} dB   <-- the absorber'.format(np.max(gam_dB[band])))
print('worst port |S11|  over 1.4-2.9 GHz: {:.2f} dB   <-- the de-embedding floor'.format(np.max(s11_dB[band])))

figure()
plot(f / 1e9, gam_dB, 'r-', linewidth=2, label='$|\\Gamma|$ time-gated (the absorber)')
plot(f / 1e9, s11_dB, 'k--', linewidth=1.5, label='$S_{11}$ from port (de-embed limited)')
plot(f / 1e9, s21_dB, 'b-', linewidth=2, label='$S_{21}$ (through, expect $\\approx$ 0 dB)')
axvline(fc_abs / 1e9, color='g', linestyle=':', linewidth=1.5, label='$f_{cut}$')
grid()
legend()
ylabel('S-Parameter (dB)')
xlabel('Frequency (GHz)')
title('Rect WG TE10 — Dispersive Modal Mur on both PEC faces')
show()
