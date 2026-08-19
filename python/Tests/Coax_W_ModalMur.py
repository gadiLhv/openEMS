"""
 Coaxial line (TEM) terminated by modal absorbers, with no open ends.

 WHY THIS EXISTS ALONGSIDE Coax_W_Modal_Absorb.py
 ------------------------------------------------
 That test puts the coax inside an airbox that extends 1 mm past each end, so
 the line has two UNTERMINATED OPEN ENDS radiating into air, with the absorber
 sheets sitting 5 and 2 cells short of them. Whatever the sheets do, the open
 ends reflect and fringe, and that fringing dominates the measurement -- it is
 clearly visible in a field dump as a bright spot hugging the conductors at the
 end of the line. No absorber can fix a geometry that is radiating past it.

 Here the coax spans the ENTIRE z domain. There is no air beyond it, no open
 end, and nothing for the sheet to miss.

 THE TWO METHODS (MODE switch below)
 -----------------------------------
 MODE = 'MUR'    Dispersive Modal Mur. The z faces are PEC and the sheets sit
                 DIRECTLY ON THEM, on the first and last z planes. The method
                 is a one-way overwrite, so it must terminate a PEC face: any
                 gap left behind it becomes a sealed cavity that fills with
                 energy it can never release. A coax TEM mode is the kc -> 0
                 limit of a TM mode, so it is declared as mode_type='TM' with
                 fc = 0; the tap generator then reduces to the pure per-cell
                 delay exp(-j*w*dz/v).

 MODE = 'SCALAR' The scalar Zw splitter. It measures direction from E = +-Zw*H
                 and therefore needs a live E field at its own plane, so it
                 CANNOT sit on a PEC face. The z faces are MUR and the sheets
                 sit SCALAR_INSET cells inside. That inset is not cosmetic: at
                 one cell the sheet and the MUR stencil overlap and the run
                 DIVERGES (energy -> inf at ~7.5k steps). Two cells is already
                 stable; four is the default, for margin and because it
                 measures best.

 WHICH ONE TO USE ON A COAX -- READ THE CONTROLS FIRST
 -----------------------------------------------------
 The two methods cannot be given the same surroundings. Modal Mur is a one-way
 overwrite and REQUIRES a PEC face; the scalar splitter needs a live E at its
 plane and CANNOT have one. So any head-to-head compares arrangements, not
 sheets, and the surroundings have to be measured too. On this bench, same
 geometry, same ruler, worst |Gamma| in band:

     PEC face + Modal Mur sheet on it .............. -13.7 dB
     MUR face, NO sheet at all  (control) .......... -22.3 dB
     MUR face + scalar sheet 4 cells in ............ -31.9 dB
     PEC face + scalar sheet 4 cells in ............ DIVERGES

 Read those four lines together, because three of them are uncomfortable:

  - Plain first-order Mur already beats the Modal Mur sheet by 8.6 dB. That is
    not a scandal, it is what Mur IS: exact for a normally incident wave at the
    velocity it was tuned to, and a coax TEM wave arriving at a flat end face is
    exactly that. There is very little here for a modal absorber to add.

  - The scalar arrangement's -31.9 dB is NOT the sheet's own figure. The sheet
    is worth about 9.6 dB on top of the Mur boundary behind it; the Mur is
    worth 22.3.

  - Put PEC behind the scalar sheet -- taking the Mur away so the sheet has to
    do all the work -- and it does not merely get worse, it DIVERGES. With
    nu = 0.026 the matched correction gain is 0.051, so the sheet corrects far
    too slowly to control the shorted stub it now sits in front of.

 The honest conclusion: on a coax, use the scalar sheet in front of a Mur
 boundary; do not read its number as the sheet's own performance, and do not
 conclude from it that the scalar splitter is the stronger absorber. The place
 to compare the two methods is a guide with a REAL CUTOFF, where Mur has no
 single correct velocity to be tuned to and the direction test E = +-Zw*H
 breaks down -- which is what forced the one-way reformulation in the first
 place. See RectWG_W_ModalMur.py.

 MEASUREMENT
 -----------
 The line is short, so a time gate does not fit in the record. Instead this
 uses the same E-only multi-plane modal fit as RectWG_W_ModalMur.py: the modal
 voltage in a source-free stretch is A*exp(-g*z) + B*exp(+g*z), and B/A at the
 fit plane is exactly the reflection looking toward +z. No current probe, so
 the half-cell E/H inconsistency that puts the port energy balance at 1.17 on
 this structure never enters.

 One difference from the waveguide: gamma is NOT known here. Staircasing the
 round conductors on a Cartesian mesh moves the effective permittivity off
 2.5, and a beta that is even slightly wrong contaminates the A/B split. So
 beta is FITTED per frequency -- chosen to minimise the two-exponential
 residual -- and the effective index it lands on is reported. That number is
 worth reading on its own: it says how much the mesh is distorting the line.

 (c) 2023-2026 Gadi Lahav <gadi@rfwithcare.com>
"""

import os, tempfile, shutil
from pylab import *

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS import utilities

# ## Which termination to test
MODE = 'MUR'          # 'MUR' (Modal Mur on PEC) or 'SCALAR' (Zw splitter + Mur face)

# ## Watching it with your own eyes
#
# VIZ trades a little accuracy for a WATCHABLE animation, and it is not a
# cosmetic switch -- with the analysis settings the picture is useless. The
# analysis pulse is 269 mm long against an 80 mm line, i.e. 3.4x the whole
# structure, so the entire coax lights up and dies down together and an
# incident wave cannot be told from a reflected one by eye.
#
# VIZ shortens the pulse and doubles the line to 160 mm, so the pulse is a
# compact blob that runs down the line, reaches the sheet on the PEC end face,
# and either disappears into it or comes back as a second blob. If it comes
# back you can SEE it, and that is the entire point.
#
# Two details that are not free choices:
#
#  - The short pulse is made by RAISING f0 with fc, to 6 +- 4 GHz. Widening fc
#    alone (2.5 +- 6 GHz) puts enormous content at DC, and DC trapped between
#    two PEC end faces never leaves: the Modal Mur delay has |D| = 1 at DC, so
#    it does not remove it either. Measured: the run stalls at -36.6 dB and sits
#    there for 150k timesteps. Same trap as the rectangular guide.
#
#  - The mesh is deliberately NOT refined to match the shorter pulse (see
#    f_mesh below); that would cost 14x the cells for a picture.
#
# VIZ numbers are indicative only. For the validated S-parameters set VIZ=False.
VIZ = True

Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_Coax_ModalMur_' + MODE)
if not os.path.exists(Sim_Path):
    os.mkdir(Sim_Path)
shutil.copy("Coax_Er.csv", Sim_Path)
shutil.copy("Coax_Hr.csv", Sim_Path)

post_proc_only = False
display_structure = False

# ## Geometry (drawing units = mm)
coax_D = 2.0                    # inner diameter of the shield
coax_shield_thick = 0.15
coax_wire_D = 0.5
coax_L = 160.0 if VIZ else 80.0  # long enough to hold ports and the fit ladder
teflon_epsR = 2.5
Airbox_Add = 1.0                # transverse only -- NONE in z, that is the point
unit = 1e-3

Zw_TEM = 238.26517157           # coax modal impedance (as in Coax_W_WG_Ports)
# EFFECTIVE INDEX -- measure it, do not assume it.
#
# sqrt(eps_r) = 1.5811 is what an ideal coax would do. This one does not: the
# round conductors are staircased onto a Cartesian mesh, and the modal fit at
# the bottom of this script reports the index the line ACTUALLY has. On this
# mesh that is 1.826, 15% slow. The Modal Mur delay filter is exp(-j*w*dz/v),
# so handing it the ideal speed builds the wrong delay for every tap.
#
# The script is self-calibrating in one iteration: run it, read the n_eff it
# prints, put that number here, run again.
N_EFF = 1.826            # measured on this mesh; 1.5811 = sqrt(eps_r), the ideal
v_ph = C0 / N_EFF

r_out = coax_D * 0.5 + coax_shield_thick        # 1.15 mm, matches Coax_Er.csv extent

f0, fc_exc = (6.0e9, 4.0e9) if VIZ else (2.5e9, 1.0e9)
# Mesh resolution is pinned to the ANALYSIS band on purpose. Deriving it from
# the excitation, as this script used to, would refine the mesh 2.4x in every
# direction under VIZ -- 14x the cells to make a picture.
f_mesh = 3.5e9

# ## FDTD setup
# VIZ: cap the run at ~5 one-way transits -- long enough to watch the pulse
# arrive and to see whether anything comes back, without waiting for a full
# energy decay. OverSampling=1 keeps the frame count sane (4 would be ~1300).
FDTD = openEMS(NrTS=(45000 if VIZ else 300000), EndCriteria=1e-5,
               OverSampling=(1 if VIZ else 4))
FDTD.SetGaussExcite(f0, fc_exc)
if MODE == 'MUR':
    # PEC z faces: the Modal Mur sheets terminate them directly.
    FDTD.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'PEC', 'PEC'])
else:
    # The scalar splitter needs a live E at its plane, so no PEC there.
    FDTD.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh_res = ((C0 / f_mesh) / unit) / 100.0

SimBox = np.array([-(r_out + Airbox_Add), (r_out + Airbox_Add),
                   -(r_out + Airbox_Add), (r_out + Airbox_Add),
                   0.0, coax_L])            # <-- z is EXACTLY the coax, no airbox
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

# center wire -- spans the whole domain
line = CSX.AddMetal('Wire_Inner')
line.AddCylinder(priority=10, start=[0, 0, 0.0], stop=[0, 0, coax_L], radius=coax_wire_D * 0.5)
mesh.AddLine('x', np.linspace(-coax_wire_D * 0.5, coax_wire_D * 0.5, 6).tolist())
mesh.AddLine('y', np.linspace(-coax_wire_D * 0.5, coax_wire_D * 0.5, 6).tolist())

# outer shield -- spans the whole domain
shield = CSX.AddMetal('Shield_Outer')
shield.AddCylindricalShell(priority=10, start=[0, 0, 0.0], stop=[0, 0, coax_L],
                           radius=(coax_D + coax_shield_thick) * 0.5,
                           shell_width=coax_shield_thick)
rad, hth = (coax_D + coax_shield_thick) * 0.5, coax_shield_thick * 0.5
mesh.AddLine('x', np.linspace(-(rad + hth), -(rad - hth), 4).tolist() +
                  np.linspace((rad - hth), (rad + hth), 4).tolist())
mesh.AddLine('y', np.linspace(-(rad + hth), -(rad - hth), 4).tolist() +
                  np.linspace((rad - hth), (rad + hth), 4).tolist())

# teflon fill -- spans the whole domain
teflon = CSX.AddMaterial('PTFE', epsilon=teflon_epsR)
teflon.AddCylindricalShell(priority=8, start=[0, 0, 0.0], stop=[0, 0, coax_L],
                           radius=(coax_wire_D + coax_D) * 0.25,
                           shell_width=(coax_D - coax_wire_D) * 0.5)
rad, hth = (coax_wire_D + coax_D) * 0.25, (coax_D - coax_wire_D) * 0.25
mesh.AddLine('x', np.linspace(-(rad + hth), -(rad - hth), 12).tolist() +
                  np.linspace((rad - hth), (rad + hth), 12).tolist())
mesh.AddLine('y', np.linspace(-(rad + hth), -(rad - hth), 12).tolist() +
                  np.linspace((rad - hth), (rad + hth), 12).tolist())

mesh.SmoothMeshLines('all', mesh_res, 1.25)

Zz = np.asarray(mesh.GetLines('z'))
nz = len(Zz)

# Sheet planes. MUR: on the PEC faces. SCALAR: inset, since it needs a live E
# at its own plane and must not overlap the MUR stencil (see header -- one cell
# diverges, two is stable, four measured best).
SCALAR_INSET = 4
idxAbs1 = 0 if MODE == 'MUR' else SCALAR_INSET
idxAbs2 = (nz - 1) if MODE == 'MUR' else (nz - 1 - SCALAR_INSET)
idxPort1 = idxAbs1 + 8
idxPort2 = idxAbs2 - 8

box_lo = [-r_out, -r_out]
box_hi = [r_out, r_out]


def sheet(idx, normal_positive):
    z = Zz.item(int(idx))
    kw = dict(E_file="Coax_Er.csv", normal_positive=normal_positive)
    if MODE == 'MUR':
        # TEM = the kc -> 0 limit of TM. fc = 0 makes the tap generator collapse
        # to the pure per-cell delay; phase_velocity carries the dielectric.
        kw.update(mode_type='TM', fc=0.0, phase_velocity=v_ph)
    else:
        kw.update(mode_type='TEM', H_file="Coax_Hr.csv", Zw=Zw_TEM)
    return FDTD.AddModalAbsorber(box_lo + [z], box_hi + [z], 'z', **kw)


abs1 = sheet(idxAbs1, True)

port1 = FDTD.AddWaveGuidePort(1, box_lo + [Zz.item(idxPort1)], box_hi + [Zz.item(idxPort1 + 1)],
                              'z', E_file="Coax_Er.csv", H_file="Coax_Hr.csv",
                              kc=0.0, excite=1, excite_type=0)
port2 = FDTD.AddWaveGuidePort(2, box_lo + [Zz.item(idxPort2)], box_hi + [Zz.item(idxPort2 - 1)],
                              'z', E_file="Coax_Er.csv", H_file="Coax_Hr.csv",
                              kc=0.0, excite=0, excite_type=0)

abs2 = sheet(idxAbs2, False)

# ## E-only modal fit ladder, in the source-free stretch between the ports
FIT_OFFSETS = np.array([0, 1, 2, 3, 5, 8, 13, 21, 34])
idxFit0 = idxPort1 + 4
fit_idx = idxFit0 + FIT_OFFSETS
assert fit_idx[-1] < idxPort2 - 1, 'fit ladder collides with port 2; lengthen coax_L'
fit_files = []
for n, idx in enumerate(fit_idx):
    zf = Zz.item(int(idx))
    fn = 'mfit_ut_{:02d}'.format(n)
    p = CSX.AddProbe(fn, p_type=10, mode_file_name="Coax_Er.csv")
    p.AddBox(box_lo + [zf], box_hi + [zf])
    fit_files.append(fn)
z_fit_m = (Zz[fit_idx] - Zz.item(int(idxFit0))) * unit

### Field export -------------------------------------------------------------
#  A y = 0 CUT through the coax axis, not the whole box. The full volume is
#  ~2.9 MB per frame and openEMS writes hundreds of frames; the axial cut is
#  ~60 kB and is the view you actually want -- it shows the wave running down
#  the line and what happens when it reaches the sheet on the PEC end face.
#
#  Open Sim_Path/Et_*.vtr as a time series in ParaView.
DUMP_FIELDS = VIZ
if DUMP_FIELDS:
    Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
    Et.AddBox([SimBox[0], 0.0, SimBox[4]], [SimBox[1], 0.0, SimBox[5]])

if display_structure:
    CSX_file = os.path.join(Sim_Path, 'coax_modal_mur.xml')
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=False)

# ## Post-processing
f = np.linspace(2.0e9, 10.0e9, 261) if VIZ else np.linspace(1.2e9, 3.8e9, 261)


def modal_split(z_off_m, Uf_row, gam):
    """LSQ split into forward/backward at the reference plane, one frequency.
    Returns A, B, condition number and relative residual."""
    M = np.stack([np.exp(-gam * z_off_m), np.exp(+gam * z_off_m)], axis=1)
    scale = np.linalg.norm(M, axis=0)
    Mn = M / scale
    sol = np.linalg.lstsq(Mn, Uf_row, rcond=None)[0]
    nrm = np.linalg.norm(Uf_row)
    res = np.linalg.norm(Mn @ sol - Uf_row) / nrm if nrm > 0 else np.inf
    A, B = sol / scale
    return A, B, np.linalg.cond(Mn), res


# Read the ladder
Uf = np.empty((len(f), len(fit_files)), dtype=complex)
for n, fn in enumerate(fit_files):
    d = np.loadtxt(os.path.join(Sim_Path, fn), comments='%')
    Uf[:, n] = utilities.DFT_time2freq(d[:, 0], d[:, 1], f)

# Fit the effective index per frequency: beta is not known a priori on a
# staircased coax, and a wrong beta leaks forward into backward.
n_scan = np.linspace(1.20, 1.90, 351)      # sqrt(2.5) = 1.581 is the ideal
Gam = np.zeros(len(f), dtype=complex)
n_eff = np.zeros(len(f))
cnd = np.zeros(len(f))
res = np.zeros(len(f))
for i, fi in enumerate(f):
    k0 = 2 * np.pi * fi / C0
    best = None
    for ne in n_scan:
        A, B, c, r = modal_split(z_fit_m, Uf[i], 1j * k0 * ne)
        if (best is None) or (r < best[0]):
            best = (r, ne, A, B, c)
    res[i], n_eff[i], A, B, cnd[i] = best
    Gam[i] = B / A
Gam_dB = 20 * np.log10(np.abs(Gam))

# Conventional port S-parameters, for continuity with the other tests
Zw = np.array([Zw_TEM])
port1.CalcPort(Sim_Path, f, ref_impedance=Zw)
port2.CalcPort(Sim_Path, f, ref_impedance=Zw)
s11 = np.asarray(port1.uf_ref / port1.uf_inc).ravel()
s21 = np.asarray(port2.uf_ref / port1.uf_inc).ravel()
s11_dB = 20 * np.log10(np.abs(s11))
s21_dB = 20 * np.log10(np.abs(s21))

print('\n===== Coax TEM, {} termination, no open ends ====='.format(MODE))
print('  ideal effective index sqrt(eps_r) = {:.4f}'.format(np.sqrt(teflon_epsR)))
print('  f[GHz]   MODAL-FIT |Gam|   n_eff   resid    port|S11|   |S21| dB')
for fi in ([2.5, 4.0, 6.0, 8.0, 9.5] if VIZ else [1.5, 2.0, 2.5, 3.0, 3.5]):
    j = np.abs(f - fi * 1e9).argmin()
    print('   {:.2f}       {:9.2f}     {:.4f}  {:.4f}    {:8.2f}   {:8.2f}'.format(
          fi, Gam_dB[j], n_eff[j], res[j], s11_dB[j], s21_dB[j]))
# Print WHERE the worst point is, not just its value. A residual threshold is
# not enough to police this: under VIZ the fit degrades over the top ~0.3 GHz
# of the band and reads -0.1 dB, but the very worst point has a residual of
# 0.002 -- it is fitting something coherent, just not the mode. Naming the
# frequency lets the reader see at a glance that it is a band edge rather than
# a reflection the absorber is failing to swallow.
jw = np.argmax(Gam_dB)
print('worst MODAL-FIT |Gamma|: {:.2f} dB at {:.2f} GHz (resid {:.4f})   <-- the absorber'.format(
      Gam_dB[jw], f[jw] / 1e9, res[jw]))
core = (f > f[0] + 0.1 * (f[-1] - f[0])) & (f < f[-1] - 0.1 * (f[-1] - f[0]))
print('   ... excluding the outer 10% of the band: {:.2f} dB'.format(np.max(Gam_dB[core])))
print('worst port |S11|        in band: {:.2f} dB   <-- de-embed limited'.format(np.max(s11_dB)))
bal = np.abs(s11) ** 2 + np.abs(s21) ** 2
print('energy balance |S11|^2+|S21|^2: min {:.3f}  max {:.3f}{}'.format(
      bal.min(), bal.max(),
      '   <-- off 1: port de-embed, not the absorber' if (bal.max() > 1.02 or bal.min() < 0.98) else ''))
print('modal-fit residual: max {:.2e}   (two-exponential model)'.format(np.max(res)))

figure()
plot(f / 1e9, Gam_dB, 'm-', linewidth=2, label=r'$|\Gamma|$ modal fit (E-only)')
plot(f / 1e9, s11_dB, 'k--', linewidth=1.5, label='$S_{11}$ from port')
plot(f / 1e9, s21_dB, 'b-', linewidth=2, label='$S_{21}$ (through)')
grid(); legend(); ylim(-70, 5)
ylabel('S-Parameter (dB)'); xlabel('Frequency (GHz)')
title('Coax TEM — {} termination, coax spans the whole domain'.format(MODE))

figure()
plot(f / 1e9, n_eff, 'k-', linewidth=2, label='fitted $n_{eff}$')
axhline(np.sqrt(teflon_epsR), color='r', linestyle='--', label=r'$\sqrt{\epsilon_r}$')
grid(); legend(); ylabel('effective index'); xlabel('Frequency (GHz)')
title('How far the staircased mesh moves the line off its ideal index')
show()
