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

 POST-PROCESSING
 ---------------
 The reported S-parameters are stabilised near and below cutoff by spending the
 waveguide's closed-form broadband response three ways: the grid's own discrete
 dispersion instead of the continuum one, an E-only multi-plane modal fit that
 needs neither a current probe nor a wave impedance nor a time gate, and a
 renormalisation of the port S-parameters to a real reference impedance so they
 stay bounded where the modal impedance is singular or reactive. See the long
 comment block above the post-processing section for why each is needed.

 (c) 2023-2026 Gadi Lahav <gadi@rfwithcare.com>
"""

# ## Import Libraries
import os, tempfile
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS import utilities

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

# --- Absorber sheets ---------------------------------------------------------
#  Same interface as Coax_W_ModalMur.py and CPW_W_ModalAbsorb.py.
#
#  Unlike those two, this guide has a REAL cutoff, so fc is a genuine mode
#  constant rather than the free-space-referencing artefact a TEM line reports.
#  fc_abs is the DISCRETE TE10 cutoff (see above), and phase_velocity is left at
#  its C0 default because the guide is air-filled: fc and v are then both
#  referenced to free space, which is the consistency modal_mur_taps.cpp:62
#  requires when it rebuilds kc = 2*pi*fc/v.
#
#  normal_positive=True: the guide lies at HIGHER index than the sheet, so the
#  read plane is MODAL_MUR_READ_CELLS cells in the +z direction.
def sheet(idx, normal_positive):
    z = Zz.item(int(idx))
    kw = dict(E_file=E_mode_file, normal_positive=normal_positive)
    kw.update(mode_type='TE', fc=fc_abs)   # -> dispersive Modal Mur
    return FDTD.AddModalAbsorber([0.0, 0.0, z], [wg_a, wg_b, z], 'z', **kw)


modal_mur_1 = sheet(idxAbs1, True)

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
modal_mur_2 = sheet(idxAbs2, False)

# ## Modal-fit probe ladder ---------------------------------------------------
#  Nine E-only mode-matched planes in the source-free stretch between the two
#  ports. They cost nine scalars per timestep, no field dump, and they are what
#  makes the near/below-cutoff numbers in the post-processing possible.
#
#  Offsets form a Fibonacci ladder because a two-plane forward/backward split is
#  blind wherever beta*d = n*pi (the two basis exponentials coincide there).
#  With spacings 1,2,3,5,8,13,21,34 cells no in-band frequency is near a null of
#  ALL pairs at once, so the joint least-squares stays well conditioned from the
#  cutoff knee up to 3 GHz. The long end (34 cells = 272 mm) is what keeps the
#  fit alive just above cutoff, where beta is small; the short end is what keeps
#  it alive at 2.9 GHz, where lambda_g/2 is only ~7 cells.
FIT_OFFSETS = np.array([0, 1, 2, 3, 5, 8, 13, 21, 34])
idxFit0 = idxPort1 + 6
fit_idx = idxFit0 + FIT_OFFSETS
fit_files = []
for n, idx in enumerate(fit_idx):
    zf = Zz.item(int(idx))
    fname = 'mfit_ut_{:02d}'.format(n)
    p_fit = CSX.AddProbe(fname, p_type=10, mode_file_name=E_mode_file)
    p_fit.AddBox([0.0, 0.0, zf], [wg_a, wg_b, zf])
    fit_files.append(fname)
z_fit_m = (np.asarray(Zz)[fit_idx] - Zz.item(int(idxFit0))) * unit

### Field export -- disabled. Uncomment to dump E(t) over the whole guide.
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
# ============================================================================
# STABILISING THE S-PARAMETERS NEAR AND BELOW CUTOFF
# ============================================================================
# A rectangular waveguide is one of the very few structures whose broadband
# response is known in closed form, and that knowledge is worth spending,
# because the standard port de-embedding falls apart exactly where this
# absorber has to prove itself.
#
# WHY IT FALLS APART. WaveguidePort splits total into incident/reflected with
#     U_inc = (U + I*ZL)/2 ,  ZL = k*Z0/beta .
# Near cutoff beta -> 0 so ZL -> infinity: U is swamped by I*ZL and |S11| is
# then a measurement of the current probe's error, amplified by ZL. Below
# cutoff ZL is purely reactive, so the "incident wave" it defines is not a wave
# at all -- |S11| is unbounded there and the energy balance is meaningless.
# Neither is a defect of the absorber; both are properties of the reference.
#
# THREE FIXES, all of them just the analytic guide model spent in the right
# place, and all of them local to this script:
#
#  (1) USE THE GRID'S OWN DISPERSION, NOT THE CONTINUUM'S. beta comes from the
#      discrete lattice relation -- the same three-branch formula the C++ tap
#      generator uses -- so the ruler and the absorber model agree exactly.
#      Continuum and discrete kc differ by 0.07% here, but beta = sqrt(k0^2 -
#      kc^2) turns that into a large relative error near cutoff, which is
#      precisely the band in question.
#
#  (2) MEASURE WITH E ONLY, ON MANY PLANES (the modal fit). The modal voltage
#      in a source-free single-mode stretch is exactly
#          U(z) = A*exp(-gamma*z) + B*exp(+gamma*z)
#      with gamma known. Sample U on a ladder of planes and least-squares for
#      (A, B). Gamma_fit = B/A is the reflection looking toward +z -- and it is
#      exact, not an approximation: everything travelling +z at that plane is
#      what hits the far absorber, everything travelling -z is what it returned,
#      so multiple bounces off the near end never contaminate it.
#      No current probe, so no half-cell E/H inconsistency (the thing that puts
#      the port's energy balance above 1). No ZL, so no cutoff singularity. And
#      no time gate, so unlike the gated number it survives just above cutoff
#      where the group velocity collapses and the round trip no longer fits in
#      the record.
#      Its own limit is honest physics: at exactly f = fc the two basis
#      exponentials become the same constant and forward cannot be told from
#      backward by ANY method. That shows up as the condition number below, and
#      it is reported rather than hidden.
#
#  (3) RENORMALISE S11 TO A REAL REFERENCE, AND READ IT AS PASSIVITY. Referred
#      to a constant real Zr instead of the singular modal ZL,
#          S11 = (Zin - Zr)/(Zin + Zr)
#      is finite for every Zin, and |S11| > 1 if and only if Re(Zin) < 0. That
#      is an exact, reference-independent statement: it says the terminated
#      guide is DELIVERING energy to the port. It is defined at every frequency
#      including DC and straight through the cutoff singularity, where no
#      reflection coefficient is. It is not an absorber figure of merit -- most
#      of what it shows is the ZL(f)/Zr mismatch, which is the guide, not the
#      sheet -- it is the monitor that says the sheet is not blowing up.
#      Above cutoff ZL -> Z0-ish and it converges back onto the modal S11.
#
# And one free consistency check: between the two ports there is nothing but
# guide, so the analytic S21 is simply exp(-gamma*d). Plotting the measurement
# against it validates gamma, the probes and the DFT in one line -- including
# below cutoff, where it is the evanescent decay exp(-alpha*d).
#
# The gated |Gamma| from before is kept: above cutoff it is an independent
# ruler that shares no machinery with the modal fit, and the two agreeing is
# worth more than either alone.
#
# WHAT "REFLECTION" EVEN MEANS BELOW CUTOFF. Not "all the power comes back".
# In a semi-infinite guide below cutoff the field is purely the DECAYING
# branch; the growing branch is forbidden by the radiation condition. The Modal
# Mur imposes exactly that, a_sheet = a_inside*exp(-alpha*dz), so the honest
# sub-cutoff score is how much growing branch survives -- |B/A| from the modal
# fit, which is a small number when the sheet is right. It is a real metric,
# not a placeholder, and it is available at frequencies where the concepts of
# incident and reflected power do not exist.
#
# CAVEAT ON BELOW-CUTOFF NUMBERS: they are only as real as the excitation
# content down there. With DC_HEAVY = False the Gaussian has essentially
# nothing below ~1.2 GHz, so those columns are noise divided by noise. They are
# masked out on that basis (SIG_FLOOR), together with anything the fit itself
# reports as untrustworthy (COND_MAX, RES_MAX). Set DC_HEAVY = True to actually
# drive the sub-cutoff band.
# ============================================================================

f = np.linspace(0.05e9, 3.0e9, 591)
COND_MAX = 50.0     # modal fit unusable above this (only happens at the knee)
SIG_FLOOR = 1e-3    # incident spectrum below this fraction of peak -> no data
RES_MAX = 0.05      # two-exponential model must explain the ladder to 5%

Zz_m = np.asarray(Zz) * unit


def wg_gamma_discrete(freq, kc, dz_m, dt_s, v=C0):
    """gamma(f) = alpha + j*beta [1/m] from the DISCRETE lattice dispersion.

    Same three branches as the C++ tap generator (modal_mur_taps.cpp), so the
    post-processing and the absorber describe the same guide:
      argd < 0        evanescent, below the modal cutoff
      sqrt(argd)*dz/2 <= 1  propagating
      sqrt(argd)*dz/2 >  1  above the GRID's own cutoff (out of band here)
    """
    w = 2.0 * np.pi * np.asarray(freq, dtype=float)
    keff = (2.0 / (v * dt_s)) * np.sin(w * dt_s / 2.0)
    argd = keff ** 2 - kc ** 2
    sarg = 0.5 * dz_m * np.sqrt(np.abs(argd))
    gam = np.zeros(len(w), dtype=complex)
    ev = argd < 0.0
    pr = (~ev) & (sarg <= 1.0)
    gc = (~ev) & (sarg > 1.0)
    gam[ev] = 2.0 * np.arcsinh(sarg[ev]) / dz_m
    gam[pr] = 2j * np.arcsin(np.clip(sarg[pr], 0.0, 1.0)) / dz_m
    gam[gc] = (2.0 * np.arccosh(sarg[gc]) + 1j * np.pi) / dz_m
    return gam


def modal_split(z_off_m, Uf, gam):
    """Least-squares split of modal voltage samples into forward/backward.

    z_off_m : (Np,)      plane offsets from the reference plane, metres
    Uf      : (Nf, Np)   complex modal voltage spectra, one column per plane
    gam     : (Nf,)      propagation constant

    Returns A, B (amplitudes AT the reference plane), the condition number of
    the equilibrated 2-column basis, and the relative fit residual. The columns
    are normalised before the solve so that the condition number reports
    separability -- can forward be told from backward -- rather than the raw
    dynamic range, which is huge and harmless deep below cutoff.

    The residual is the quality control: two exponentials with the KNOWN gamma
    must explain the sampled field to the noise floor. A residual that is not
    small means the ladder is seeing something the single-mode model does not
    contain (a higher-order mode, source leakage, or a wrong gamma), and the
    A/B split it produced should not be believed.
    """
    Nf = len(gam)
    A = np.zeros(Nf, dtype=complex)
    B = np.zeros(Nf, dtype=complex)
    cnd = np.zeros(Nf)
    res = np.zeros(Nf)
    for i in range(Nf):
        M = np.stack([np.exp(-gam[i] * z_off_m), np.exp(+gam[i] * z_off_m)], axis=1)
        scale = np.linalg.norm(M, axis=0)
        Mn = M / scale
        cnd[i] = np.linalg.cond(Mn)
        sol = np.linalg.lstsq(Mn, Uf[i], rcond=None)[0]
        nrm = np.linalg.norm(Uf[i])
        res[i] = np.linalg.norm(Mn @ sol - Uf[i]) / nrm if nrm > 0 else np.inf
        A[i], B[i] = sol / scale
    return A, B, cnd, res


def gated_reflection(sim_path, port_file, f_eval, t_gate):
    """|Gamma|(f) from one voltage record, split at t_gate into incident/reflected."""
    d = np.loadtxt(os.path.join(sim_path, port_file), comments='%')
    t, u = d[:, 0], d[:, 1]
    dt_s = t[1] - t[0]
    W = np.exp(-2j * np.pi * np.outer(f_eval, t)) * dt_s
    inc = W @ (u * (t < t_gate))
    ref = W @ (u * (t >= t_gate))
    return np.abs(ref / inc)


# --- the guide model, from the grid the solver actually ran ------------------
#  The timestep must come from the EXCITATION dump 'et', which openEMS writes
#  every single step. Probe files are decimated (3:1 here), so their sample
#  interval is NOT dt, and feeding that to the lattice dispersion mis-sizes
#  sin(w*dt/2) by ~1.4% at 2.9 GHz -- small, but the ladder accumulates it over
#  15 rad of phase and the fit residual jumps from 1e-3 to 7e-2.
u1 = np.loadtxt(os.path.join(Sim_Path, 'port_ut_1'), comments='%')
et = np.loadtxt(os.path.join(Sim_Path, 'et'), comments='%')
dt_s = et[1, 0] - et[0, 0]
gam = wg_gamma_discrete(f, kc_num, dz * unit, dt_s)
alpha, beta = gam.real, gam.imag
print('\ntimestep {:.4f} ps, discrete cutoff {:.4f} GHz'.format(dt_s * 1e12, fc_abs / 1e9))

# --- (2) the E-only modal fit -----------------------------------------------
Uf_fit = np.empty((len(f), len(fit_files)), dtype=complex)
for n, fname in enumerate(fit_files):
    d = np.loadtxt(os.path.join(Sim_Path, fname), comments='%')
    Uf_fit[:, n] = utilities.DFT_time2freq(d[:, 0], d[:, 1], f)
A_fit, B_fit, cnd, res = modal_split(z_fit_m, Uf_fit, gam)

# |Gamma| AT THE FIT PLANE, and deliberately not referred anywhere else.
#
# Above cutoff that costs nothing: the guide is lossless, exp(2*gamma*d) has
# unit magnitude, so |Gamma| is the same at every plane including the sheet.
#
# Below cutoff it is the only normalisation that means anything. In a
# semi-infinite guide below cutoff the solution is purely the DECAYING branch,
# B = 0, and that is exactly what the Modal Mur imposes: a_sheet = a_inside *
# exp(-alpha*dz). So |B/A| is a real sub-cutoff absorber metric -- the fraction
# of the field that is the growing branch the boundary should have forbidden --
# and it should be small, not unity. Referring it forward to the sheet instead
# divides by exp(-2*alpha*d), which at 1 GHz here is 1e-12: that ratio is noise
# over noise and reads +180 dB no matter what the sheet does.
Gam_fit = B_fit / A_fit
d_fit_abs = Zz_m[idxAbs2] - Zz_m[idxFit0]
decay_1way_dB = 20.0 * np.log10(np.exp(-gam.real * d_fit_abs))

# --- ports: modal reference (as before) and real reference (stabilised) ------
port1.CalcPort(Sim_Path, f)
port2.CalcPort(Sim_Path, f)

s11 = np.asarray(port1.uf_ref / port1.uf_inc).ravel()
s21 = np.asarray(port2.uf_ref / port1.uf_inc).ravel()
s11_dB = 20.0 * np.log10(np.abs(s11))
s21_dB = 20.0 * np.log10(np.abs(s21))

# (3) same de-embedding, real reference impedance. Zr = Z0 is the guide's own
# high-frequency wave impedance, so the two references merge well above cutoff.
Zr = Z0
U1 = np.asarray(port1.uf_tot).ravel()
I1 = np.asarray(port1.if_tot).ravel()
U2 = np.asarray(port2.uf_tot).ravel()
I2 = np.asarray(port2.if_tot).ravel()
inc1_r = 0.5 * (U1 + Zr * I1)
ref1_r = 0.5 * (U1 - Zr * I1)
trn2_r = 0.5 * (U2 - Zr * I2)          # same combination the modal S21 uses
s11r = ref1_r / inc1_r
s21r = trn2_r / inc1_r
s11r_dB = 20.0 * np.log10(np.abs(s11r))
s21r_dB = 20.0 * np.log10(np.abs(s21r))

# excitation content mask -- below this there is nothing to measure with
sig = np.abs(inc1_r) / np.abs(inc1_r).max()
have_signal = sig > SIG_FLOOR
fit_ok = have_signal & (cnd < COND_MAX) & (res < RES_MAX)

# the free consistency check: port-to-port is pure guide
d12 = Zz_m[idxPort2 - 1] - Zz_m[idxPort1 + 1]
s21_ana_dB = 20.0 * np.log10(np.abs(np.exp(-gam * d12)))

# --- (1)-independent ruler: the time gate -----------------------------------
d_round = 2.0 * (Zz.item(idxAbs2) - Zz.item(idxPort1)) * unit
t_peak = u1[np.argmax(np.abs(u1[:, 1])), 0]
t_gate = t_peak + 0.95 * d_round / C0
gam_dB = 20.0 * np.log10(gated_reflection(Sim_Path, 'port_ut_1', f, t_gate))
print('time gate at {:.2f} ns (incident peak {:.2f} ns + round trip {:.2f} ns)'.format(
      t_gate * 1e9, t_peak * 1e9, d_round / C0 * 1e9))

# Octave reference (WG_ModalMur_BoundaryTermination.m, pocket 0, off-DC source).
# Those numbers come from that script's own de-embedding at +-12 cells, so they
# are an order-of-magnitude orientation, not a bit-for-bit target.
f_tab = np.array([1.3, 1.5, 1.7, 1.9, 2.11, 2.30, 2.50, 2.70, 2.90]) * 1e9
oct_offdc = [-35.4, -46.8, -54.9, -63.2, -66.3, -70.4, -78.6, -75.0, -64.9]
oct_dc = [-37.1, -44.1, -63.4, -51.4, -52.6, -51.6, -49.9, -46.6, -42.3]
ref = oct_dc if DC_HEAVY else oct_offdc


def dbm(x, ok=True):
    return '     ---' if not ok else '{:8.2f}'.format(x)


print('\n===== ABOVE CUTOFF: absorber performance, {} source ====='.format(
      'DC-heavy' if DC_HEAVY else 'off-DC'))
print('  f[GHz]  octave-ref    GATED   MODAL-FIT   cond   port|S11|  |S21|dB')
for i, fi in enumerate(f_tab):
    j = np.abs(f - fi).argmin()
    print('   {:.2f}    {:8.2f}  {}  {}  {:5.1f}  {}  {}'.format(
          fi / 1e9, ref[i], dbm(gam_dB[j]),
          dbm(20 * np.log10(np.abs(Gam_fit[j])), fit_ok[j]), cnd[j],
          dbm(s11_dB[j]), dbm(s21_dB[j])))

print('\n===== BROADBAND, INCLUDING THE CUTOFF KNEE AND BELOW =====')
print('  |ZL|/Z0 shows the reference blowing up at cutoff -- that column is why')
print('  S11(ZL) cannot be trusted there, and why the other columns exist.')
print('  MODAL-FIT is |B/A| at the ladder: above cutoff the reflection, below')
print('  cutoff the growing-branch contamination the sheet should forbid.')
print('  S11(Zr) is referred to a real {:.1f} Ohm: bounded by 1 for any passive'.format(Zr))
print('  termination, so it is a passivity monitor, not an absorber figure.')
print('\n  f[GHz]  |ZL|/Z0  MODAL-FIT  cond  resid   1way    S11(ZL)   S11(Zr)   S21(Zr)  S21 ana')
k0_tab = 2.0 * np.pi * f / C0
for fi in [0.2, 0.6, 1.0, 1.15, 1.2, 1.25, 1.3, 1.5, 2.0, 2.5, 2.9]:
    j = np.abs(f - fi * 1e9).argmin()
    ZLr = np.abs(k0_tab[j] * Z0 / (gam[j] / 1j)) / Z0 if np.abs(gam[j]) > 0 else np.inf
    print('   {:.2f}  {:7.2f}  {} {:5.1f} {:6.3f} {:7.1f}  {}  {}  {}  {}'.format(
          fi, ZLr,
          dbm(20 * np.log10(np.abs(Gam_fit[j])), fit_ok[j]), cnd[j], res[j],
          decay_1way_dB[j],
          dbm(s11_dB[j], have_signal[j]), dbm(s11r_dB[j], have_signal[j]),
          dbm(s21r_dB[j], have_signal[j]), dbm(s21_ana_dB[j])))
if not DC_HEAVY:
    print('  (--- = masked: no excitation content, or the fit is not trustworthy')
    print('        there; set DC_HEAVY = True to actually drive the sub-cutoff band)')

band = (f >= 1.4e9) & (f <= 2.9e9)
print('\nworst GATED     |Gamma| over 1.4-2.9 GHz: {:.2f} dB   <-- independent ruler'.format(
      np.max(gam_dB[band])))
sel = band & fit_ok
if sel.any():
    print('worst MODAL-FIT |Gamma| over 1.4-2.9 GHz: {:.2f} dB   <-- E-only, no gate'.format(
          np.max(20 * np.log10(np.abs(Gam_fit[sel])))))
else:
    print('worst MODAL-FIT |Gamma| over 1.4-2.9 GHz: no usable points (cond > {:.0f})'.format(COND_MAX))
print('worst port |S11| (modal ZL) over 1.4-2.9 GHz: {:.2f} dB   <-- de-embedding floor'.format(
      np.max(s11_dB[band])))

# --- the three broadband checks ---------------------------------------------
#
# (a) Energy balance against the modal reference, propagating band only. Below
#     cutoff the reference is reactive and this is not defined.
bal = (np.abs(s11) ** 2 + np.abs(s21) ** 2)[band]
print('\nenergy balance, modal ZL ref, 1.4-2.9 GHz : min {:.3f}  max {:.3f}{}'.format(
      bal.min(), bal.max(), '   <-- >1: de-embed limited' if bal.max() > 1.02 else ''))
#
# (b) Passivity, defined everywhere. |S11| against a REAL reference impedance
#     exceeds 1 if and only if Re(Zin) < 0, i.e. if and only if the terminated
#     guide is delivering energy to the port. That is the one thing a boundary
#     condition must never do, it is a yes/no question, and unlike a reflection
#     coefficient it stays meaningful through the cutoff singularity and all
#     the way down to DC. (There is no matching statement for S21(Zr): with a
#     real reference the far port is not matched, so a2 != 0 and S21(Zr) is not
#     a proper S-parameter. It is printed only as a rough continuity check --
#     believe the analytic column next to it instead.)
#     Report it as the resistive FRACTION Re(Zin)/|Zin| rather than a raw
#     verdict. Below cutoff Zin is very nearly a pure reactance, so its real
#     part is a small difference of large numbers and the half-cell E/H offset
#     is enough to push it slightly negative; a fraction of -0.06 is three
#     degrees of phase on a 90-degree load, i.e. measurement noise. A boundary
#     that is genuinely pumping energy shows up as a fraction of order -1.
Zin1 = U1 / I1
resfrac = np.real(Zin1) / np.abs(Zin1)
#     Scanned over fit_ok, not merely have_signal: at 0.79 GHz the excitation
#     is 1e-3 of peak and the ladder residual is 0.5, i.e. what the port sees
#     there is not the TE10 mode at all. Reading a resistance out of that is
#     reading noise, and it duly reports -0.12.
jact = np.argmin(np.where(fit_ok, resfrac, np.inf))
print('passivity  min Re(Zin)/|Zin|: {:+.3f} at {:.3f} GHz  (|S11| vs real {:.0f} Ohm = {:.3f}){}'.format(
      resfrac[jact], f[jact] / 1e9, Zr, np.abs(s11r)[jact],
      '   <-- ACTIVE' if resfrac[jact] < -0.1 else '   (numerically passive)'))
#
# (c) The free one. Between the two port planes there is nothing but guide, so
#     S21 is analytically exp(-gamma*d) with no unknowns at all. Matching it
#     validates gamma, the mode files, the probes and the DFT in one number --
#     and it does so across the whole band, evanescent part included.
print('S21 vs analytic exp(-gamma*d), 1.4-2.9 GHz: max deviation {:.2f} dB'.format(
      np.max(np.abs(s21r_dB[band] - s21_ana_dB[band]))))
print('modal-fit residual over 1.4-2.9 GHz       : max {:.2e}  (two-exponential model)'.format(
      np.max(res[band])))

# ## Plots
figure()
plot(f / 1e9, gam_dB, 'r-', linewidth=2, label='$|\\Gamma|$ time-gated (independent ruler)')
gf = np.where(fit_ok, 20 * np.log10(np.abs(Gam_fit)), np.nan)
plot(f / 1e9, gf, 'm-', linewidth=2, label='$|\\Gamma|$ modal fit (E-only, stabilised)')
plot(f / 1e9, s11_dB, 'k--', linewidth=1.5, label='$S_{11}$ port, modal $Z_L$ ref')
plot(f / 1e9, s11r_dB, 'k-', linewidth=1.0, label='$S_{11}$ port, real $Z_r$ ref')
plot(f / 1e9, s21r_dB, 'b-', linewidth=2, label='$S_{21}$ (real $Z_r$ ref)')
plot(f / 1e9, s21_ana_dB, 'c:', linewidth=2, label='$S_{21}$ analytic $e^{-\\gamma d}$')
axvline(fc_abs / 1e9, color='g', linestyle=':', linewidth=1.5, label='$f_{cut}$')
grid()
legend(loc='lower right', fontsize=8)
ylim(-90, 10)
ylabel('S-Parameter (dB)')
xlabel('Frequency (GHz)')
title('Rect WG TE10 — Dispersive Modal Mur, stabilised post-processing')

figure()
semilogy(f / 1e9, cnd, 'k-', linewidth=2)
axhline(COND_MAX, color='r', linestyle='--', label='usable limit')
axvline(fc_abs / 1e9, color='g', linestyle=':', linewidth=1.5, label='$f_{cut}$')
grid()
legend()
ylabel('condition number of the forward/backward basis')
xlabel('Frequency (GHz)')
title('Where forward and backward can be told apart (analytic, data-free)')
show()
