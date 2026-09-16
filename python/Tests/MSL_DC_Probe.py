"""Does the static buildup at the excitation survive a PML-terminated line?

This is a single-purpose diagnostic, not an absorber test. It answers one
question: is the linearly-growing static E field at the port plane caused by
the absorber sheets, or by the excitation itself?

THE SETUP IS CHOSEN TO EXONERATE EVERYTHING EXCEPT THE PORT
  * NO absorber sheets at all.
  * PML_8 on both z faces instead of PEC. That removes the DC short that the
    PEC end walls put across the line, and PML absorbs down to DC.
  * The line is twice as long (240 mm) and the excitation sits on the EXACT
    middle plane, so the two ends are symmetric and far away.
  * The cross-section, the x/y mesh and dz are byte-identical to
    MSL_W_ModalMur.py, so MSL_E.csv / MSL_H.csv apply unchanged.

If the ramp is still there, the sheets are innocent.

It also settles a second question for free. In MSL_W_ModalMur.py the field dump
lands on the Nyquist interval (140 timesteps = 24.02 GHz), so anything
oscillating at exactly that rate would alias to a constant and masquerade as a
static field. MSL_DC_OVS (default 3) oversamples the dump, which moves the alias
frequency. A real DC component does not care; an aliased one does.

Usage:
    python3 MSL_DC_Probe.py                 # PML ends, port excitation
    MSL_DC_ENDS=PEC python3 MSL_DC_Probe.py # A/B the end condition
    MSL_DC_OVS=1 python3 MSL_DC_Probe.py    # back to the Nyquist dump rate
"""

import os, tempfile, shutil, json
import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

ENDS = os.environ.get('MSL_DC_ENDS', 'PML')      # 'PML' or 'PEC'
OVS = int(os.environ.get('MSL_DC_OVS', '3'))     # dump oversampling vs Nyquist
NRTS = int(os.environ.get('MSL_DC_NRTS', '60000'))
ENDCRIT = float(os.environ.get('MSL_DC_ENDCRIT', '1e-30'))   # off by default: a tiny
#  DC field carries no energy, so the -50 dB criterion stops the run before the
#  ramp is visible. This diagnostic needs wall-clock time, not convergence.

Sim_Path = os.path.join(tempfile.gettempdir(), 'Test_MSL_DC_%s_ovs%d_zm%d' % (ENDS, OVS, int(os.environ.get('MSL_DC_ZMEAN','1'))))
shutil.rmtree(Sim_Path, ignore_errors=True)
os.makedirs(Sim_Path)

# ## Cross-section -- IDENTICAL to MSL_W_ModalMur.py, so the mode files fit
microstrip_W = 1.875
substrate_epsR = 4.5
substrate_width = 15.0
substrate_thickness = 1.0
substrate_cells = 8
trace_cells = 7
cu_thick = 0.1
port_w_fact = 6.0
port_h_fact = 5.5
Airbox_Add = 12.5
unit = 1e-3

substrate_length = 240.0          # twice as long; the excitation goes in the middle

f0, fc_exc = 1.55e9, 1.45e9

FDTD = openEMS(NrTS=NRTS, EndCriteria=ENDCRIT)
FDTD.SetGaussExcite(f0, fc_exc)
ZMEAN = bool(int(os.environ.get('MSL_DC_ZMEAN', '1')))
FDTD.SetExciteZeroMean(ZMEAN)
FDTD.SetOverSampling(OVS)

zbc = 'PML_8' if ENDS == 'PML' else 'PEC'
FDTD.SetBoundaryCond(['PEC', 'PEC', 'PEC', 'MUR', zbc, zbc])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)
mesh_res = ((C0 / (f0 + fc_exc)) / unit) / 75

SimBox = np.array([-substrate_width * 0.5, substrate_width * 0.5,
                   -cu_thick, substrate_thickness * (1.0 + port_h_fact) + Airbox_Add,
                   0.0, substrate_length])
mesh.AddLine('x', SimBox[0:2])
mesh.AddLine('y', SimBox[2:4])
mesh.AddLine('z', SimBox[4:6])

line = CSX.AddMetal('cu_top')
line.AddBox(priority=20, start=[-microstrip_W / 2, substrate_thickness, 0.0],
            stop=[microstrip_W / 2, substrate_thickness + cu_thick, substrate_length])
mesh.AddLine('x', [-microstrip_W / 2, microstrip_W / 2])
mesh.AddLine('y', [substrate_thickness, substrate_thickness + cu_thick])
mesh.AddLine('x', np.linspace(-0.5 * microstrip_W, 0.5 * microstrip_W, trace_cells))

p_x = microstrip_W * (1.0 + port_w_fact) * 0.5
p_y1 = substrate_thickness * (1.0 + port_h_fact)
mesh.AddLine('x', [-p_x, p_x])
mesh.AddLine('y', [0.0, p_y1])

sub = CSX.AddMaterial('FR4', epsilon=substrate_epsR)
sub.AddBox(priority=2, start=[-substrate_width / 2, 0.0, 0.0],
           stop=[substrate_width / 2, substrate_thickness, substrate_length])
mesh.AddLine('y', np.linspace(0, substrate_thickness, substrate_cells + 1))

gnd = CSX.AddMetal('cu_bot')
gnd.AddBox(priority=10, start=[-substrate_width / 2, -cu_thick, 0.0],
           stop=[substrate_width / 2, 0.0, substrate_length])
mesh.AddLine('y', [-cu_thick, 0.0])

# Same dz as MSL_W_ModalMur.py (0.66667 mm), so the cell at the port is the same
dz_target = mesh_res / 2.0
nz_cells = int(round(substrate_length / dz_target))
mesh.AddLine('z', np.linspace(0.0, substrate_length, nz_cells + 1).tolist())
mesh.SmoothMeshLines('x', mesh_res, 1.4)
mesh.SmoothMeshLines('y', mesh_res, 1.4)

Zz = np.unique(np.asarray(mesh.GetLines('z')))
nz = len(Zz)
assert nz_cells % 2 == 0, 'need an even cell count so the middle is a mesh line'
idxExc = nz_cells // 2                               # the EXACT middle plane
print('z mesh: %d lines, dz = %.5f mm, excitation on z[%d] = %.4f mm (middle)'
      % (nz, np.diff(Zz).max(), idxExc, Zz[idxExc]))

with open('MSL_mode_params.json') as fh:
    _mp = json.load(fh)
assert abs(_mp['microstrip_W'] - microstrip_W) < 1e-9
shutil.copy('MSL_E.csv', os.path.join(Sim_Path, 'MSL_E.csv'))
shutil.copy('MSL_H.csv', os.path.join(Sim_Path, 'MSL_H.csv'))

box_lo = [-p_x, 0.0]
box_hi = [p_x, p_y1]

port = FDTD.AddWaveGuidePort(1, box_lo + [Zz.item(idxExc)], box_hi + [Zz.item(idxExc + 1)],
                             'z', E_file="MSL_E.csv", H_file="MSL_H.csv",
                             kc=0.0, excite=1, excite_type=0)

# Dump a slab around the excitation only -- the buildup was localised to two
# cells, and the full box at this length would be ~2 GB.
HALF = 12.0   # mm either side of the excitation plane
Et = CSX.AddDump('Et', file_type=0, dump_type=0, dump_mode=1)
Et.AddBox([SimBox[0], SimBox[2], Zz[idxExc] - HALF],
          [SimBox[1], SimBox[3], Zz[idxExc] + HALF])

print('ends = %s   oversampling = %d   NrTS = %d   Sim_Path = %s'
      % (zbc, OVS, NRTS, Sim_Path))
FDTD.Run(Sim_Path, verbose=3, cleanup=False)
