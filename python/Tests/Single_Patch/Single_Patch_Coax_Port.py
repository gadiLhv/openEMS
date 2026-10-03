"""Single patch antenna (Blit's "Single Patch.step"), fed by a coax through the
ground plate, with a waveguide port + invisible PML on the coax input.

Geometry, straight from the STEP (mm). The solids are read in the STEP's order:
  1  substrate, z -1.575..0, 140 x 140, PLUS a plug of radius 4.0 at (0, 9)
     that fills the hole in the plate down to z = -2.8: the coax dielectric
  2  ground plate, z -2.8..-1.575, with a radius-4.0 hole at (0, 9): the coax
     outer conductor is that hole's wall
  3  centre pin, radius 1.12, z -2.8..0, up to the patch
  4  patch 39 x 39, two opposite corners cut 3.7 mm, z 0..0.035
  5-12 parasitic strips on the patch layer, z 0..0.035
The cylinders are 32-gons in the STEP itself. Coax: Z0 ~ 60/sqrt(2.2) ln(4/1.12)
~ 51.5 Ohm.

Every solid is exported to an STL and imported as a polyhedron, as Blit hands
openEMS its geometry. The only change: metal thinner than THIN_METAL_MAX
(the 35 um patch layer) becomes a zero-thickness polygon on its bottom face,
so the mesh does not need a 35 um cell.

Port 1: on the coax opening, z = -2.8 (bottom of the plate), guide on +z. A
square window, PORT_MARGIN wider than the coax on every side. An IPML sheet
covers the window, with a PEC block of the same window one mesh cell deep behind
it, at least one cell clear of the boundary PML.

The port mode is solved here, Blit's way: the STLs are sliced inside the coax
section (stls2geo.genPortGeoFromSTLs), solved with sparselizard
(blit_mode_solver), and sampled on this script's own mesh lines
(yee_consistent_export). Solver by-products go to ./work.

Configuration: the block below. No command-line arguments.
"""
import os, sys, json, shutil
import numpy as np

# ============================== CONFIG =======================================
STEP_FILE = '/home/gadi/Desktop/Vbox_Share/blit_port_mode_solver/Single Patch.step'
MODE_SOLVER_DIR = '/home/gadi/Repositories/modeSolver'  # blit_mode_solver, stls2geo, yee_consistent_export

DIELECTRIC_SOLIDS = [1]  # STEP solid numbers that are substrate; every other solid is PEC
SUBSTRATE_EPSR = 2.2
SUBSTRATE_TAND = 0.0009  # turned into a conductivity at F_RES (see below)
THIN_METAL_MAX = 0.1  # mm; thinner metal solids become sheets on their bottom face

F_RES = 2.495e9  # expected patch resonance; also the mode-solve frequency
F_MIN, F_MAX = 1.5e9, 3.5e9  # excitation band

PORT_Z = -2.8  # the coax opening (bottom of the plate)
COAX_CENTER = (0.0, 9.0)  # (x, y)
COAX_R_OUTER = 4.0
PORT_MARGIN = 2.0  # mm, window beyond the coax outer conductor, each side

# mesh
PORT_MESH = 'refined'  # 'original' (1.0 mm: the coax carries backward-wave modes that the PML amplifies, diverges at ~33 GHz) or 'refined'
ORIGINAL_RES = 1.0  # mm, uniform cells across the port window (x and y) for PORT_MESH = 'original'
REFINED_RES = 0.5  # mm, the same for PORT_MESH = 'refined' (stable; see debug_pml_backward/)
assert PORT_MESH in ('original', 'refined')
PORT_RES = ORIGINAL_RES if PORT_MESH == 'original' else REFINED_RES
N_COAX_CELLS = 5  # z cells in the coax section (-2.8..-1.575)
N_SUB_CELLS = 4  # z cells in the substrate (-1.575..0)
AIR_MARGIN = 30.0  # mm of air beyond the board on every side (~lambda/4 at F_RES)
MAX_RATIO = 1.4  # mesh grading ratio
N_BOUNDARY_PML = 8  # PML_8 on all six faces

# mode solve
MODE_MESH_SIZE = 0.25  # mm, genPortGeoFromSTLs characteristic length
MODE_MESH_FACTOR = 0.5  # gmsh MeshSizeFactor on top of that (DEFAULT_MESH_FACTOR)
NUM_MODES = 1
SHOW_MODE_GUI = False  # open gmsh on the solved mode

# outputs
SHOW_CSXCAD = False  # write the XML and open AppCSXCAD
DUMP_FIELDS = False  # E-field dumps (only written by an FDTD run)
DUMP_PATCH_BOX = True  # time-domain E dump of the patch's bounding box (independent of DUMP_FIELDS)
PATCH_DUMP_MARGIN = 10.0  # mm beyond the patch in x, y and above it; below, it reaches the port plane
PATCH_DUMP_CELLS_BELOW = 6  # mesh cells below the port plane (through the PEC block, into the air under the plate)
RUN_FDTD = True  # the mode solver is under test, not the FDTD run
# =============================================================================

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, 'work')
SIM_PATH = os.path.join(WORK, 'sim')
E_FILE = os.path.join(HERE, 'coax_port_E.csv')
H_FILE = os.path.join(HERE, 'coax_port_H.csv')
PARAMS_FILE = os.path.join(HERE, 'coax_port_mode_params.json')
C0 = 299792458.0
EPS0 = 8.8541878128e-12
os.makedirs(WORK, exist_ok=True)


# ---- STEP -> one STL per solid, and sheets for the thin metal ---------------
def step_to_solids(step_file, out_dir):
    import gmsh, trimesh
    gmsh.initialize()
    gmsh.option.setNumber('General.Verbosity', 1)
    gmsh.model.occ.importShapes(step_file)
    gmsh.model.occ.synchronize()
    # Triangulate with the facet corners only: every face of this STEP is planar,
    # and vertices inserted along an edge put whole rows of STL vertices at
    # arbitrary heights, where a slice plane can land on them.
    gmsh.option.setNumber('Mesh.MeshSizeMin', 1e3)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 1e3)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.option.setNumber('Mesh.MinimumCurveNodes', 2)  # default 3 puts a node mid-edge
    gmsh.model.mesh.generate(2)
    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    xyz = dict(zip(ntags, ncoords.reshape(-1, 3)))

    solids = []
    for _, v in gmsh.model.getEntities(3):
        faces = [t for _, t in gmsh.model.getBoundary([(3, v)], oriented=False)]
        tris = []
        for f in faces:
            etypes, _, enodes = gmsh.model.mesh.getElements(2, f)
            for et, en in zip(etypes, enodes):
                assert et == 2, 'expected triangles'
                tris.append(np.asarray(en).reshape(-1, 3))
        tris = np.vstack(tris)
        tm = trimesh.Trimesh(vertices=np.array([xyz[n] for n in tris.ravel()]),
                             faces=np.arange(tris.size).reshape(-1, 3), process=True)
        trimesh.repair.fix_normals(tm)
        assert tm.is_watertight, 'solid %d is not watertight' % v
        stl = os.path.join(out_dir, 'solid_%02d.stl' % v)
        tm.export(stl)
        # gmsh pads its bounding boxes; the STL's own vertices do not
        bb = list(tm.bounds[0]) + list(tm.bounds[1])

        sheet = None
        if bb[5] - bb[2] < THIN_METAL_MAX:
            # the bottom face's outline, in order (gmsh pads bounding boxes,
            # so: the flat face with the lowest z)
            flat = [(gmsh.model.getBoundingBox(2, f), f) for f in faces]
            flat = [(b[2], f) for b, f in flat if b[5] - b[2] < 1e-6]
            bottom = min(flat)[1]
            loop = gmsh.model.getBoundary([(2, bottom)], oriented=True)
            pts = []
            for _, c in loop:
                p = [t for _, t in gmsh.model.getBoundary([(1, abs(c))], oriented=False)]
                if c < 0:
                    p = p[::-1]
                pts.append(gmsh.model.getValue(0, p[0], []))
            pts = np.array(pts)
            sheet = dict(z=float(np.mean(pts[:, 2])), xy=pts[:,:2])
        solids.append(dict(nr=v, bbox=bb, stl=stl, sheet=sheet,
                           dielectric=v in DIELECTRIC_SOLIDS))
    gmsh.finalize()
    return solids


solids = step_to_solids(STEP_FILE, WORK)
print('STEP solids:')
for s in solids:
    bb = s['bbox']
    print('  %2d  %-10s x %+8.3f..%+8.3f  y %+8.3f..%+8.3f  z %+7.3f..%+7.3f%s'
          % (s['nr'], 'dielectric' if s['dielectric'] else 'PEC', bb[0], bb[3], bb[1], bb[4], bb[2], bb[5],
             '  -> sheet at z=%g' % s['sheet']['z'] if s['sheet'] is not None else ''))
assert not any(s['dielectric'] and s['sheet'] is not None for s in solids), 'a thin dielectric?'

# ---- mesh -------------------------------------------------------------------
from CSXCAD import ContinuousStructure, AppCSXCAD_BIN
from CSXCAD.CSProperties import ABCtype
from CSXCAD.SmoothMeshLines import SmoothMeshLines
from openEMS import openEMS

bb_all = np.array([s['bbox'] for s in solids])
board = (bb_all[:, 0].min(), bb_all[:, 3].max(), bb_all[:, 1].min(), bb_all[:, 4].max(),
         bb_all[:, 2].min(), bb_all[:, 5].max())
half_win = COAX_R_OUTER + PORT_MARGIN
win_x = (COAX_CENTER[0] - half_win, COAX_CENTER[0] + half_win)
win_y = (COAX_CENTER[1] - half_win, COAX_CENTER[1] + half_win)

lam_min = C0 / F_MAX * 1e3  # mm
max_res = lam_min / np.sqrt(SUBSTRATE_EPSR) / 20.0  # everywhere: the board fills most of it
print('max cell %.2f mm (lambda/20 in the substrate at %.2f GHz)' % (max_res, F_MAX / 1e9))


def xy_lines(axis, win):
    k = 0 if axis == 'x' else 1
    lines = [board[2 * k] - AIR_MARGIN, board[2 * k + 1] + AIR_MARGIN]
    for s in solids:  # box edges of every solid
        lines += [s['bbox'][k], s['bbox'][k + 3]]
        if s['sheet'] is not None:  # and the sheet outlines (patch corners)
            lines += list(s['sheet']['xy'][:, k])
    n = int(round((win[1] - win[0]) / PORT_RES))
    fine = list(np.linspace(win[0], win[1], n + 1))
    # inside the window only the uniform port grid
    lines = [l for l in lines if not (win[0] - 0.5 * PORT_RES < l < win[1] + 0.5 * PORT_RES)] + fine
    lines = np.unique(np.round(lines, 6))
    return SmoothMeshLines(lines, max_res, MAX_RATIO)


coax_top = -1.575
for s in solids:  # the plate's top is the coax/substrate interface
    if not s['dielectric'] and s['bbox'][2] < PORT_Z + 1e-9 and s['bbox'][3] - s['bbox'][0] > 50:
        coax_top = s['bbox'][5]
dz_coax = (coax_top - PORT_Z) / N_COAX_CELLS
z_fixed = (list(np.linspace(PORT_Z, coax_top, N_COAX_CELLS + 1))
           +list(np.linspace(coax_top, 0.0, N_SUB_CELLS + 1))
           +[PORT_Z - dz_coax,  # the PEC block behind the port
              board[4] - AIR_MARGIN, board[5] + AIR_MARGIN])
z_fixed += [s['bbox'][5] for s in solids if s['sheet'] is None]
z_lines = SmoothMeshLines(np.unique(np.round(z_fixed, 6)), max_res, MAX_RATIO)

MESH_X, MESH_Y, MESH_Z = xy_lines('x', win_x), xy_lines('y', win_y), np.asarray(z_lines)
print('mesh %d x %d x %d = %.2f M cells' % (len(MESH_X), len(MESH_Y), len(MESH_Z),
                                            len(MESH_X) * len(MESH_Y) * len(MESH_Z) / 1e6))


def line_index(lines, v):
    return int(np.argmin(np.abs(np.asarray(lines) - v)))


for lines, v, nm in ((MESH_X, win_x[0], 'x0'), (MESH_X, win_x[1], 'x1'), (MESH_Y, win_y[0], 'y0'),
                     (MESH_Y, win_y[1], 'y1'), (MESH_Z, PORT_Z, 'port z')):
    assert abs(lines[line_index(lines, v)] - v) < 1e-6, 'port window %s is not a mesh line' % nm

# ---- port mode: slice the STLs inside the coax, solve, sample on the mesh ---
sys.path.insert(0, MODE_SOLVER_DIR)
from stls2geo import genPortGeoFromSTLs
from stls2geo.geomod import mergeGeoChains
from blit_mode_solver import MatDef
from blit_mode_solver import WaveguideBetaSolver as wgs
from yee_consistent_export import port_grid_coords, port_mode_arrays

# The slice is taken inside the coax section, where the cross-section is uniform,
# and never on a row of STL vertices: there the slice yields degenerate 2-point
# loops, and stls2geo drops the WHOLE STL with only a printed message -- the
# mode is then solved without it. So: the middle of the widest vertex-free gap.
import trimesh
cut = [s for s in solids if s['bbox'][2] < 0.5 * (PORT_Z + coax_top) < s['bbox'][5]]
zv = np.concatenate([trimesh.load(s['stl'], force='mesh').vertices[:, 2] for s in cut])
zv = np.unique(np.concatenate([zv[(zv > PORT_Z) & (zv < coax_top)], [PORT_Z, coax_top]]))
k = int(np.argmax(np.diff(zv)))
z_slice = 0.5 * (zv[k] + zv[k + 1])
stls = {os.path.basename(s['stl']): ('Dielectric' if s['dielectric'] else 'PEC') for s in cut}
print('mode solve: slice at z = %.4f (%.3f mm from the nearest STL vertex) through %s'
      % (z_slice, 0.5 * (zv[k + 1] - zv[k]), ', '.join(sorted(stls))))
materials = {"outer": MatDef(physIdx=0, matType="bounds"),
             "background": MatDef(physIdx=1, matType="normal"),
             "PEC": MatDef(physIdx=2, matType="pec"),
             "Dielectric": MatDef(physIdx=3, matType="normal", er=SUBSTRATE_EPSR)}
portDims = [win_x[0], win_x[1], win_y[0], win_y[1], z_slice, z_slice]

cwd = os.getcwd()
os.chdir(WORK)
genPortGeoFromSTLs(stls, materials, portDims, 3, minAth=1e-7, unitScale=1e-3,
                   geoFileName='coax_port', meshSize=MODE_MESH_SIZE)
# merge the sliced polygon chains (the 32-gon facets) to cut the FEM mesh down
mergeGeoChains('coax_port.geo', 'coax_port_geomod.geo')
# propDir +z: the geo's X / Y are openEMS x / y, the mode file's u / v
gridCoords, yee_meta = port_grid_coords(MESH_X, MESH_Y, win_x, win_y, unit=1e-3)
solver = wgs(geoFile='coax_port_geomod.geo', matDict=materials, gridCoords=gridCoords)
solver.generateMesh(meshFactor=MODE_MESH_FACTOR)
solver.setupSolver(f0=F_RES)
Er, Hr, f_cutoff, kz, modeType, Zw, Zl = solver.solveFields(
    propDir=2, epsr_eff=SUBSTRATE_EPSR, numModes=NUM_MODES, saveEHfieldFiles=True)
solver.finalize()
os.chdir(cwd)

# what an ideal coax of these radii would give, for the circles and the 32-gons
r_in = 1.12
n_ideal = np.sqrt(SUBSTRATE_EPSR)
Zl_circ = 60.0 / n_ideal * np.log(COAX_R_OUTER / r_in)
for i in range(len(kz)):
    n_eff = float(kz[i]) * C0 / (2 * np.pi * F_RES)
    print('\nmode %d: %s  f_cutoff %.4g GHz  n_eff %.5f (ideal TEM %.5f, %+.3f %%)'
          % (i, modeType[i], f_cutoff[i] / 1e9, n_eff, n_ideal, 100 * (n_eff / n_ideal - 1)))
    print('        Zw %.3f Ohm (TEM %.3f)   Zl %.3f Ohm (ideal circular coax %.3f)'
          % (Zw[i], 376.730313 / n_ideal, Zl[i], Zl_circ))

E_rows, H_rows = port_mode_arrays(Er, Hr, yee_meta, mode=0)
np.savetxt(E_FILE, E_rows, delimiter=',')
np.savetxt(H_FILE, H_rows, delimiter=',')
params = dict(mode_type=str(modeType[0]), f_cutoff=float(f_cutoff[0]), beta=float(kz[0]),
              n_eff=float(kz[0]) * C0 / (2 * np.pi * F_RES), Zw=float(Zw[0]), Zl=float(Zl[0]),
              f_solve=F_RES, z_slice=z_slice, win_x=win_x, win_y=win_y, port_z=PORT_Z,
              x_lines=[float(v) for v in MESH_X if win_x[0] - 1e-9 <= v <= win_x[1] + 1e-9],
              y_lines=[float(v) for v in MESH_Y if win_y[0] - 1e-9 <= v <= win_y[1] + 1e-9],
              mode_mesh_size=MODE_MESH_SIZE, mode_mesh_factor=MODE_MESH_FACTOR)
with open(PARAMS_FILE, 'w') as fh:
    json.dump(params, fh, indent=2)
print('wrote %s, %s (%d rows), %s' % (os.path.basename(E_FILE), os.path.basename(H_FILE),
                                      len(E_rows), os.path.basename(PARAMS_FILE)))

if SHOW_MODE_GUI:
    import gmsh
    gmsh.initialize()
    gmsh.open(os.path.join(WORK, 'coax_port_geomod.geo'))
    gmsh.merge(os.path.join(WORK, 'coax_port_geomod.msh'))
    for nm in ('E_r_coax_port_geomod_0.pos', 'H_r_coax_port_geomod_0.pos'):
        if os.path.exists(os.path.join(WORK, nm)):
            gmsh.merge(os.path.join(WORK, nm))
    gmsh.fltk.run()
    gmsh.finalize()

# ---- openEMS model ----------------------------------------------------------
FDTD = openEMS(EndCriteria=1e-4)
FDTD.SetGaussExcite(0.5 * (F_MIN + F_MAX), 0.5 * (F_MAX - F_MIN))
FDTD.SetExciteZeroMean(True)  # zero time-integral: the pulse leaves no static charge at the source
FDTD.SetBoundaryCond(['PML_%d' % N_BOUNDARY_PML] * 6)
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
grid = CSX.GetGrid()
grid.SetDeltaUnit(1e-3)
grid.SetLines('x', MESH_X); grid.SetLines('y', MESH_Y); grid.SetLines('z', MESH_Z)

# loss tangent as a conductivity at the resonance: kappa = 2 pi f eps0 epsr tand
kappa = 2 * np.pi * F_RES * EPS0 * SUBSTRATE_EPSR * SUBSTRATE_TAND
sub = CSX.AddMaterial('Substrate', epsilon=SUBSTRATE_EPSR, kappa=kappa)
pec = CSX.AddMetal('PEC')
for s in solids:
    if s['dielectric']:
        sub.AddPolyhedronReader(s['stl'], priority=1).ReadFile()
    elif s['sheet'] is None:
        pec.AddPolyhedronReader(s['stl'], priority=10).ReadFile()
    else:
        pec.AddPolygon(s['sheet']['xy'].T, 'z', s['sheet']['z'], priority=10)


def clear_of_boundary_pml(start, stop):
    # At least one mesh cell between the box and the boundary PML, on every face.
    for lines, lo, hi, axis in zip((MESH_X, MESH_Y, MESH_Z), start, stop, 'xyz'):
        i_lo, i_hi = line_index(lines, min(lo, hi)), line_index(lines, max(lo, hi))
        if i_lo < N_BOUNDARY_PML + 1 or i_hi > len(lines) - 1 - (N_BOUNDARY_PML + 1):
            raise ValueError('box %s..%s touches the boundary PML along %s' % (start, stop, axis))


# port 1: IPML sheet on the window, PEC block one cell deep behind it (-z),
# excitation one cell in front of the sheet, measurement two cells in
i0 = line_index(MESH_Z, PORT_Z)
ipml = CSX.AddAbsorbingBC('port_1_ipml', NormalSignPositive=True, AbsorbingBoundaryType=ABCtype.PML_8)
ipml.AddBox([win_x[0], win_y[0], PORT_Z], [win_x[1], win_y[1], PORT_Z], priority=60)
block = ([win_x[0], win_y[0], MESH_Z[i0 - 1]], [win_x[1], win_y[1], PORT_Z])
clear_of_boundary_pml(*block)
pec.AddBox(*block, priority=50)
port = FDTD.AddWaveGuidePort(1, [win_x[0], win_y[0], MESH_Z[i0 + 1]], [win_x[1], win_y[1], MESH_Z[i0 + 2]], 'z',
                             kc=0.0, excite=1.0, excite_type=0, E_file=E_FILE, H_file=H_FILE, mode_type='TEM')

if DUMP_FIELDS:
    # time domain on two planes (a full-volume time-domain dump would be GBs):
    # through the coax axis, and the middle of the substrate
    xz = CSX.AddDump('Et_xz_coax', dump_type=0, file_type=0)
    xz.AddBox([MESH_X[0], COAX_CENTER[1], MESH_Z[0]], [MESH_X[-1], COAX_CENTER[1], MESH_Z[-1]])
    xy = CSX.AddDump('Et_xy_substrate', dump_type=0, file_type=0)
    xy.AddBox([MESH_X[0], MESH_Y[0], 0.5 * coax_top], [MESH_X[-1], MESH_Y[-1], 0.5 * coax_top])
    # frequency domain, whole volume, at the expected resonance
    ef = CSX.AddDump('Ef_res', dump_type=10, file_type=1, frequency=[F_RES])
    ef.AddBox([MESH_X[0], MESH_Y[0], MESH_Z[0]], [MESH_X[-1], MESH_Y[-1], MESH_Z[-1]])

if DUMP_PATCH_BOX:
    # the patch: the thin metal sheet above the coax feed
    patch = [s for s in solids if s['sheet'] is not None
             and s['bbox'][0] <= COAX_CENTER[0] <= s['bbox'][3] and s['bbox'][1] <= COAX_CENTER[1] <= s['bbox'][4]]
    assert len(patch) == 1, 'expected one sheet above the coax, found %d' % len(patch)
    pb = patch[0]['bbox']
    z_bot = MESH_Z[max(line_index(MESH_Z, PORT_Z) - PATCH_DUMP_CELLS_BELOW, 0)]
    pd = CSX.AddDump('Et_patch_box', dump_type=0, file_type=0)
    pd.AddBox([pb[0] - PATCH_DUMP_MARGIN, pb[1] - PATCH_DUMP_MARGIN, z_bot],
              [pb[3] + PATCH_DUMP_MARGIN, pb[4] + PATCH_DUMP_MARGIN, pb[5] + PATCH_DUMP_MARGIN])
    print('patch dump box: x %+.2f..%+.2f  y %+.2f..%+.2f  z %+.3f..%+.3f mm'
          % (pb[0] - PATCH_DUMP_MARGIN, pb[3] + PATCH_DUMP_MARGIN, pb[1] - PATCH_DUMP_MARGIN,
             pb[4] + PATCH_DUMP_MARGIN, z_bot, pb[5] + PATCH_DUMP_MARGIN))

if os.path.exists(SIM_PATH):
    shutil.rmtree(SIM_PATH)
os.makedirs(SIM_PATH)
CSX_file = os.path.join(SIM_PATH, 'single_patch.xml')
CSX.Write2XML(CSX_file)
print('wrote %s' % CSX_file)
if SHOW_CSXCAD:
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if RUN_FDTD:
    FDTD.Run(SIM_PATH, cleanup=False)
