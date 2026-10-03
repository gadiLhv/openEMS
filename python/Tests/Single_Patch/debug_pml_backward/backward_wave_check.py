"""Does a discretised closed waveguide carry BACKWARD waves (phase and group
velocity opposite)? A PML -- openEMS's boundary UPML and the IPML alike --
amplifies those, and a short guide closed by a PML then grows without bound.

The guide's Yee lattice (as openEMS meshes it: an E edge is PEC if its centre is
in metal, the edge permittivity the mean over the quarter-cells around it) gets a
Bloch phase theta = kz*dz along z. Each mode's d(omega^2)/d(theta) follows from
Hellmann-Feynman (within degenerate subspaces); a negative one at theta > 0 is a
backward wave. Seconds per case.

Configuration: the CASES list below. No command-line arguments.
"""
import numpy as np
from yee_dispersion import edges_and_mask, curl_matrix, C0

# ============================== CONFIG =======================================
# (shape, transverse cell mm, eps_r); shapes: 'coax' (pin r 1.12, outer r 4.0),
# 'tube' (no pin), 'rect' (8 x 8 mm). All in a 12 x 12 mm PEC box.
CASES = [('coax', 1.0, 2.2), ('coax', 0.75, 2.2), ('coax', 0.6, 2.2), ('coax', 0.5, 2.2),
         ('tube', 1.0, 2.2), ('rect', 2.0, 2.2), ('rect', 1.0, 2.2)]
DZ = 0.245                  # mm, cell along the guide
THETAS = (0.02, 0.05, 0.1)  # kz*dz
F_MAX = 150e9
# =============================================================================


def eps_at(x, y, shape, er):
    if shape == 'rect':
        return er if (abs(x) < 4e-3 and abs(y) < 4e-3) else 1.0
    return er if np.hypot(x, y) < 4e-3 else 1.0     # metal is not a material: the pin region is dielectric


def edge_eps(E, xs, h, shape, er):
    q = h / 4; out = []
    for c, i, j in E:
        if c == 'x': pts = [(0.5 * (xs[i] + xs[i + 1]), xs[j] + s * q) for s in (-1, 1)]
        elif c == 'y': pts = [(xs[i] + s * q, 0.5 * (xs[j] + xs[j + 1])) for s in (-1, 1)]
        else: pts = [(xs[i] + sx * q, xs[j] + sy * q) for sx in (-1, 1) for sy in (-1, 1)]
        out.append(np.mean([eps_at(x, y, shape, er) for x, y in pts]))
    return np.array(out)


def backward_modes(shape, h, dz, theta, er, half=6.0):
    xs, n, E = edges_and_mask(h * 1e-3, half * 1e-3, shape, 1.12e-3, 4.0e-3)
    M = np.diag(1 / np.sqrt(edge_eps(E, xs, h * 1e-3, shape, er)))
    C = curl_matrix(n, E, h * 1e-3, dz * 1e-3, theta); d = 1e-6
    dC = (curl_matrix(n, E, h * 1e-3, dz * 1e-3, theta + d) - curl_matrix(n, E, h * 1e-3, dz * 1e-3, theta - d)) / (2 * d)
    A = M @ C.conj().T @ C @ M; dA = M @ (dC.conj().T @ C + C.conj().T @ dC) @ M
    lam, V = np.linalg.eigh(A); scale = (C0 / (2 * np.pi)) ** 2; out = []; k = 0
    while k < len(lam):
        k2 = k + 1
        while k2 < len(lam) and abs(lam[k2] - lam[k]) <= 1e-7 * max(lam[k], 1e-30) + 1e-9 * lam.max(): k2 += 1
        f = np.sqrt(max(lam[k], 0) * scale)
        if lam[k] > 1e-8 * lam.max() and f < F_MAX:
            Vs = V[:, k:k2]
            for gi in np.linalg.eigvalsh(Vs.conj().T @ dA @ Vs):
                if gi < 0: out.append((f, gi / lam[k]))
        k = k2
    return out


for shape, h, er in CASES:
    for th in THETAS:
        b = backward_modes(shape, h, DZ, th, er)
        print('%-4s h=%.2f er=%.1f theta=%.2f: %d backward%s' % (shape, h, er, th, len(b),
              ''.join('  [%.2f GHz, dln(w2)/dth %.3g]' % (f / 1e9, g) for f, g in sorted(b)[:5])))
    print()
