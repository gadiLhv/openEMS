"""Dispersion of the DISCRETE (Yee, continuous-time) closed waveguide, air-filled.

One z-period of the Yee lattice with a Bloch phase theta = kz*dz. E edges whose
centre lies in metal are PEC (as openEMS's CalcPEC: the Yee coordinate of the
edge). omega^2 = eig(C^H C)/(mu0 eps0). A branch with d(omega)/d(theta) < 0 at
small theta > 0 is a backward wave -- the case where a PML amplifies.
"""
import sys, numpy as np

C0 = 299792458.0


def edges_and_mask(h, half, shape, r_in=1.12, r_out=4.0):
    xs = np.arange(-half, half + 1e-9, h); n = len(xs)
    def metal(x, y):
        r = np.hypot(x, y)
        if shape == 'coax':
            return (r < r_in) or (r > r_out * (1 + 1e-12))
        if shape == 'tube':
            return r > r_out * (1 + 1e-12)
        if shape == 'rect':
            return (abs(x) > r_out + 1e-12) or (abs(y) > r_out + 1e-12)
    E = []          # (comp, i, j) ; Ex at (i+1/2, j), Ey at (i, j+1/2), Ez at (i, j)
    for i in range(n - 1):
        for j in range(n):
            if not metal(0.5 * (xs[i] + xs[i + 1]), xs[j]): E.append(('x', i, j))
    for i in range(n):
        for j in range(n - 1):
            if not metal(xs[i], 0.5 * (xs[j] + xs[j + 1])): E.append(('y', i, j))
    for i in range(n):
        for j in range(n):
            if not metal(xs[i], xs[j]): E.append(('z', i, j))
    return xs, n, E


def curl_matrix(n, E, h, dz, theta):
    idx = {e: k for k, e in enumerate(E)}
    ph = np.exp(-1j * theta) - 1.0          # forward difference along z, Bloch
    rows = []
    def ev(c, i, j): return idx.get((c, i, j))
    # Hx at (i, j+1/2, k+1/2): dEz/dy - dEy/dz
    for i in range(n):
        for j in range(n - 1):
            r = {}
            for k, s in ((ev('z', i, j + 1), 1 / h), (ev('z', i, j), -1 / h), (ev('y', i, j), -ph / dz)):
                if k is not None: r[k] = r.get(k, 0) + s
            rows.append(r)
    # Hy at (i+1/2, j, k+1/2): dEx/dz - dEz/dx
    for i in range(n - 1):
        for j in range(n):
            r = {}
            for k, s in ((ev('x', i, j), ph / dz), (ev('z', i + 1, j), -1 / h), (ev('z', i, j), 1 / h)):
                if k is not None: r[k] = r.get(k, 0) + s
            rows.append(r)
    # Hz at (i+1/2, j+1/2, k): dEy/dx - dEx/dy
    for i in range(n - 1):
        for j in range(n - 1):
            r = {}
            for k, s in ((ev('y', i + 1, j), 1 / h), (ev('y', i, j), -1 / h), (ev('x', i, j + 1), -1 / h), (ev('x', i, j), 1 / h)):
                if k is not None: r[k] = r.get(k, 0) + s
            rows.append(r)
    C = np.zeros((len(rows), len(E)), complex)
    for a, r in enumerate(rows):
        for k, s in r.items(): C[a, k] = s
    return C


def omegas(h, dz, theta, shape, half=6.0):
    xs, n, E = edges_and_mask(h * 1e-3, half * 1e-3, shape, 1.12e-3, 4.0e-3)
    C = curl_matrix(n, E, h * 1e-3, dz * 1e-3, theta)
    lam = np.linalg.eigvalsh(C.conj().T @ C)
    lam = lam[lam > 1e-6 * lam.max()]
    return np.sort(np.sqrt(lam) * C0 / (2 * np.pi))       # Hz


if __name__ == '__main__':
    shape = sys.argv[1]; h = float(sys.argv[2]); dz = float(sys.argv[3])
    fmax = float(sys.argv[4]) if len(sys.argv) > 4 else 80e9
    th = np.linspace(0, 0.6, 61)
    F = [omegas(h, dz, t, shape) for t in th]
    nb = min(len(f) for f in F)
    F = np.array([f[:nb] for f in F])
    # "backward": a branch (by index, sorted) whose frequency falls as theta rises from 0
    print('%s h=%.3g dz=%.3g: %d branches below %.0f GHz at theta=0' % (shape, h, dz, (F[0] < fmax).sum(), fmax / 1e9))
    for b in range(nb):
        if F[0, b] > fmax: break
        d = np.diff(F[:12, b])
        tag = 'BACKWARD' if d[0] < -1e-6 * F[0, b] else ''
        print('  branch %2d  f(theta=0) %7.3f GHz  f(theta=0.1) %7.3f GHz  %s' % (b, F[0, b] / 1e9, F[10, b] / 1e9, tag))
