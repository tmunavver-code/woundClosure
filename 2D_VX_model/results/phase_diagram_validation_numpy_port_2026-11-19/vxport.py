"""
NumPy port of the 2D_VX_model mechanics (no wound terms), for a
fixed-topology mesh in which every cell has the same number of vertices.

Mirrors, formula for formula:
  getPolygonalCellArea.m  -> area = sum_k |triangle(R_k, R_k+1, Rp)|, Rp = vertex mean
  getAreaDerivative.m     -> the same centroid-triangulation derivative, with S/|S|
  getPerimeterDerivative.m
  getVertexForcesClassical.m  (bulk part only)
  getTissueEnergyClassical.m  (bulk part only)
  updateVertexPositions.m     -> r = r0 + dt*f/eta, then wrap into the box
  modifyVerticesForPeriodicity.m -> unwrap each cell relative to its first vertex
"""
import numpy as np

C_HEX = np.sqrt(8.0 * np.sqrt(3.0))  # shape index of a regular hexagon, 3.7224


def unwrap(r, conn, Lx, Ly):
    X = r[conn]                                   # (N, n, 2)
    d = X - X[:, :1, :]
    # modifyVerticesForPeriodicity: shift by L when beyond +-L/2 of vertex 1
    d[..., 0] = np.where(d[..., 0] > Lx / 2, d[..., 0] - Lx, d[..., 0])
    d[..., 0] = np.where(d[..., 0] < -Lx / 2, d[..., 0] + Lx, d[..., 0])
    d[..., 1] = np.where(d[..., 1] > Ly / 2, d[..., 1] - Ly, d[..., 1])
    d[..., 1] = np.where(d[..., 1] < -Ly / 2, d[..., 1] + Ly, d[..., 1])
    return X[:, :1, :] + d


def z(a, b):
    """z-component of (a,0) x (b,0)."""
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def zcross_vec(sig, v):
    """(sig*zhat) x (v,0)  ->  sig*(-v_y, v_x)."""
    return np.stack([-v[..., 1], v[..., 0]], axis=-1) * sig[..., None]


def vec_zcross(v, sig):
    """(v,0) x (sig*zhat)  ->  sig*(v_y, -v_x)."""
    return np.stack([v[..., 1], -v[..., 0]], axis=-1) * sig[..., None]


def geometry(r, conn, Lx, Ly):
    X = unwrap(r, conn, Lx, Ly)
    n = conn.shape[1]
    Rp = X.mean(axis=1, keepdims=True)
    Xn = np.roll(X, -1, axis=1)   # next
    Xp = np.roll(X, 1, axis=1)    # prev
    # getTriangleArea(p1=R_k, p2=R_k+1, p3=Rp) = 0.5*|(p2-p1) x (p3-p1)|
    A = 0.5 * np.abs(z(Xn - X, Rp - X)).sum(axis=1)
    e = np.linalg.norm(Xn - X, axis=-1)
    P = e.sum(axis=1)
    return X, Xn, Xp, Rp, A, P, e


def area_derivative(X, Xn, Xp, Rp):
    """Literal port of getAreaDerivative.m for every (cell, vertex)."""
    n = X.shape[1]
    # current & next: S = 0.5*(Rp - Rnext) x (Rcurr - Rnext); dA = 0.5*(S/|S|) x (Rp - Rnext)
    s1 = np.sign(0.5 * z(Rp - Xn, X - Xn))
    dA = 0.5 * zcross_vec(s1, Rp - Xn)
    # current & prev: S = 0.5*(Rcurr - Rprev) x (Rp - Rprev); dA += 0.5*(Rp - Rprev) x S/|S|
    s2 = np.sign(0.5 * z(X - Xp, Rp - Xp))
    dA = dA + 0.5 * vec_zcross(Rp - Xp, s2)
    # sum over all triangles of the Rp-dependence: 0.5/N * (Rv - Rvnext) x S_v/|S_v|
    sv = np.sign(0.5 * z(Rp - Xn, X - Xn))                      # per triangle v
    term = (0.5 / n) * vec_zcross(X - Xn, sv)                   # (N, n, 2)
    dA = dA + term.sum(axis=1, keepdims=True)
    return dA


def perimeter_derivative(X, Xn, Xp):
    a = X - Xn
    b = X - Xp
    return a / np.linalg.norm(a, axis=-1, keepdims=True) + \
        b / np.linalg.norm(b, axis=-1, keepdims=True)


def forces(r, conn, Lx, Ly, ka, rstiff, p0):
    X, Xn, Xp, Rp, A, P, e = geometry(r, conn, Lx, Ly)
    dA = area_derivative(X, Xn, Xp, Rp)
    dP = perimeter_derivative(X, Xn, Xp)
    fc = -ka * 2.0 * (A - 1.0)[:, None, None] * dA \
         - 2.0 / rstiff * (P - p0)[:, None, None] * dP
    f = np.zeros_like(r)
    np.add.at(f, conn.ravel(), fc.reshape(-1, 2))
    return f, A, P, e


def energy(r, conn, Lx, Ly, ka, rstiff, p0):
    _, _, _, _, A, P, _ = geometry(r, conn, Lx, Ly)
    return np.sum(ka * (A - 1.0) ** 2 + (1.0 / rstiff) * (P - p0) ** 2)


def update(r, f, dt, eta, Lx, Ly):
    r = r + dt * f / eta
    r[:, 0] = np.where(r[:, 0] > Lx, r[:, 0] - Lx, r[:, 0])
    r[:, 0] = np.where(r[:, 0] < 0, r[:, 0] + Lx, r[:, 0])
    r[:, 1] = np.where(r[:, 1] > Ly, r[:, 1] - Ly, r[:, 1])
    r[:, 1] = np.where(r[:, 1] < 0, r[:, 1] + Ly, r[:, 1])
    return r


def hex_lattice(nx, ny):
    """Periodic pointy-top hexagonal tiling, every cell area exactly 1.
    Returns r (V,2), conn (N,6) counter-clockwise, Lx, Ly."""
    assert ny % 2 == 0
    a = np.sqrt(2.0 / (3.0 * np.sqrt(3.0)))
    w, h = np.sqrt(3.0) * a, 1.5 * a
    Lx, Ly = nx * w, ny * h
    keys, verts, conn = {}, [], []
    ang = np.deg2rad(30 + 60 * np.arange(6))      # CCW from 30 deg
    for j in range(ny):
        for i in range(nx):
            cx = (i + 0.5 * (j % 2) + 0.5) * w
            cy = (j + 0.5) * h
            ids = []
            for t in ang:
                vx = (cx + a * np.cos(t)) % Lx
                vy = (cy + a * np.sin(t)) % Ly
                k = (round(vx / Lx * 1e6) % 1000000, round(vy / Ly * 1e6) % 1000000)
                if k not in keys:
                    keys[k] = len(verts)
                    verts.append([vx, vy])
                ids.append(keys[k])
            conn.append(ids)
    return np.array(verts), np.array(conn), Lx, Ly
