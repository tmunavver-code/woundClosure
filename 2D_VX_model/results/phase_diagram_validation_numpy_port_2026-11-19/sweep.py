import sys, csv, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "/home/claude/phase")
from vxport import *

KA, DT, ETA = 1.0, 0.01, 1.0          # model defaults (SingleCellContractionParameters.m)
SOFT_TOL = 1e-3                        # soft if mean shape index exceeds hexagon's by > 1e-3


def relax(task):
    tag, Lbar, Gbar, nx, ny, nmax, seed, save = task
    p0 = -Lbar / (2.0 * Gbar)
    rstiff = 1.0 / (KA * Gbar)
    r, conn, Lx, Ly = hex_lattice(nx, ny)
    rng = np.random.default_rng(seed)
    r = r + 0.02 * rng.standard_normal(r.shape)
    r[:, 0] %= Lx; r[:, 1] %= Ly
    t0 = time.time()
    status = "maxsteps"
    for s in range(1, nmax + 1):
        f, A, P, e = forces(r, conn, Lx, Ly, KA, rstiff, p0)
        if not np.all(np.isfinite(f)):
            status = "nonfinite"; break
        fmax = np.abs(f).max()
        if fmax < 2e-6:
            status = "converged"; break
        r = update(r, f, DT, ETA, Lx, Ly)
    f, A, P, e = forces(r, conn, Lx, Ly, KA, rstiff, p0)
    q = np.mean(P / np.sqrt(np.maximum(A, 1e-12)))
    cv = e.std() / e.mean()
    epc = np.mean(KA * (A - 1) ** 2 + (P - p0) ** 2 / rstiff)
    if A.min() < 0.3 or status == "nonfinite":
        phase = "collapse"
    elif q - C_HEX > SOFT_TOL:
        phase = "soft"
    else:
        phase = "hexagonal"
    row = dict(set=tag, Lbar=Lbar, Gbar=Gbar, p0=p0, ka=KA, rstiff=rstiff, ncells=conn.shape[0],
               steps=s, status=status, q_mean=q, q_minus_c=q - C_HEX, edge_cv=cv,
               energy_per_cell=epc, A_min=A.min(), A_max=A.max(), max_force=np.abs(f).max(),
               phase=phase, seconds=time.time() - t0)
    if save:
        np.savez(f"/home/claude/phase/out/cfg_{save}.npz", r=r, conn=conn, Lx=Lx, Ly=Ly,
                 Lbar=Lbar, Gbar=Gbar, p0=p0)
    return row


def tasks():
    T = []
    # 1) along the line: dense p0 band at six Gbar values, 8x8 lattice
    for G in [0.02, 0.05, 0.08, 0.11, 0.14, 0.17]:
        for p0 in [3.55, 3.62, 3.67, 3.70, 3.71, 3.715, 3.72, 3.725, 3.73, 3.74, 3.76, 3.80, 3.88]:
            T.append(("line", -2 * p0 * G, G, 8, 8, 80000, 1, None))
    # 2) labelled points read off the figure
    for name, L, G in [("I", 0.12, 0.04), ("II", 0.0, 0.10), ("III", -0.85, 0.10),
                       ("IV", -0.32, 0.04), ("V", 0.0, 0.04)]:
        T.append(("point_" + name, L, G, 8, 8, 80000, 1, name))
    # 3) coarse grid over the whole diagram, 6x6 lattice
    for G in np.round(np.arange(0.01, 0.18, 0.02), 3):
        for L in np.round(np.arange(-1.4, 0.61, 0.1), 3):
            T.append(("grid", float(L), float(G), 6, 6, 30000, 1, None))
    # 4) fixed-box area instability: theory says Lbar > 8/c = 2.149, any Gbar
    for G in [0.05, 0.12]:
        for L in [1.6, 1.9, 2.05, 2.25, 2.4, 2.7]:
            T.append(("fixedbox_collapse", L, G, 6, 6, 40000, 1, None))
    # example networks for insets
    T.append(("inset_soft", -0.60, 0.05, 8, 8, 80000, 2, "soft_example"))
    return T


if __name__ == "__main__":
    import os
    os.makedirs("/home/claude/phase/out", exist_ok=True)
    T = tasks()
    # longest jobs first for better load balance
    T.sort(key=lambda t: -t[5] * t[3] * t[4])
    fields = None
    with open("/home/claude/phase/out/sweep_results.csv", "w", newline="") as fh, Pool(2) as pool:
        for i, row in enumerate(pool.imap_unordered(relax, T), 1):
            if fields is None:
                fields = list(row.keys())
                w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
            w.writerow(row); fh.flush()
            print(f"[{i}/{len(T)}] {row['set']:18s} Lbar={row['Lbar']:+.3f} Gbar={row['Gbar']:.3f} "
                  f"p0={row['p0']:+.3f} -> {row['phase']:9s} q-c={row['q_minus_c']:+.2e} "
                  f"{row['status']} {row['steps']} steps {row['seconds']:.0f}s", flush=True)
    print("DONE", flush=True)
