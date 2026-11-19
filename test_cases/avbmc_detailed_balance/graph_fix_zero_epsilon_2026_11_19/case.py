#!/usr/bin/env python3
"""Prepare and analyze a short zero-epsilon GPU AVBMC diagnostic.

The production LJ kernels retain a hard overlap cutoff of 0.5 even when
epsilon is zero, so this is a near-ideal runtime check, not an exact ideal gas.
"""
import argparse
import math
import random
from pathlib import Path

N = 80
L = 12.0
RBIND = 2.4
ROOT = Path(__file__).resolve().parent


def distance(a, b):
    d = [abs(x - y) for x, y in zip(a, b)]
    return math.sqrt(sum(min(x, L - x) ** 2 for x in d))


def prepare():
    if int(L / (1.10 * 2.5)) < 4:
        raise ValueError("Box too small for four checkerboard cells per axis")
    rng = random.Random(20261109)
    points = []
    while len(points) < N:
        candidate = tuple(rng.random() * L for _ in range(3))
        if all(distance(candidate, p) > 0.6 for p in points):
            points.append(candidate)
    with (ROOT / "data.atoms").open("w") as f:
        f.write("Zero-epsilon AVBMC diagnostic\n\n")
        f.write(f"{N} atoms\n1 atom types\n\n")
        for axis in "xyz":
            f.write(f"{-L / 2:.8f} {L / 2:.8f} {axis}lo {axis}hi\n")
        f.write("\nAtoms # atomic\n\n")
        for i, point in enumerate(points, 1):
            f.write(f"{i} 1 " + " ".join(f"{x - L / 2:.8f}" for x in point) + "\n")
    (ROOT / "datos.nml").write_text("""&Control_Params
 istep_ini=0, istep_fin=5000, Neq=1000, Nmove=1,
 Nsave=100, Ndump=10, Nrestart=5000,
 traj_format='plain', data_file='data.atoms', Npart_types=1,
 ncluster=1, rcl=2.1, minPts=3, cluster_threshold=6,
 asym_threshold=0.5, cluster_types=1,
 avbmc=T, avbmc_k_trials=1, avbmc_max_trials=20,
 avbmc_p_in=0.35, avbmc_rbind=2.4, lammps=F
/
&MC_Params
 hmax=0.05, omax=0.5, vmax=0 0 0,
 displ_update=F, temp0=1, npt=F, seed=20261109,
 swap_moves=F
/
&Potential_Params
 model='LJ', rangepp=2.5, table_mc=F
/
--- LJ SIGMA MATRIX ---
 1.000
--- LJ EPSILON MATRIX ---
 0.000
""")


def analyze():
    lines = (ROOT / "movie.xyz").read_text().splitlines()
    counts_one = []
    counts_all = []
    for start in range(0, len(lines), N + 2):
        if int(lines[start]) != N:
            raise RuntimeError("Unexpected trajectory frame size")
        # movie.xyz stores reduced coordinates in [0,1); the binding radius is physical.
        pos = [tuple(float(x) * L for x in lines[start + 2 + i].split()[1:4]) for i in range(N)]
        counts_one.append(sum(distance(pos[0], pos[j]) < RBIND for j in range(1, N)))
        pairs = sum(distance(pos[i], pos[j]) < RBIND for i in range(N) for j in range(i))
        counts_all.append(2 * pairs / N)
    expected = (N - 1) * (4 * math.pi * RBIND**3 / 3) / L**3
    mean_one = sum(counts_one) / len(counts_one)
    mean_all = sum(counts_all) / len(counts_all)
    print(f"frames={len(counts_all)} one_target_Nin={mean_one:.5f} all_target_Nin={mean_all:.5f} ideal_limit={expected:.5f}")
    print("Note: the production LJ kernel still excludes r<0.5 at epsilon=0.")
    if abs(mean_all - expected) > 0.15:
        raise SystemExit("Near-ideal occupancy outside diagnostic tolerance")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["prepare", "analyze"])
    action = parser.parse_args().action
    prepare() if action == "prepare" else analyze()
