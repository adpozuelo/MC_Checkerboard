#!/usr/bin/env python3
"""
Comparison script: MC_Checkerboard (corrected dual-patch CSW) vs LAMMPS MD benchmarks.
Generates comprehensive comparative plots.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    mc_dir_54k = "/home/e.lomba/MC_Checkerboard/test_cases/sb_avbmc_50_50_3x3x3"
    mc_prev_54k = "/home/e.lomba/MC_Checkerboard/test_cases/sb_avbmc_50_50_3x3x3/prev_unlike_only"
    lmp_dir = "/home/e.lomba/sanchez_burgos_20260915/sanchez_burgos_lammps/runs/paper_a/seed_712345"
    out_dir = "/home/e.lomba/.gemini/antigravity-cli/brain/94fa1bcd-3919-4b15-ba00-2516a47d2736"

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # -------------------------------------------------------------
    # Panel 1: Energy Evolution
    # -------------------------------------------------------------
    ax = axes[0, 0]

    # Load Previous MC 54k (unlike-only)
    prev_data = os.path.join(mc_prev_54k, "run-data.dat")
    if os.path.exists(prev_data):
        rows = []
        with open(prev_data) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 8:
                    try:
                        rows.append([float(x) for x in parts])
                    except ValueError:
                        continue
        if rows:
            d_prev = np.array(rows)
            stride = max(1, len(d_prev) // 500)
            ax.plot(d_prev[::stride, 0], d_prev[::stride, 7], label=r"MC Unlike-only CSW ($N=54,000$)", color="gray", ls=":", lw=1.5)

    # Load New MC 54k (corrected dual-patch CSW)
    data_54k_file = os.path.join(mc_dir_54k, "run-data.dat")
    if os.path.exists(data_54k_file):
        rows = []
        with open(data_54k_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 8:
                    try:
                        rows.append([float(x) for x in parts])
                    except ValueError:
                        continue
        if rows:
            d54k = np.array(rows)
            step_54k = d54k[:, 0]
            e_54k = d54k[:, 7]
            stride = max(1, len(step_54k) // 1000)
            ax.plot(step_54k[::stride], e_54k[::stride], label=r"MC Dual CSW ($N=54,000$, corrected)", color="#1f77b4", lw=2.2)

    # LAMMPS reference levels
    ax.axhline(-15.27, color="#d62728", ls="--", lw=1.8, label=r"LAMMPS MD continuation ($\approx -15.27\,k_BT$)")
    ax.axhline(-15.12, color="#ff7f0e", ls="-.", lw=1.5, label=r"LAMMPS MD slab audit ($\approx -15.12\,k_BT$)")
    ax.axhline(-15.60, color="#1f77b4", ls="--", alpha=0.7, label=r"MC Dual CSW final ($\approx -15.60\,k_BT$)")

    ax.set_xlabel("MC Step", fontsize=12)
    ax.set_ylabel(r"Potential Energy per Molecule $\langle U/N \rangle$ [$k_B T$]", fontsize=11)
    ax.set_title("Energy Evolution: MC vs LAMMPS Benchmarks", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8.5, loc="center right")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-16.5, -13.5)

    # -------------------------------------------------------------
    # Panel 2: Radial Distribution Function g(r) Comparison
    # -------------------------------------------------------------
    ax = axes[0, 1]

    # MC 54k gmixsim.dat
    g_mc_file = os.path.join(mc_dir_54k, "gmixsim.dat")
    if os.path.exists(g_mc_file):
        g_mc = np.loadtxt(g_mc_file)
        r_mc = g_mc[:, 0]
        mask = (r_mc >= 0.8) & (r_mc <= 2.5)
        ax.plot(r_mc[mask], g_mc[mask, 4], label=r"MC: $g_{12}(r)$ (Scaffold-Client)", color="#1f77b4", lw=2.2)
        ax.plot(r_mc[mask], g_mc[mask, 3], label=r"MC: $g_{11}(r)$ (Scaffold-Scaffold, active)", color="#1f77b4", ls="--", lw=1.8)
        ax.plot(r_mc[mask], g_mc[mask, 5], label=r"MC: $g_{22}(r)$ (Client-Client, zero)", color="#1f77b4", ls=":", lw=1.5)

    # LAMMPS gmixsim.dat
    g_lmp_file = os.path.join(lmp_dir, "gmixsim.dat")
    if os.path.exists(g_lmp_file):
        g_lmp = np.loadtxt(g_lmp_file)
        r_lmp = g_lmp[:, 0]
        mask_lmp = (r_lmp >= 0.8) & (r_lmp <= 2.5)
        ax.plot(r_lmp[mask_lmp], g_lmp[mask_lmp, 4], label=r"LAMMPS: $g_{12}(r)$", color="#d62728", lw=1.8, alpha=0.8)
        ax.plot(r_lmp[mask_lmp], g_lmp[mask_lmp, 3], label=r"LAMMPS: $g_{11}(r)$ (Active CSW)", color="#ff7f0e", ls="--", lw=1.8, alpha=0.8)
        ax.plot(r_lmp[mask_lmp], g_lmp[mask_lmp, 5], label=r"LAMMPS: $g_{22}(r)$", color="#2ca02c", ls=":", lw=1.5, alpha=0.8)

    ax.set_xlabel(r"Separation $r / \sigma$", fontsize=12)
    ax.set_ylabel(r"$g_{ij}(r)$", fontsize=12)
    ax.set_title("Pair Correlations: Corrected Dual CSW vs LAMMPS", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8.5, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.8, 2.5)
    ax.set_ylim(0, 52)

    # -------------------------------------------------------------
    # Panel 3: Structure Factor S(Q)
    # -------------------------------------------------------------
    ax = axes[1, 0]

    sq_mc_file = os.path.join(mc_dir_54k, "sq.dat")
    if os.path.exists(sq_mc_file):
        sq_mc = np.loadtxt(sq_mc_file)
        q_mc = sq_mc[:, 0]
        s_nn_mc = sq_mc[:, 1]
        mask_q = (q_mc >= 0.2) & (q_mc <= 8.0)
        ax.semilogy(q_mc[mask_q], s_nn_mc[mask_q], "o-", ms=3, label=r"MC Dual CSW $N=54,000$ $S_{NN}(Q)$", color="#1f77b4")

    sq_lmp_file = os.path.join(lmp_dir, "sq.dat")
    if os.path.exists(sq_lmp_file):
        sq_lmp = np.loadtxt(sq_lmp_file)
        q_lmp = sq_lmp[:, 0]
        s_nn_lmp = sq_lmp[:, 1]
        mask_ql = (q_lmp >= 0.2) & (q_lmp <= 8.0)
        ax.semilogy(q_lmp[mask_ql], s_nn_lmp[mask_ql], "s--", ms=3, label=r"LAMMPS MD $N=1,820$ $S_{NN}(Q)$", color="#d62728", alpha=0.8)

    ax.set_xlabel(r"Wavevector $Q \sigma$", fontsize=12)
    ax.set_ylabel(r"Total Structure Factor $S_{NN}(Q)$ (log scale)", fontsize=11)
    ax.set_title("Structure Factor $S_{NN}(Q)$: Phase Separation Divergence", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9.5)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(0.2, 8.0)

    # -------------------------------------------------------------
    # Panel 4: Cluster Fraction & Composition Evolution
    # -------------------------------------------------------------
    ax = axes[1, 1]

    cl_mc_file = os.path.join(mc_dir_54k, "clusevol_mc.dat")
    if os.path.exists(cl_mc_file):
        cl_mc = np.loadtxt(cl_mc_file)
        step_cl = cl_mc[:, 0]
        pct_cl = cl_mc[:, 4]
        ax.plot(step_cl, pct_cl, label=r"MC Dual CSW $N=54,000$: % in Clusters", color="#1f77b4", lw=2)

    # LAMMPS clusevol
    cl_lmp_file = os.path.join(lmp_dir, "clusevol.dat")
    if os.path.exists(cl_lmp_file):
        cl_lmp = np.loadtxt(cl_lmp_file)
        conf_lmp = cl_lmp[:, 0]
        pct_lmp = cl_lmp[:, 2]
        stride = max(1, len(conf_lmp) // 200)
        ax.plot(conf_lmp[::stride] * 5000, pct_lmp[::stride], label=r"LAMMPS MD: % in Clusters", color="#d62728", ls="--", lw=2)

    ax.axhline(99.68, color="#1f77b4", ls=":", alpha=0.8, label=r"MC Dual CSW Plateau ($99.68\%$)")
    ax.axhline(94.6, color="#d62728", ls=":", alpha=0.8, label=r"LAMMPS Plateau ($\sim 94.6\%$)")

    ax.set_xlabel("Simulation Step / Scale", fontsize=12)
    ax.set_ylabel(r"Condensed Particle Fraction [%]", fontsize=12)
    ax.set_title("Condensate Formation: Cluster Size Fraction", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(85, 102)

    plt.tight_layout()
    plot_path = os.path.join(out_dir, "mc_vs_lammps_corrected.png")
    plt.savefig(plot_path, dpi=180)
    print(f"Saved figure to: {plot_path}")

    local_path = "/home/e.lomba/MC_Checkerboard/test_cases/mc_vs_lammps_corrected.png"
    plt.savefig(local_path, dpi=180)
    print(f"Saved figure to: {local_path}")

if __name__ == "__main__":
    main()
