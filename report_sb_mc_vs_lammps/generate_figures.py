#!/usr/bin/env python3
"""
Generate publication-quality figures for the LaTeX comparison report.
Outputs both PDF and PNG formats.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'lines.linewidth': 1.8,
    'axes.linewidth': 1.0,
})

def main():
    mc_dir_54k = "/home/e.lomba/MC_Checkerboard/test_cases/sb_avbmc_50_50_3x3x3"
    mc_prev_54k = "/home/e.lomba/MC_Checkerboard/test_cases/sb_avbmc_50_50_3x3x3/prev_unlike_only"
    lmp_dir = "/home/e.lomba/sanchez_burgos_20260915/sanchez_burgos_lammps/runs/paper_a/seed_712345"
    out_dir = "/home/e.lomba/MC_Checkerboard/report_sb_mc_vs_lammps"

    os.makedirs(out_dir, exist_ok=True)

    # -----------------------------------------------------------------
    # Figure 1: Energy Evolution & Stabilization
    # -----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.8))

    # Previous unlike-only
    prev_file = os.path.join(mc_prev_54k, "run-data.dat")
    if os.path.exists(prev_file):
        rows = []
        with open(prev_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 8:
                    try:
                        rows.append([float(x) for x in parts])
                    except ValueError:
                        continue
        if rows:
            d_prev = np.array(rows)
            stride = max(1, len(d_prev) // 600)
            ax.plot(d_prev[::stride, 0], d_prev[::stride, 7], label=r"MC Unlike-Only CSW ($N=54,000$)", color="#7f7f7f", ls=":", lw=1.6)

    # Corrected dual CSW
    cur_file = os.path.join(mc_dir_54k, "run-data.dat")
    if os.path.exists(cur_file):
        rows = []
        with open(cur_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 8:
                    try:
                        rows.append([float(x) for x in parts])
                    except ValueError:
                        continue
        if rows:
            d_cur = np.array(rows)
            stride = max(1, len(d_cur) // 800)
            ax.plot(d_cur[::stride, 0], d_cur[::stride, 7], label=r"MC Dual CSW ($N=54,000$, Corrected)", color="#1f77b4", lw=2.2)

    ax.axhline(-15.27, color="#d62728", ls="--", lw=1.8, label=r"LAMMPS MD Continuation ($\sim -15.27\,k_BT$)")
    ax.axhline(-15.12, color="#ff7f0e", ls="-.", lw=1.6, label=r"LAMMPS MD Slab Audit ($\sim -15.12\,k_BT$)")
    ax.axhline(-15.60, color="#1f77b4", ls="--", alpha=0.6, label=r"MC Dual CSW Plateau ($\sim -15.60\,k_BT$)")

    ax.set_xlabel("Monte Carlo Sweep", fontsize=12)
    ax.set_ylabel(r"Potential Energy per Molecule $\langle U/N \rangle$ [$k_B T$]", fontsize=12)
    ax.set_title("Potential Energy Evolution & Stabilization", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="center right", framealpha=0.95)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 150000)
    ax.set_ylim(-16.4, -13.6)
    plt.tight_layout()

    fig.savefig(os.path.join(out_dir, "fig1_energy.pdf"))
    fig.savefig(os.path.join(out_dir, "fig1_energy.png"), dpi=250)
    plt.close(fig)
    print("Saved Figure 1 (Energy Evolution).")

    # -----------------------------------------------------------------
    # Figure 2: Radial Distribution Functions g(r)
    # -----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.8))

    g_mc_file = os.path.join(mc_dir_54k, "gmixsim.dat")
    if os.path.exists(g_mc_file):
        g_mc = np.loadtxt(g_mc_file)
        r_mc = g_mc[:, 0]
        mask = (r_mc >= 0.85) & (r_mc <= 2.5)
        ax.plot(r_mc[mask], g_mc[mask, 4], label=r"MC: $g_{12}(r)$ (Scaffold--Surfactant)", color="#1f77b4", lw=2.2)
        ax.plot(r_mc[mask], g_mc[mask, 3], label=r"MC: $g_{11}(r)$ (Scaffold--Scaffold)", color="#1f77b4", ls="--", lw=1.8)
        ax.plot(r_mc[mask], g_mc[mask, 5], label=r"MC: $g_{22}(r)$ (Surfactant--Surfactant)", color="#1f77b4", ls=":", lw=1.6)

    g_lmp_file = os.path.join(lmp_dir, "gmixsim.dat")
    if os.path.exists(g_lmp_file):
        g_lmp = np.loadtxt(g_lmp_file)
        r_lmp = g_lmp[:, 0]
        mask_lmp = (r_lmp >= 0.85) & (r_lmp <= 2.5)
        ax.plot(r_lmp[mask_lmp], g_lmp[mask_lmp, 4], label=r"LAMMPS: $g_{12}(r)$", color="#d62728", lw=1.6, alpha=0.85)
        ax.plot(r_lmp[mask_lmp], g_lmp[mask_lmp, 3], label=r"LAMMPS: $g_{11}(r)$", color="#ff7f0e", ls="--", lw=1.6, alpha=0.85)
        ax.plot(r_lmp[mask_lmp], g_lmp[mask_lmp, 5], label=r"LAMMPS: $g_{22}(r)$", color="#2ca02c", ls=":", lw=1.4, alpha=0.85)

    ax.set_xlabel(r"Separation $r / \sigma$", fontsize=12)
    ax.set_ylabel(r"Radial Distribution Function $g_{ij}(r)$", fontsize=12)
    ax.set_title("Site-Center Pair Correlations $g_{ij}(r)$", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right", framealpha=0.95)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.85, 2.5)
    ax.set_ylim(0, 52)
    plt.tight_layout()

    fig.savefig(os.path.join(out_dir, "fig2_rdf.pdf"))
    fig.savefig(os.path.join(out_dir, "fig2_rdf.png"), dpi=250)
    plt.close(fig)
    print("Saved Figure 2 (RDF).")

    # -----------------------------------------------------------------
    # Figure 3: Structure Factor S(Q)
    # -----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.8))

    sq_mc_file = os.path.join(mc_dir_54k, "sq.dat")
    if os.path.exists(sq_mc_file):
        sq_mc = np.loadtxt(sq_mc_file)
        q_mc = sq_mc[:, 0]
        s_nn_mc = sq_mc[:, 1]
        mask_q = (q_mc >= 0.2) & (q_mc <= 8.0)
        ax.semilogy(q_mc[mask_q], s_nn_mc[mask_q], "o-", ms=3.5, label=r"MC Dual CSW $N=54,000$ ($L=73.5\,\sigma$)", color="#1f77b4")

    sq_lmp_file = os.path.join(lmp_dir, "sq.dat")
    if os.path.exists(sq_lmp_file):
        sq_lmp = np.loadtxt(sq_lmp_file)
        q_lmp = sq_lmp[:, 0]
        s_nn_lmp = sq_lmp[:, 1]
        mask_ql = (q_lmp >= 0.2) & (q_lmp <= 8.0)
        ax.semilogy(q_lmp[mask_ql], s_nn_lmp[mask_ql], "s--", ms=3.5, label=r"LAMMPS MD $N=1,820$ Slab ($L_x=43.7\,\sigma$)", color="#d62728", alpha=0.85)

    ax.set_xlabel(r"Wavevector $Q \sigma$", fontsize=12)
    ax.set_ylabel(r"Total Structure Factor $S_{NN}(Q)$ [log scale]", fontsize=12)
    ax.set_title(r"Structure Factor $S_{NN}(Q)$: Phase Separation Divergence", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9.5, loc="upper right", framealpha=0.95)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(0.2, 8.0)
    ax.set_ylim(0.1, 2e4)
    plt.tight_layout()

    fig.savefig(os.path.join(out_dir, "fig3_sq.pdf"))
    fig.savefig(os.path.join(out_dir, "fig3_sq.png"), dpi=250)
    plt.close(fig)
    print("Saved Figure 3 (Structure Factor).")

    # -----------------------------------------------------------------
    # Figure 4: Cluster Fraction & Composition
    # -----------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Cluster percentage
    cl_mc_file = os.path.join(mc_dir_54k, "clusevol_mc.dat")
    if os.path.exists(cl_mc_file):
        cl_mc = np.loadtxt(cl_mc_file)
        ax1.plot(cl_mc[:, 0], cl_mc[:, 4], label=r"MC Dual CSW ($N=54,000$)", color="#1f77b4", lw=2.0)

    cl_lmp_file = os.path.join(lmp_dir, "clusevol.dat")
    if os.path.exists(cl_lmp_file):
        cl_lmp = np.loadtxt(cl_lmp_file)
        conf_lmp = cl_lmp[:, 0]
        pct_lmp = cl_lmp[:, 2]
        stride = max(1, len(conf_lmp) // 200)
        ax1.plot(conf_lmp[::stride] * 5000, pct_lmp[::stride], label=r"LAMMPS MD ($N=1,820$)", color="#d62728", ls="--", lw=1.8)

    ax1.axhline(99.68, color="#1f77b4", ls=":", alpha=0.8, label=r"MC Plateau ($99.68\%$)")
    ax1.axhline(94.6, color="#d62728", ls=":", alpha=0.8, label=r"LAMMPS Plateau ($\sim 94.6\%$)")
    ax1.set_xlabel("Simulation Step / Scale", fontsize=11)
    ax1.set_ylabel(r"Condensed Molecule Fraction [\%]", fontsize=11)
    ax1.set_title("Condensed Fraction Evolution", fontsize=12, fontweight="bold")
    ax1.legend(fontsize=8.5, loc="lower right", framealpha=0.95)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(88, 101)

    # Scaffold fraction x1 in cluster
    trj_cl_file = os.path.join(mc_dir_54k, "clusevol.dat")
    if os.path.exists(trj_cl_file):
        d_trj_cl = np.loadtxt(trj_cl_file)
        # cols: conf, no. clusters, % particles, x_av(1)
        ax2.plot(d_trj_cl[:, 0], d_trj_cl[:, 3], label=r"MC Dual CSW ($N=54,000$)", color="#1f77b4", lw=2.0)

    if os.path.exists(cl_lmp_file):
        ax2.plot(conf_lmp[::stride], cl_lmp[::stride, 3], label=r"LAMMPS MD (DBSCAN)", color="#d62728", ls="--", lw=1.8)

    ax2.axhline(0.585, color="black", ls="-.", lw=1.5, label=r"Published Dense Target ($0.585 \pm 0.005$)")
    ax2.axhline(0.485, color="#1f77b4", ls=":", lw=1.5, label=r"MC Cluster Plateau ($\approx 0.485$)")

    ax2.set_xlabel("Analyzed Production Frame", fontsize=11)
    ax2.set_ylabel(r"Scaffold Mole Fraction in Cluster $x_1$", fontsize=11)
    ax2.set_title("Dense Phase Scaffold Fraction", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=8.5, loc="lower right", framealpha=0.95)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0.35, 0.65)

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig4_cluster.pdf"))
    fig.savefig(os.path.join(out_dir, "fig4_cluster.png"), dpi=250)
    plt.close(fig)
    print("Saved Figure 4 (Cluster Properties).")

if __name__ == "__main__":
    main()
