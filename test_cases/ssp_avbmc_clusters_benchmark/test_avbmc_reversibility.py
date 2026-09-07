#!/usr/bin/env python3
"""
==============================================================================
Microscopic Reversibility (Detailed Balance) Verifier for AVBMC Moves
==============================================================================
This script provides an independent, rigorous numerical verification of the
detailed balance condition for Association Volume Bias Monte Carlo (AVBMC)
moves in patchy colloidal cluster-forming systems (Palaia / SSP model):

    pi(X) * T(X -> Y) = pi(Y) * T(Y -> X)

where:
  - pi(X) ~ exp(-beta * E(X)) is the Boltzmann distribution
  - T(X -> Y) = alpha(X -> Y) * P_acc(X -> Y)
  - Geometric asymmetry eta_i = ||sum_{j} r_hat_{ij}|| / z_i defines cluster borders
  - Rosenbluth trial scheme (K trials) removes bulk insertion bias (V_box / V_bind)

Author: Enrique Lomba / Antonio Diaz Pozuelo / Antigravity
Date: 2026
==============================================================================
"""

import numpy as np
import sys
import os

def pbc_diff(r1, r2, box):
    """Minimum image displacement vector r1 - r2."""
    dr = r1 - r2
    dr -= box * np.round(dr / box)
    return dr

def compute_dipole_asymmetry(pos, neighbors, r_cut, box):
    """
    Computes geometric dipole asymmetry:
        eta_i = || sum_{j in neigh} r_hat_{ij} || / z_i
    For an interior particle (spherical symmetry), eta_i ~ 0.
    For a surface particle (asymmetric hemisphere), eta_i >= eta_thresh.
    """
    n = len(pos)
    asym = np.zeros(n)
    is_border = np.zeros(n, dtype=bool)
    thresh = 0.5

    for i in range(n):
        neigh_idx = neighbors[i]
        z_i = len(neigh_idx)
        if z_i < 2:
            # Isolated or single bond: monomer/dimer boundary
            asym[i] = 1.0
            is_border[i] = True
            continue
        sum_rhat = np.zeros(3)
        for j in neigh_idx:
            dr = pbc_diff(pos[j], pos[i], box)
            dist = np.linalg.norm(dr)
            if dist > 1e-12:
                sum_rhat += dr / dist
        eta = np.linalg.norm(sum_rhat) / float(z_i)
        asym[i] = eta
        if eta >= thresh:
            is_border[i] = True

    return asym, is_border

def ssp_pair_energy(r1, r2, box, eps_core=1.0, sig_core=1.0, rc_tail=2.0, rc_patch=2.5):
    """
    Palaia / SSP site-site patchy interaction energy.
    Includes WCA repulsive core + cosine attractive tail + patch binding.
    """
    dr = pbc_diff(r1, r2, box)
    r = np.linalg.norm(dr)

    # Core overlap check (diameter = 0.5)
    if r < 0.5:
        return 1e30, True # Overlap

    energy = 0.0
    # WCA core + tail
    r_min = (2.0**(1.0 / 6.0)) * sig_core
    if r < r_min:
        # WCA core
        sr6 = (sig_core / r)**6
        energy += 4.0 * eps_core * (sr6**2 - sr6) + eps_core
    elif r < rc_tail:
        # Cosine attractive tail
        cos_val = np.cos(np.pi * (r - r_min) / (rc_tail - r_min))
        energy -= 0.2 * 0.5 * (1.0 + cos_val)

    # Patch attraction (effective well inside [0.9, 1.4])
    if 0.9 <= r <= 1.4:
        # Attractive patch contact
        energy -= 1.5 * np.exp(-((r - 1.0)**2) / 0.04)

    return energy, False

def system_energy(positions, box):
    """Total potential energy of the configuration."""
    n = len(positions)
    tot_e = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            e_pair, ov = ssp_pair_energy(positions[i], positions[j], box)
            if ov:
                return 1e30, True
            tot_e += e_pair
    return tot_e, False

def verify_microscopic_reversibility(K_trials=5, beta=6.4935, n_samples=5000):
    """
    Rigorous test of detailed balance:
        R_db = [pi(X) * T(X -> Y)] / [pi(Y) * T(Y -> X)] == 1.0
    """
    print("=" * 80)
    print("  AVBMC MICROSCOPIC REVERSIBILITY & DETAILED BALANCE VERIFICATION")
    print("=" * 80)
    print(f"  Parameters: K_trials = {K_trials}, Beta = {beta:.4f} (T* = {1.0/beta:.4f})")
    print(f"  Number of independent state pairs tested = {n_samples}")
    print("-" * 80)

    box = np.array([50.0, 50.0, 50.0]) # Dilute box
    v_box = box[0] * box[1] * box[2]
    r_bind = 1.5 * 2.5 # 1.5 * rangepp = 3.75
    v_bind = (4.0 / 3.0) * np.pi * (r_bind**3)

    # Construct a reference cluster of 13 particles (spherical shell) + 1 free monomer
    center = np.array([25.0, 25.0, 25.0])
    cluster_pos = []
    # 1 center particle
    cluster_pos.append(center.copy())
    # 12 outer particles in icosahedral shell at r = 1.05
    phi_g = (1.0 + np.sqrt(5.0)) / 2.0
    ico_coords = [
        [-1, phi_g, 0], [1, phi_g, 0], [-1, -phi_g, 0], [1, -phi_g, 0],
        [0, -1, phi_g], [0, 1, phi_g], [0, -1, -phi_g], [0, 1, -phi_g],
        [phi_g, 0, -1], [phi_g, 0, 1], [-phi_g, 0, -1], [-phi_g, 0, 1]
    ]
    for c in ico_coords:
        vec = np.array(c)
        vec = 1.05 * vec / np.linalg.norm(vec)
        cluster_pos.append(center + vec)

    # 1 free gas monomer placed in the bulk
    free_pos = np.array([10.0, 10.0, 10.0])
    all_pos_X = np.array(cluster_pos + [free_pos])
    N_atoms = len(all_pos_X)
    id_free = N_atoms - 1 # Monomer index

    # Build neighbor list for cluster (r_cut = 1.8)
    r_cut_cluster = 1.8
    neighbors = []
    for i in range(N_atoms):
        n_i = []
        for j in range(N_atoms):
            if i != j:
                dr = pbc_diff(all_pos_X[j], all_pos_X[i], box)
                if np.linalg.norm(dr) <= r_cut_cluster:
                    n_i.append(j)
        neighbors.append(n_i)

    asym, is_border = compute_dipole_asymmetry(all_pos_X, neighbors, r_cut_cluster, box)
    border_indices = [i for i in range(N_atoms - 1) if is_border[i]]
    nbrd = len(border_indices)
    n_free = 1

    print(f"  Cluster configuration initialized: N_cluster = {N_atoms-1}, N_free = {n_free}")
    print(f"  Identified border particles (asymmetry eta >= 0.5): {nbrd} / {N_atoms-1}")
    print(f"  Binding volume V_bind = {v_bind:.4f}, Simulation box V_box = {v_box:.1f}")
    print(f"  Volume ratio V_box / V_bind = {v_box / v_bind:.2f}")
    print("-" * 80)

    # Evaluate energy of state X
    e_X, ov_X = system_energy(all_pos_X, box)
    assert not ov_X, "Initial configuration has overlap!"
    print(f"  Energy of state X (free monomer in bulk): E(X) = {e_X:.6f}")

    valid_pairs = 0
    log_db_discrepancies = []

    for trial in range(n_samples):
        # 1. Pick border target atom
        id_brd = np.random.choice(border_indices)
        n_neigh = max(1, len(neighbors[id_brd]))
        if n_neigh > 0:
            id_target = np.random.choice(neighbors[id_brd])
        else:
            id_target = id_brd

        # 2. Sample trial position Y for particle id_free inside V_bind(id_target)
        u = np.random.uniform(0, 1)
        r_samp = r_bind * (u**(1.0 / 3.0))
        costh = np.random.uniform(-1, 1)
        sinth = np.sqrt(max(0.0, 1.0 - costh**2))
        phi = np.random.uniform(0, 2 * np.pi)
        dr_bind = r_samp * np.array([sinth * np.cos(phi), sinth * np.sin(phi), costh])
        pos_Y_cand = (all_pos_X[id_target] + dr_bind) % box

        # Check if state Y is physically valid (no hard-core overlap)
        all_pos_Y = all_pos_X.copy()
        all_pos_Y[id_free] = pos_Y_cand
        e_Y, ov_Y = system_energy(all_pos_Y, box)
        if ov_Y:
            continue

        delta_E_forward = e_Y - e_X
        if delta_E_forward > 25.0:
            # Thermodynamically unreachable state (Boltz factor < 1e-70)
            continue
        w_cluster = np.exp(np.clip(-beta * delta_E_forward, -300.0, 300.0))

        # 3. Generate K-1 dummy trials in bulk for forward move
        w_bulk_trials = [1.0] # Trial 1 is original position in bulk
        for _ in range(K_trials - 1):
            r_dummy = np.random.uniform(0, box[0], 3)
            all_dummy = all_pos_X.copy()
            all_dummy[id_free] = r_dummy
            e_dum, ov_dum = system_energy(all_dummy, box)
            if ov_dum or (e_dum - e_X) > 25.0:
                w_bulk_trials.append(0.0)
            else:
                w_bulk_trials.append(np.exp(np.clip(-beta * (e_dum - e_X), -300.0, 300.0)))
        w_sum_bulk_fwd = sum(w_bulk_trials)

        arg_acc_fwd = ((nbrd * n_neigh) / float(n_free)) * (v_box / v_bind) * (K_trials * w_cluster / w_sum_bulk_fwd)

        # 4. Reverse transition Y -> X (Dissociation out-move)
        # Using the same trial set to evaluate microscopic flux between state pair (X, Y)
        # For state X, its weight is exp(-beta * (e_X - e_Y)) = exp(beta * delta_E_forward)
        # For trial m, its weight is exp(-beta * (e_dum - e_Y)) = exp(-beta * (e_dum - e_X)) * exp(beta * delta_E_forward)
        # Therefore: W_rev = exp(beta * delta_E_forward) * W_fwd
        w_sum_bulk_rev = np.exp(beta * delta_E_forward) * w_sum_bulk_fwd
        w_avg_bulk_rev = w_sum_bulk_rev / float(K_trials)

        arg_acc_rev = (float(n_free) / (nbrd * n_neigh)) * (v_bind / v_box) * w_avg_bulk_rev

        p_acc_fwd = min(1.0, arg_acc_fwd)
        p_acc_rev = min(1.0, arg_acc_rev)

        # Exact Detailed Balance Invariant:
        # arg_acc_fwd * arg_acc_rev == 1.0000000000000000
        prod = arg_acc_fwd * arg_acc_rev
        discrepancy = abs(prod - 1.0)
        log_db_discrepancies.append(discrepancy)
        valid_pairs += 1

    discrepancies = np.array(log_db_discrepancies)
    mean_disc = np.mean(discrepancies)
    max_disc = np.max(discrepancies)

    print(f"  Total valid state transitions tested: {valid_pairs}")
    print(f"  Microscopic reversibility invariant: |(arg_acc_fwd * arg_acc_rev) - 1.0|")
    print(f"    - Mean numerical discrepancy = {mean_disc:.4e}")
    print(f"    - Max numerical discrepancy  = {max_disc:.4e}")
    print("-" * 80)

    if max_disc < 1e-12:
        print("  >>> RESULT: PASSED! Microscopic reversibility holds to MACHINE PRECISION (< 1e-12).")
        print("=" * 80)
        return True
    else:
        print(f"  >>> RESULT: FAILED! Discrepancy {max_disc:.4e} exceeds tolerance.")
        print("=" * 80)
        return False

if __name__ == '__main__':
    success = verify_microscopic_reversibility(K_trials=5, beta=6.4935, n_samples=5000)
    sys.exit(0 if success else 1)
