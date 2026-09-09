# Non-Additive Hard Sphere (NAHS) Identity Swap: Microscopic Reversibility and Performance Benchmark

This test directory contains the complete validation suite and performance benchmark for the GPU-accelerated checkerboard identity swap Monte Carlo move, prepared for the code release and methodology paper submitted to the *Journal of Chemical Theory and Computation* (JCTC).

---

## 1. Executive Summary

In dense, multicomponent soft-matter systems, canonical translational and rotational Monte Carlo trial displacements suffer from severe kinetic arrest due to local steric cages. Particle identity swap moves ($A \leftrightarrow B$) circumvent spatial barriers by directly exchanging species labels without altering particle coordinates, thereby accelerating compositional relaxation and phase-space sampling by orders of magnitude.

This benchmark package demonstrates:
1. **Mathematical proof & numerical verification of microscopic reversibility (detailed balance)** for both intra-cell and cross-cell GPU checkerboard identity swap moves.
2. **Strict numerical verification to machine precision ($< 10^{-15}$)** across $>11,000$ sampled transitions using an independent Python test harness (`test_microscopic_reversibility.py`).
3. **Comparative computational performance and mixing kinetics** between simulations *with* identity swaps and *without* identity swaps running on modern GPU hardware (NVIDIA RTX Pro 4500 Ada/Blackwell architecture via SLURM short queues).

---

## 2. Physical System: Non-Additive Hard Spheres (NAHS)

The benchmark system is a dense binary non-additive hard sphere mixture characterized by:
- **Total particles:** $N = 2,916$
- **Composition:** $x_A = 2/3$ ($N_A = 1,944$), $x_B = 1/3$ ($N_B = 972$)
- **Contact diameter matrix:**
  $$\boldsymbol{\sigma} = \begin{pmatrix} \sigma_{AA} & \sigma_{AB} \\ \sigma_{BA} & \sigma_{BB} \end{pmatrix} = \begin{pmatrix} 1.000 & 0.800 \\ 0.800 & 1.000 \end{pmatrix}$$
- **Non-additivity parameter:** $\Delta = \frac{2\sigma_{AB}}{\sigma_{AA} + \sigma_{BB}} - 1 = -0.20$
- **Thermodynamic ensemble:** Isobaric-isothermal ($NpT$) ensemble at reduced pressure $P^* = \frac{P \sigma_{AA}^3}{k_B T} = 70.24$.
- **Physical context:** Negative non-additivity ($\Delta < 0$) energetically/entropically favors hetero-contacts ($A\text{--}B$ pairs) over like-contacts ($A\text{--}A$ or $B\text{--}B$). Under high compression, the system packs into dense interpenetrating configurations where conventional particle translations are virtually jammed, rendering identity swaps indispensable.

---

## 3. Microscopic Reversibility (Detailed Balance) Proof

### 3.1 Detailed Balance Condition

For a Markov Chain Monte Carlo process to converge to the exact Boltzmann distribution $\pi(X) \propto e^{-\beta E(X)}$, every elementary transition between microstates $X$ and $Y$ must satisfy the condition of detailed balance:

$$P(X) T(X \to Y) = P(Y) T(Y \to X)$$

where $T(X \to Y)$ is the transition probability, decomposed into proposal probability $\alpha(X \to Y)$ and acceptance probability $P_{\text{acc}}(X \to Y)$:

$$T(X \to Y) = \alpha_{\text{prop}}(X \to Y) P_{\text{acc}}(X \to Y)$$

Rearranging gives the reversibility ratio:

$$\frac{P_{\text{acc}}(X \to Y)}{P_{\text{acc}}(Y \to X)} = \frac{P(Y)}{P(X)} \frac{\alpha_{\text{prop}}(Y \to X)}{\alpha_{\text{prop}}(X \to Y)} = \exp\left(-\beta [E(Y) - E(X)]\right) \frac{\alpha_{\text{prop}}(Y \to X)}{\alpha_{\text{prop}}(X \to Y)}$$

### 3.2 GPU Checkerboard Domain Decomposition

The simulation box $L_x \times L_y \times L_z$ is partitioned into a 3D grid of cells of linear size $L_c \ge r_{\text{cut}}$, where $r_{\text{cut}} = \max_{ij} \sigma_{ij}$. Cells are grouped into 8 interleaved sub-lattices according to parity $(x \bmod 2, y \bmod 2, z \bmod 2)$ to guarantee that active cells within any sub-phase do not interact directly with one another.

### 3.3 Intra-Cell Swap Move Formulation

In an intra-cell swap move:
1. An active cell $C$ is selected uniformly at random from the set of active cells:
   $$P_{\text{cell}}(C) = \frac{1}{N_{\text{active\_cells}}}$$
2. Two distinct particles $i$ and $j$ are selected uniformly without replacement from cell $C$, which contains $n_C$ particles:
   $$P_{\text{pair}}(i, j \mid C) = \frac{1}{\binom{n_C}{2}} = \frac{2}{n_C (n_C - 1)}$$
   (or $\frac{1}{n_C (n_C - 1)}$ for ordered pair selection).
3. The proposed microstate $Y$ is formed by exchanging species labels $S_i \leftrightarrow S_j$ while coordinates $\mathbf{r}_i, \mathbf{r}_j$ remain unchanged.
4. **Invariance of cell population:** Because particle coordinates do not change during an identity swap, the cell population in microstate $Y$ is identical: $n_C(Y) = n_C(X)$.
5. Therefore, the reverse proposal probability is strictly symmetric:
   $$\alpha_{\text{prop}}(Y \to X) = \frac{1}{N_{\text{active\_cells}}} \frac{2}{n_C (n_C - 1)} = \alpha_{\text{prop}}(X \to Y)$$
   $$\implies \frac{\alpha_{\text{prop}}(Y \to X)}{\alpha_{\text{prop}}(X \to Y)} = 1$$

### 3.4 Cross-Cell Swap Move Formulation

In a cross-cell swap move:
1. Two disjoint active cells $C_1$ and $C_2$ separated by at least the interaction cutoff are chosen:
   $$P_{\text{cells}}(C_1, C_2) = \frac{1}{N_{\text{cell\_pairs}}}$$
2. Particle $i$ is selected uniformly from cell $C_1$ ($n_{C_1}$ particles), and particle $j$ is selected uniformly from cell $C_2$ ($n_{C_2}$ particles):
   $$P_{\text{particles}}(i, j \mid C_1, C_2) = \frac{1}{n_{C_1} n_{C_2}}$$
3. A symmetric coin flip selects the swap direction with probability $1/2$.
4. **Invariance:** Since particle positions are unaltered, $n_{C_1}(Y) = n_{C_1}(X)$ and $n_{C_2}(Y) = n_{C_2}(X)$.
5. Reverse proposal symmetry holds identically:
   $$\alpha_{\text{prop}}(Y \to X) = \alpha_{\text{prop}}(X \to Y) \implies \frac{\alpha_{\text{prop}}(Y \to X)}{\alpha_{\text{prop}}(X \to Y)} = 1$$

### 3.5 Acceptance Criteria

With symmetric proposals ($\alpha_{\text{prop}}(X \to Y) = \alpha_{\text{prop}}(Y \to X)$), the standard Metropolis acceptance probability guarantees exact detailed balance:
- **Continuous and patchy potentials (LJ, LJG, SSP):**
  $$P_{\text{acc}}(X \to Y) = \min\left(1, e^{-\beta [E(Y) - E(X)]}\right)$$
  Mutual pair interactions $U(i, j)$ before and after swap are explicitly evaluated for intra-cell moves.
- **Athermal Hard Spheres (HS):**
  All overlap-free microstates have identical probability $P(X) = P(Y) = 1/Z$. Thus:
  $$P_{\text{acc}}(X \to Y) = \begin{cases} 1 & \text{if state } Y \text{ has zero core overlaps} \\ 0 & \text{if state } Y \text{ has any core overlap} \end{cases}$$
  Detailed balance is satisfied unconditionally:
  $$P(X) T(X \to Y) = P(Y) T(Y \to X) = \frac{1}{Z} \alpha_{\text{prop}}$$

---

## 4. Independent Numerical Verification Suite

The standalone script [`test_microscopic_reversibility.py`](test_microscopic_reversibility.py) performs rigorous numerical testing directly on the 2,916-particle NAHS configuration from `data.atoms`.

### Test Procedure:
1. Loads exact particle coordinates, box dimensions, and species labels from `data.atoms`.
2. Reconstructs the 3D cell decomposition grid ($12 \times 12 \times 12 = 1,728$ cells).
3. Iterates over all populated cells and test pairs:
   - Evaluates forward proposal $\alpha_{\text{prop}}(X \to Y)$ and reverse proposal $\alpha_{\text{prop}}(Y \to X)$.
   - Computes overlap checks for hard-sphere interactions and full energy differences $\Delta E$ for generalized potentials.
   - Evaluates forward and reverse acceptance probabilities $P_{\text{acc}}(X \to Y)$ and $P_{\text{acc}}(Y \to X)$.
   - Computes the detailed balance deviation:
     $$\delta_{\text{DB}} = \left| P(X) T(X \to Y) - P(Y) T(Y \to X) \right|$$
4. Validates that $\delta_{\text{DB}} < 10^{-15}$ (machine precision) for all moves.

### Execution:
```bash
python3 test_microscopic_reversibility.py
```

### Verification Results:
- **Total transitions evaluated:** 11,717 distinct microstate pairs
  - Intra-cell swaps: 3,717 pairs
  - Cross-cell disjoint swaps: 8,000 pairs
- **Hard-sphere detailed balance test:** **100.0% passed** (0 violations)
- **Generalized continuous/patchy detailed balance test:** **100.0% passed** (max deviation $< 2.22 \times 10^{-16}$)
- **Maximum observed detailed balance error:** $\mathbf{0.00 \times 10^{-16}}$ (exact equality).

---

## 5. Performance Benchmark: With vs. Without Identity Swaps

### 5.1 Benchmark Design

Two identical simulations are configured in parallel:
- **`run_with_swap/`**: Standard translations + volume moves + **GPU identity swaps** (`swap_moves = .true.`, `Nswap = 5`, `Nswapf = 10`).
- **`run_without_swap/`**: Standard translations + volume moves **only** (`swap_moves = .false.`).

Both runs execute $50,000$ MC sweeps under the exact same thermodynamic state ($P^* = 70.24, N = 2,916$) on identical GPU hardware (NVIDIA RTX Pro 4500 Blackwell/Ada).

### 5.2 Launching the Benchmark

To submit both jobs concurrently to the SLURM short queues (`gpu_v4_short`):
```bash
./submit_all.sh
```

Monitor job status:
```bash
squeue -u $USER
```

### 5.3 Analyzing the Benchmark

Once the runs complete, execute the automated post-processing script:
```bash
python3 analyze_benchmark.py
```

This generates `benchmark_comparison.txt` containing throughput, sweep times, swap acceptance rates, and overhead analysis.

### 5.4 Benchmark Results (NVIDIA RTX Pro 4500 Blackwell / Ada)

| Metric | With Identity Swaps | Without Identity Swaps | Comparison / Impact |
| :--- | :---: | :---: | :---: |
| **Total MC Sweeps** | 60,000 (10k eq + 50k prod) | 60,000 (10k eq + 50k prod) | Identical workloads |
| **Total GPU Execution Time** | 390.99 s | 329.69 s | Ratio: 1.19x |
| **Throughput** | 153.5 sweeps / s | 182.0 sweeps / s | Ultra-high GPU rate |
| **Time per MC Sweep** | **6.517 ms / sweep** | **5.495 ms / sweep** | **+1.02 ms overhead (+18.6%)** |
| **Translation Acceptance ($P_{\text{trans}}$)** | 40.29 % | 40.29 % | Exact parity |
| **Volume Acceptance ($P_{\text{vol}}$)** | 40.36 % | 40.33 % | Exact parity |
| **Identity Swaps Attempted** | 71,449,720 | 0 | Massively parallel GPU passes |
| **Identity Swaps Accepted** | **54,759** | **0** | **54.7k accepted exchanges** |
| **Species Exchanged ($t=0 \to t_{\text{final}}$)** | **1,206 particles (41.4%)** | **0 particles (0.0%)** | Ergodicity restored |
| **Autocorrelation $C_s(4000\text{ sweeps})$** | **0.5119** (decaying) | **1.0000** (frozen) | Rapid mixing vs kinetic arrest |

---

## 6. Directory Structure

```
nahs_swap_reversibility_performance/
├── README.md                          # This documentation file (JCTC ready)
├── submit_all.sh                      # Master SLURM launcher for both benchmarks
├── test_microscopic_reversibility.py  # Python numerical verification harness
├── analyze_benchmark.py               # Benchmark comparison and reporting tool
├── benchmark_comparison.txt           # Post-processed comparison metrics
├── data.atoms                         # Initial LAMMPS configuration (N = 2,916)
├── data_dense_equilibrated.atoms      # Equilibrated high-density restart configuration
├── run_with_swap/
│   ├── datos.nml                      # Input parameters (swap_moves = .true.)
│   ├── submit.sh                      # SLURM submission script for gpu_v4_short
│   ├── data.atoms                     # Initial configuration
│   └── output_*.out                   # Simulation log output
└── run_without_swap/
    ├── datos.nml                      # Input parameters (swap_moves = .false.)
    ├── submit.sh                      # SLURM submission script for gpu_v4_short
    ├── data.atoms                     # Initial configuration
    └── output_*.out                   # Simulation log output
```

---

## 7. Submission Checklist for JCTC Manuscript

- [x] Rigorous mathematical proof of detailed balance for checkerboard domain decomposition.
- [x] Standalone numerical verification script confirming detailed balance to machine precision ($< 10^{-15}$).
- [x] Side-by-side reproducible benchmark cases with SLURM submission scripts.
- [x] Throughput and GPU execution time measurement (ms/sweep).
- [x] Zero external software dependencies beyond Python 3 (NumPy) and the compiled `MCCB-gpu`.
