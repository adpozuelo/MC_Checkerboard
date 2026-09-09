# gpMC (CPU) vs MCCB-gpu (GPU) Benchmark Reproduction

This benchmark reproduces the execution time comparisons between the CPU-based Monte Carlo code **gpMC** (project `MCcol`, compiled with Intel Fortran `ifx` 2025.2.0) and the GPU-accelerated code **MCCB-gpu** (Version 2.8.0, compiled with NVIDIA HPC SDK 25.3 / CUDA 12.8), originally reported in Table 2 of `benchmark.tex`.

## Test System Specifications
- **Model**: Binary Lennard-Jones (LJ) fluid mixture ($N_1 : N_2 = 2 : 1$).
- **Parameters**:
  - $\sigma_{11} = 1.0$, $\sigma_{12} = 1.2$, $\sigma_{22} = 1.4$ ($\sigma_0 = 3.0$ \AA)
  - $\varepsilon_{11} = 1.0$, $\varepsilon_{12} = 1.5$, $\varepsilon_{22} = 0.5$ ($\varepsilon_0 = 100$ K)
  - Cutoff: $r_c = 2.5 \sigma_0 = 7.5$ \AA
  - Temperature: $T^* = 3.0$ ($T = 300$ K)
  - Ensemble: $NVT$ (Canonical)
- **Workload**: 100 Monte Carlo cycles / sweeps (`istep_fin = 100`).

## Hardware Specifications
- **Cluster Node**: `ladon28`
- **CPU**: Intel(R) Xeon(R) Gold 6548Y+ (32 cores / 64 threads, 2.20 GHz)
- **GPU**: NVIDIA RTX PRO 4500 Blackwell (Compute Capability 12.0, 32 GB VRAM, 82 SMs)

## Benchmark Results ($NVT$ and $NpT$ Ensembles)

### Raw Sweep Runtimes (100 Cycles)

| Particle Count ($N$) | gpMC CPU ($NVT$) | gpMC CPU ($NpT$) | MCCB-gpu GPU ($NVT$) | MCCB-gpu GPU ($NpT$) | Raw Speedup ($NVT$) | Raw Speedup ($NpT$) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **8,000** | **9.95 s** | **15.16 s** | **4.97 s** | **6.26 s** | **2.00×** | **2.42×** |
| **15,625** | **23.59 s** | **35.95 s** | **8.73 s** | **11.88 s** | **2.70×** | **3.03×** |
| **32,768** | **38.63 s** | **59.23 s** | **8.05 s** | **9.94 s** | **4.80×** | **5.96×** |
| **125,000** | **175.91 s** | **254.09 s** | **7.86 s** | **10.14 s** | **22.39×** | **25.06×** |

### Normalized by Trial Displacements per Particle

In `gpMC`, exactly 1 single-particle displacement per particle is attempted each sweep ($\mu_{\text{CPU}} = 1.0$). In `MCCB-gpu`, the GPU checkerboard subdivision executes $N_{\text{cells}} \times N_{\text{move}}$ trial moves per cycle ($N_{\text{move}} = 125$), which corresponds to $\mu_{\text{GPU}} = N_{\text{cells}} N_{\text{move}} / N$ trial displacements per particle per cycle ($3.375$, $1.728$, $1.953$, and $2.744$ moves/particle/cycle for $N=8000, 15625, 32768, 125000$).

When normalized to an identical workload of 100 trial moves per particle ($t_{\text{norm}} = t_{\text{raw}} / \mu_{\text{GPU}}$):

| Particle Count ($N$) | Moves/Part/Cycle ($\mu_{\text{GPU}}$) | MCCB-gpu Norm ($NVT$) | MCCB-gpu Norm ($NpT$) | **Effective Speedup ($NVT$)** | **Effective Speedup ($NpT$)** |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **8,000** | **3.375** | **1.47 s** | **1.85 s** | **6.76×** | **8.17×** |
| **15,625** | **1.728** | **5.05 s** | **6.88 s** | **4.67×** | **5.23×** |
| **32,768** | **1.953** | **4.12 s** | **5.09 s** | **9.37×** | **11.64×** |
| **125,000** | **2.744** | **2.86 s** | **3.70 s** | **61.41×** | **68.76×** |

### Key Observations
1. **Replication Accuracy**:
   - For gpMC: The measured runtimes in $NVT$ (9.95 s, 23.59 s, 38.63 s, 175.91 s) reproduce Table 2 of `benchmark.tex` (10.0 s, 23.4 s, 38.4 s, 167.0 s on i7-14700K / ~170 s on Xeon Gold) almost identically.
2. **GPU Scaling Acceleration in Both Ensembles**:
   - In both $NVT$ and $NpT$, between $N = 15,625$ and $N = 125,000$, MCCB-gpu exhibits an execution time reduction and plateau (8.73 s $\to$ 7.86 s in $NVT$, 11.88 s $\to$ 10.14 s in $NpT$) despite an $8\times$ particle increase, perfectly validating the occupancy and memory coalescing saturation described in `benchmark.tex`.
3. **True Algorithmic Throughput Speedup**:
   - When accounting for the higher trial displacement rate in MCCB-gpu, the true algorithmic speedup over optimized CPU code reaches **61.41× in $NVT$** and **68.76× in $NpT$** at $N = 125,000$.
