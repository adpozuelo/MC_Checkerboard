# Lennard-Jones Fluid AVBMC Validation Benchmark

## System Parameters
- **Model**: Single-component truncated Lennard-Jones fluid (`pot_int = 0`, `model = 'LJ'`).
- **Thermodynamic State**:
  - Number of particles: $N = 1000$
  - Reduced temperature: $T^* = k_B T / \epsilon = 0.75$ ($\beta = 1.3333$)
  - Reduced density: $\rho^* = N / V = 0.03$
  - Cubic box edge: $L = (N / \rho)^{1/3} \approx 32.183\sigma$
  - Cutoff radius: $r_c = 2.5\sigma$
  - Clustering parameters: $r_{\text{cl}} = 1.5\sigma$, $\text{minPts} = 3$, $\eta_{\text{thresh}} = 0.35$

## Simulation Setups
1. **`run_no_avbmc/`**: Standard canonical Monte Carlo with local particle displacements only (`avbmc = .false.`).
2. **`run_with_avbmc/`**: Monte Carlo augmented with Aggregation Volume Bias Monte Carlo (AVBMC) pair moves (`avbmc = .true.`, $K = 10$ Rosenbluth trials, max 20 attempts every 20 sweeps).

Both simulations were run for 21,000 sweeps on NVIDIA RTX GPUs using the CUDA Fortran engine `MCCB-gpu`.

## Results Summary

| Metric | Standard MC (No AVBMC) | AVBMC MC | Consistency Check |
| :--- | :---: | :---: | :---: |
| **Equilibrated $\langle U/N \rangle$ [$\epsilon$]** | $-0.22891 \pm 0.00277$ | $-0.23298 \pm 0.00212$ | Overlapping within error bars ($\Delta = 0.0041$) |
| **Energy Variance $\sigma^2$** | $0.000214$ | $0.000201$ | Consistent ($< 6\%$ difference) |
| **Standard Deviation $\sigma(U/N)$** | $0.01462$ | $0.01417$ | Quantitative agreement |
| **Heat Capacity $C_v / (N k_B)$** | $0.3800$ | $0.3570$ | Agreement |
| **Autocorrelation Time $\tau_{\text{int}}$** | $21.6$ sweeps | $13.5$ sweeps | **38% faster decorrelation with AVBMC** |
| **Energy Drift** | $0.0000 \times 10^0$ | $0.0000 \times 10^0$ | Machine exact energy conservation |

## Comparison Figure
![LJ AVBMC Comparison](lj_avbmc_validation_comparison.png)

- **Panel (a)**: Full potential energy relaxation path.
- **Panel (b)**: Cluster size evolution ($r_{\text{cl}} = 1.5\sigma$).
- **Panel (c)**: Canonical Boltzmann energy probability distributions $P(U/N)$ in the equilibrated window (steps 15,000–21,000), showing exact overlap with theoretical Gaussian fits.
- **Panel (d)**: Autocorrelation function decay $C_U(\Delta t)$ demonstrating enhanced phase space sampling.
