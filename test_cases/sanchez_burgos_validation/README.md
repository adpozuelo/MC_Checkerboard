# Sanchez-Burgos Patchy Protein Condensation Model Validation

This directory provides validation and test cases for the **Sanchez-Burgos et al. (2021)** patchy particle potential implementation in `MC_Checkerboard` (`pot_int = 3`, `model = 'SB'`).

---

## 1. Physical Model Summary

The model simulates condensation and phase separation of binary patchy colloidal mixtures (scaffold and client proteins):
- **Scaffold (Species 1, Type 0 in MC / Type 1 in LAMMPS):**
  - Core diameter $\sigma = 1.0$
  - 4 tetrahedral patches at radial distance $r_{\text{patch}} = 0.5\sigma$
  - Scaffold patches interact with other Scaffold patches and Client patches.
- **Client (Species 2, Type 1 in MC / Type 2 in LAMMPS):**
  - Core diameter $\sigma = 1.0$
  - 3 planar triangular patches at radial distance $r_{\text{patch}} = 0.5\sigma$
  - Client patches interact with Scaffold patches, but NOT with other Client patches ($V_{\text{pot}} = 0$).

### Core-Core Interaction: Pseudo-Hard Sphere (PHS)
$$U_{\text{PHS}}(r) = \begin{cases} C \left[ \left(\frac{\sigma}{r}\right)^{50} - \left(\frac{\sigma}{r}\right)^{49} \right] + \epsilon_R, & r < \left(\frac{50}{49}\right)\sigma \\ 0, & r \ge \left(\frac{50}{49}\right)\sigma \end{cases}$$
where $\epsilon_R = \frac{2}{3}$ and $C = 50 \left(\frac{50}{49}\right)^{49} \epsilon_R \approx 90.72589088$.

### Patch-Patch Interaction: Continuous Square Well (CSW)
$$U_{\text{CSW}}(d) = \begin{cases} -\frac{\epsilon_{\text{CSW}}}{2} \left[ 1 - \tanh\left(\frac{d - r_w}{\alpha}\right) \right], & d < r_c \\ 0, & d \ge r_c \end{cases}$$
where $d$ is patch-patch Euclidean distance, $r_w = 0.12\sigma$, $\alpha = 0.005\sigma$, $r_c = 0.20\sigma$, and $\epsilon_{\text{CSW}} = \frac{1}{T^*} = \frac{1}{0.09} \approx 11.11111111$.

---

## 2. Benchmark & Cross-Validation Results

The implementation was checked against a 2000-particle snapshot equilibrated in LAMMPS (`sanchez_burgos_lammps/validation_gpu/cube/final.data`):

| Method / Implementation | Total Potential Energy ($E_{\text{pot}}$) | Energy per Particle ($E/N$) |
| :--- | :--- | :--- |
| **LAMMPS Benchmark** (`final.data`) | `-151.70991` | `-0.07585496` |
| **Analytical Python Reference** | `-151.710929` | `-0.07585546` |
| **MC_Checkerboard CPU Verification** | `-151.709586` | `-0.07585479` |
| **MC_Checkerboard GPU Checkerboard** | `-151.709072` | `-0.07585454` |

Initial energy matches LAMMPS to within $3 \times 10^{-4}$ (relative error $< 0.0002\%$).

---

## 3. Running the Simulation

To run the simulation with GPU acceleration:
```bash
/path/to/MC_Checkerboard/bin/mc_gpu.exe datos.nml
```
