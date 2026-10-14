# MCCB-gpu Manuscript: Figure Reproduction and Publication Package

This archive contains all the source files, plotting scripts, raw/processed datasets, and documentation required to **independently reproduce all figures and compile the manuscript** for:

> **"Massively Parallel Checkerboard Monte Carlo on CUDA Architectures: Accelerating Dense Fluids, Hybrid MC/MD, and Aggregation Kinetics"**  
> *Enrique Lomba, Antonio Díaz-Pozuelo, and Eva G. Noya*  
> *Computer Physics Communications* (CPC), 2026.

---

## 1. Quick Start: Reproducing All Figures

### Prerequisites
The plotting suite requires **Python 3.8+** with `numpy` and `matplotlib`. No proprietary or GUI dependencies are required (all scripts use the headless `Agg` backend).

Install dependencies:
```bash
pip install -r requirements.txt
```

### One-Click Reproduction
To re-generate all publication figures in high-resolution vector PDF and 300 DPI PNG formats, simply run:
```bash
bash reproduce_all_figures.sh
```
or directly via Python:
```bash
python3 reproduce_all_figures.py
```

The script will automatically parse all datasets in `data/`, generate the figures, verify output integrity, and store the resulting files directly in `figures/`.

---

## 2. Package Structure

```
Manuscript/
├── README.md                      # This comprehensive documentation guide
├── requirements.txt               # Minimal Python dependencies (numpy, matplotlib)
├── reproduce_all_figures.sh       # One-click shell runner to generate all figures
├── reproduce_all_figures.py       # Python master reproduction runner
├── compile_manuscript.sh          # Helper script to compile the LaTeX manuscript
├── mccb_gpu_cpc.tex               # Complete LaTeX source of the manuscript (elsarticle)
├── references.bib                 # BibTeX database with complete references
├── mccb_gpu_cpc.bbl               # Pre-compiled bibliography
├── mccb_gpu_cpc.pdf               # Pre-compiled 12-page camera-ready manuscript
│
├── figures/                       # Destination folder for all publication figures
│   ├── consistency_all.pdf        # Fig 2: 6-panel thermodynamic consistency
│   ├── consistency_hs.pdf         # Fig 2a: Hard spheres individual panel
│   ├── consistency_lj.pdf         # Fig 2b: Lennard-Jones individual panel
│   ├── consistency_ssp.pdf        # Fig 2c: Site-site patchy individual panel
│   ├── scaling_benchmark.pdf      # Fig 3: CPU vs GPU scaling and throughput speedup
│   ├── hmc_overhead_benchmark_comparison.pdf # Fig 4: In-memory HMC profiling
│   ├── avbmc_benchmark_comparison.pdf        # Fig 5: AVBMC kinetics and acceptance
│   ├── avbmc_csd_distribution.pdf # Fig 6: Final cluster size distribution (CSD)
│   └── avbmc_benchmark_all.pdf    # Unified 3-panel AVBMC composite figure
│
├── scripts/                       # Modular, fully independent plotting scripts
│   ├── plot_fig2_consistency.py   # Script for Figure 2 (HS, LJ, SSP consistency)
│   ├── plot_fig3_scaling.py       # Script for Figure 3 (gpMC CPU vs MCCB GPU)
│   ├── plot_fig4_hmc_overhead.py  # Script for Figure 4 (LAMMPS HMC overhead)
│   └── plot_fig5_fig6_avbmc.py    # Script for Figures 5 & 6 (AVBMC kinetics & CSD)
│
└── data/                          # Standalone raw and processed benchmark datasets
    ├── consistency/               # Thermodynamic consistency trajectory and run logs
    │   ├── hs/ (npt, nvt)         # Hard spheres run-data.dat and output.out
    │   ├── lj/ (npt, nvt)         # Lennard-Jones run-data.dat
    │   └── ssp/ (npt, nvt)        # Site-site patchy run-data.dat
    ├── scaling/                   # gpMC vs MCCB-gpu runtime benchmarks across ensembles
    │   └── scaling_benchmark_data.csv
    ├── hmc/                       # Direct in-memory C-API HMC profiling CSV tables
    │   ├── data_transfer_breakdown.csv
    │   ├── trajectory_scaling_analytic.csv
    │   └── trajectory_scaling_table.csv
    └── avbmc/                     # Patchy colloidal cluster aggregation trajectories
        ├── run_with_avbmc/        # clusevol_mc.dat, output.out, mclast_clconf.lammpstrj
        └── run_without_avbmc/     # clusevol_mc.dat, output.out, mclast_clconf.lammpstrj
```

---

## 3. Mapping: Manuscript Figures vs Scripts and Datasets

| Manuscript Reference | Description | Plotting Script | Input Data Sources | Output Figures (`figures/`) |
| :--- | :--- | :--- | :--- | :--- |
| **Figure 2** | Thermodynamic consistency across $NpT$ and $NVT$ ensembles for Hard Spheres, Lennard-Jones, and Site-Site Patchy systems ($N = 7{,}695$). | `scripts/plot_fig2_consistency.py` | `data/consistency/hs/`, `data/consistency/lj/`, `data/consistency/ssp/` | `consistency_all.pdf`<br>`consistency_all.png` |
| **Figure 3** | Wallclock execution time and throughput speedup factor comparing serial CPU (`gpMC`) against GPU (`MCCB-gpu`) from $N = 8\text{k}$ to $125\text{k}$. | `scripts/plot_fig3_scaling.py` | `data/scaling/scaling_benchmark_data.csv` | `scaling_benchmark.pdf`<br>`scaling_benchmark.png` |
| **Figure 4** | In-memory HMC coupling profiling on NVIDIA RTX PRO 4500 (287× RAM vs disk transfer, trajectory scaling, task overhead). | `scripts/plot_fig4_hmc_overhead.py` | `data/hmc/*.csv` | `hmc_overhead_benchmark_comparison.pdf`<br>`hmc_overhead_benchmark_comparison.png` |
| **Figure 5** | Patchy colloidal cluster sampling kinetics: maximum cluster size evolution $S_{\max}(t)$ and move acceptance rates with vs without AVBMC. | `scripts/plot_fig5_fig6_avbmc.py` | `data/avbmc/run_with_avbmc/`<br>`data/avbmc/run_without_avbmc/` | `avbmc_benchmark_comparison.pdf`<br>`avbmc_benchmark_comparison.png` |
| **Figure 6** | Final cluster size distribution (CSD) $s \cdot n(s)/N$ comparing canonical MC translations against accelerated AVBMC cluster exchanges. | `scripts/plot_fig5_fig6_avbmc.py` | `data/avbmc/*/mclast_clconf.lammpstrj` | `avbmc_csd_distribution.pdf`<br>`avbmc_csd_distribution.png` |

*(Note: Figure 1 is rendered natively within the LaTeX source via TikZ and requires no external graphic file).*

---

## 4. Running Individual Scripts Independently

Each script in `scripts/` can be executed completely independently from any terminal directory. Default paths are automatically resolved relative to the script location, or can be customized via command-line arguments:

```bash
# Example: Generate Figure 2 (Thermodynamic Consistency)
python3 scripts/plot_fig2_consistency.py --data-dir data/consistency --out-dir figures

# Example: Generate Figure 3 (Scaling Benchmarks)
python3 scripts/plot_fig3_scaling.py --data-dir data/scaling --out-dir figures

# Example: Generate Figure 4 (HMC Overhead Breakdown)
python3 scripts/plot_fig4_hmc_overhead.py --data-dir data/hmc --out-dir figures

# Example: Generate Figures 5 & 6 (AVBMC Kinetics & CSD)
python3 scripts/plot_fig5_fig6_avbmc.py --data-dir data/avbmc --out-dir figures
```

---

## 5. Compiling the LaTeX Manuscript

To compile `mccb_gpu_cpc.tex` into the camera-ready 12-page PDF:
```bash
bash compile_manuscript.sh
```
or manually with `pdflatex`:
```bash
pdflatex mccb_gpu_cpc.tex
bibtex mccb_gpu_cpc
pdflatex mccb_gpu_cpc.tex
pdflatex mccb_gpu_cpc.tex
```

---

## 6. Portability and Independence

- **No Absolute Paths**: All scripts resolve data and output paths dynamically using relative offsets from the script file or command line arguments.
- **Headless Compatibility**: All scripts explicitly enforce `matplotlib.use('Agg')` before graphical imports, guaranteeing successful execution on remote clusters, compute nodes without X11 displays, and Docker containers.
- **Publication Styling**: Figures use serif fonts, high-contrast accessible color palettes, and dual format export (vector PDF for LaTeX embedding and 300 DPI PNG for web/presentations).
