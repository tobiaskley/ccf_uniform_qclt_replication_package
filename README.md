# Stein's method near the unit root — Replication package

Code to reproduce the simulation figures (Figures 1–4) of the paper. There
are two independent parts: an R simulation study (Figures 1 & 2) and a
Python data example / Monte Carlo study (Figures 3 & 4).

## Repository structure

```
R_code/sim/1_Wasserstein_distance/
  sim_params.R                 simulation grid (600 configs: model x schedule x innovation x repetition)
  sim.R                        SLURM worker: simulates W1/KS distance for one config
  sim.sh                       SLURM array submission script (array 1-600), calls sim.R
  sim_merge.R                  merges the 600 per-task outputs into output/W1.Rdata
  sim_plot_endpoint_slope.R    reads output/W1.Rdata, produces Figures 1 & 2
  output/W1.Rdata              merged simulation results (included, so plotting doesn't require rerunning the simulation)

python_code/
  data/TCU_levels.csv          FRED series TCU (capacity utilization), monthly levels
  ar1_model.py                 AR(1) process classes (Gaussian / Student-t innovations)
  wasserstein_estimation.py    Wasserstein-1 distance estimator used by the notebooks
  consts.py                    moment/kurtosis constants for the Gaussian and t-distributions
  01_data_example.ipynb        Figure 3: TCU data example and AR(1) fit
  02_monte_carlo_distributions.ipynb   Figure 4: 2D grid of standardized-statistic histograms
```

## Requirements

**R** : `ggplot2`, `patchwork`.

**Python** (tested with Python 3.11+): `numpy`, `scipy`, `pandas`, `matplotlib`,
`statsmodels`, `tqdm`, `jupyter`.

## Figure 1 & 2 — MA(1) / AR(1) Wasserstein simulation (R)

Figure 1 (MA(1)) and Figure 2 (AR(1)) each combine the Gaussian- and
Student-t-innovation results as sub-panels.

The merged simulation output is already
included at `R_code/sim/1_Wasserstein_distance/output/W1.Rdata`
(600 configs x 20 sample sizes, `R = 10^7` replications, `B = 50` Monte
Carlo repetitions each) To reproduce the figures from this
file:

```bash
cd R_code/sim/1_Wasserstein_distance
Rscript sim_plot_endpoint_slope.R
```

This writes `plots/ma1_gaussian.pdf`, `plots/ma1_t.pdf` (→ Figure 1 panels),
and `plots/ar1_gaussian.pdf`, `plots/ar1_t.pdf` (→ Figure 2 panels).

**Full path (reruns the simulation):** requires a SLURM cluster.

```bash
cd R_code/sim/1_Wasserstein_distance
sbatch sim.sh              # runs sim.R as a 600-task array job, writes output/W1_<id>.Rdata
Rscript sim_merge.R        # merges into output/W1.Rdata
Rscript sim_plot_endpoint_slope.R
```

`sim.sh` is configured for our cluster (partition, module
versions, memory/walltime); adjust these `#SBATCH` settings for your own
environment before submitting.

## Figure 3 — TCU data example (Python)

```bash
cd python_code
jupyter notebook 01_data_example.ipynb
```

Fits an AR(1) model to U.S. capacity utilization (FRED `TCU`) to obtain a
realistic persistence level for the simulation study. Figure 3's two panels
are produced by the cells plotting the full series (1967–present) and the
selected estimation window (2010–2019).

## Figure 4 — Monte Carlo grid (Python)

```bash
cd python_code
jupyter notebook 02_monte_carlo_distributions.ipynb
```

Sweeps a grid of sample sizes and persistence levels, and for each cell
plots a histogram of the standardized statistic
`sqrt(n) * (gamma_hat(k) - gamma(k)) / sqrt(Sigma(k))` against the standard
normal density, alongside its estimated Wasserstein-1 distance.
