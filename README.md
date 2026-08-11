# Stab-Knock CFPS Analysis

Reproducible R code for the Stab-Knock method, empirical analyses of
China Family Panel Studies (CFPS) data, robustness and mediation analyses, and
simulation experiments.

> **Data availability:** The participant-level dataset is not included in this
> repository and is not licensed for redistribution. Only analysis code is
> released. See [`data/README.md`](data/README.md) for the expected local file.

## Repository structure

```text
.
├── R/                       # Reusable Stab-Knock and DS/MDS functions
├── analysis/                # Numbered empirical and simulation scripts
├── data/                    # Private local data; ignored by Git
├── results/                 # Generated outputs; ignored by Git
├── install_dependencies.R   # Installs required R packages
└── run_all.R                # Runs empirical analyses in dependency order
```

## Requirements

- R 4.2 or newer is recommended.
- R packages: `ggplot2`, `glmnet`, `hdi`, `knockoff`, `RColorBrewer`,
  `readxl`, and `Rfast`.

Install the dependencies from the repository root:

```bash
Rscript install_dependencies.R
```

## Private data setup

Obtain the analysis dataset through the data provider's authorized access
procedure and save the prepared workbook as:

```text
data/data.xlsx
```

The first worksheet must contain 27 columns in this order:

```text
AG, GE, HRS, HEA, EHEA, ES, JS, IRS, SWB, LS, FC, SSES, SSS, HS,
DD, FA, FMO, FHEA, FPP, MA, MMO, MHEA, MPP, FZ, NFI, REL_GDP,
REL_GDP_Per_Cap
```

The scripts assign these names programmatically, so column order is part of
the input contract. The stage-specific script additionally assumes 13,128
rows ordered by survey wave: 2,965 (2010), 3,653 (2014), 3,712 (2018), and
2,798 (2020).

## Reproducing the analyses

Run all empirical analyses from the repository root:

```bash
Rscript run_all.R
```

This executes the numbered scripts in their required order. Generated tables,
figures, and R objects are written under `results/`.

The simulation study is computationally intensive and is intentionally not
part of `run_all.R`. Run it separately with:

```bash
Rscript analysis/06_simulation.R
```

Random seeds are fixed in the scripts for reproducibility. Some numerical
results may still vary slightly across operating systems, R versions, or
package versions.

## Analysis scripts

| Script | Purpose | Dependency |
|---|---|---|
| `01_cfps_analysis.R` | CFPS dataset analysis | Private data |
| `02_stage_specific_analysis.R` | Stage-specific Stab-Knock analysis | Private data |
| `03_sensitivity_analysis.R` | IPW, placebo, and E-value sensitivity checks | Output of script 01 |
| `04_mediation_analysis.R` | Standardized-covariate mediation analysis | Private data |
| `05_additional_metrics.R` | Comparable effect and balance metrics | Output of script 01 |
| `06_simulation.R` | FDR, power, AIC, and model-size simulations | None |

## Privacy and version-control safeguards

The `.gitignore` excludes everything under `data/` except its README, along
with generated results and R session files. Before publishing, verify the Git
staging area with:

```bash
git status --short
git ls-files data
```

The second command should list only `data/README.md`.

## License

The code is released under the MIT License. The license applies to the source
code only and does not grant rights to the underlying CFPS data.

## Citation

If this repository accompanies a paper, please cite the paper and archive the
release (for example, with Zenodo) after replacing the generic repository title
and author metadata with the final bibliographic information.
