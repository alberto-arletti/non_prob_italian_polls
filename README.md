# Adjusting Non-Probability Samples Under Weak Auxiliary Information

[![DOI](https://zenodo.org/badge/1052735401.svg)](https://doi.org/10.5281/zenodo.21200506)

Replication code and shareable data for the paper
**"Adjusting Non-Probability Samples Under Weak Auxiliary Information: An Empirical Benchmark from the 2022 Italian Elections"**
by Alberto Arletti, Omar Paccagnella, and Maria Letizia Tanturri.

This repository benchmarks seven adjustment methods for non-probability survey
samples against a known population target — the official results of the 2022
Italian national elections. Each method is applied to five polls across four
census-based post-stratification tables, and evaluated on bias, MSE, and
interval coverage.

> **Note on data.** The individual-level survey data are proprietary. Only
> **dataset 5** is redistributed here, with the permission of Demetra Opinioni
> s.r.l. See [Data availability](#data-availability) below.

## Repository structure

```
non_prob_italian_polls/
├── data_rev/
│   └── sample.Rdata     # study variables used for adjustment (dataset 5 only, shared with permission)
├── census/
│   ├── edu.csv          # region × gender × age × education
│   ├── emp.csv          # region × gender × age × employment
│   ├── pastvote.csv     # region × 2018 coalition vote
│   ├── approx.csv       # region × gender × age × education × employment (Multiscopo approximation)
│   └── GT.csv           # official 2022 national election results (country-wide)
├── helper.R             # core functions used for adjustment
├── main.Rmd             # end-to-end workflow that runs the adjustment and produces results
└── README.md
```

The four post-stratification tables correspond to the `edu`, `emp`, `pastvote`,
and `approx` tables described in the paper. Please confirm the exact file names
in `census/` match your local copy.

## Data availability

- **Dataset 5** (`data_rev/sample.Rdata`) is included for replication, shared
  with the permission of Demetra Opinioni s.r.l.
- **Datasets 1–4** are proprietary and are **not** distributed here. 
- **Census tables** in `census/` are derived from publicly available Istat
  dissemination systems (I.Stat, IstatData, Permanent Population and Housing
  Census) and from the Italian Ministry of the Interior electoral archive
  (Eligendo). `GT.csv` reports the official 2022 election result.

## Requirements

The analysis runs in R (developed under R ≥ 4.2 — please pin your version).
Key packages by role:

- **Adjustment / modelling:** `survey`, `nnet`, `xgboost`, `nonprobsvy`,
  `rstanarm`, `torch`
- **Data handling and plotting:** `tidyverse` (incl. `dplyr`, `ggplot2`)

Install the core packages with:

```r
# packages needed to run the adjustment workflow
install.packages(c(
  "survey", "nnet", "xgboost", "nonprobsvy",
  "rstanarm", "torch", "tidyverse"
))
# torch also needs its backend downloaded on first use:
torch::install_torch()
```

For a fully reproducible environment, consider committing an `renv.lock` file
(via [`renv`](https://rstudio.github.io/renv/)) or a `sessionInfo()` dump.

## Reproducing the analysis

1. Clone the repository and open it as an R project.
2. Ensure the packages above are installed.
3. Open `main.Rmd` and knit it (or run the chunks in order). It sources
   `helper.R`, loads `data_rev/sample.Rdata` and the tables in `census/`, runs
   the seven adjustment methods, and reproduces the metrics and figures.

Because only dataset 5 is public, the shipped workflow reproduces the results
for that dataset; the full multi-dataset results in the paper require the
restricted data.

## Citation

If you use this code or data, please cite the paper:

```bibtex
@software{arletti2025code,
  author  = {Arletti, Alberto and Paccagnella, Omar and Tanturri, Maria Letizia},
  title   = {non\_prob\_italian\_polls: Replication code and data},
  year    = {2026},
  version = {1.0.0},
  doi     = {10.5281/zenodo.21200507},
  url     = {https://doi.org/10.5281/zenodo.21200507}
}
```

Update the entry with the journal, volume, and DOI once published.

## License

No license is currently specified, which by default means the code is **not**
openly licensed for reuse. To allow reuse, add a license file — for example
[MIT](https://choosealicense.com/licenses/mit/) for the code. Note that the
data are shared under the separate permission described above and are not
covered by a code license.

## Contact

Alberto Arletti — Department of Statistical Sciences, University of Padova.
For questions about the code, please open an issue on this repository.
