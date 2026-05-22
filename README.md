# Randkluft

**Randkluft** is an interactive Shiny application for automated, reproducible gating of marker expression in multiplexed tissue imaging data (CyCIF, CODEX, CosMX, and similar platforms).

> *Manuscript under review in **Bioinformatics**.*
> Live app: [https://irscope.shinyapps.io/Randkluft/](https://irscope.shinyapps.io/Randkluft/)

---

## Overview

In multiplexed imaging, true-positive cells (high marker intensity) are frequently obscured by a large background population of true-negative cells with low intensity. Defining a reliable, objective gating threshold is therefore non-trivial and highly operator-dependent.

Randkluft addresses this by analysing the **density profile of each marker distribution**. Using a stochastic search strategy inspired by the **Robbins–Monro algorithm**, it iteratively locates a discriminative threshold within the right-hand shoulder of the distribution — the region where rare positive cells reside — without requiring manual initialisation or distributional assumptions.

The algorithm converges when the sub-distribution below the candidate threshold reaches a skewness of zero, meaning the background population is approximately symmetric and the positive tail has been isolated.

---

## Features

| Feature | Description |
|---|---|
| Automated gating | Skewness-based threshold estimation for each marker × patient combination |
| Visual inspection | Toggle histogram overlays with estimated and manually adjusted gates |
| Manual refinement | Override any gate via numeric input; update propagates instantly |
| Bivariate gating | Scatter plot view with simultaneous two-marker gates and spatial overlay |
| Trivariate gating | Interactive 3-D scatter view for three-marker combinations |
| Phenotyping | Define Boolean marker combinations to assign cell phenotypes |
| PDF export | Download the current 4-panel plot or a full multi-page PDF for all markers |
| CSV export | Download gate estimates and phenotyped cell tables |

---

## Repository structure

```
Randkluft/
├── app.R                              # Main Shiny application
├── utils.R                            # Utility functions (sourced by app.R)
├── .renvignore                        # renv configuration
├── exemplar-001--unmicst_cell.csv     # Example cell-marker input table
├── tuulia_data_GT.csv                 # Reference gate values for benchmarking
├── phenotype_table_help.csv           # Example phenotype definition table
├── example_workflow_crevasse.csv      # Example phenotyping workflow
└── documents/                         # In-app help pages (Markdown)
    ├── home.md
    ├── help.md
    ├── faq.md
    └── contact.md
```

---

## Input format

Randkluft expects a **cell-by-marker CSV table** with the following columns:

| Column | Type | Description |
|---|---|---|
| `X_centroid` | numeric | X spatial coordinate of the cell |
| `Y_centroid` | numeric | Y spatial coordinate of the cell |
| `<MarkerName>` | numeric | One column per marker (raw or log-transformed intensity) |

- DNA / Hoechst / DAPI channels are **automatically detected and excluded** from gating.
- If values are raw (max > 20), they are **log-transformed automatically**.
- An example file (`exemplar-001--unmicst_cell.csv`) is provided and can also be downloaded from within the app.

---

## Running locally

### Prerequisites

```r
# R >= 4.3.0 recommended
install.packages("BiocManager")
BiocManager::install(c("Biobase", "BiocGenerics", "S4Vectors", "ComplexHeatmap"))

install.packages(c(
  "shiny", "shinyjs", "shinyvalidate", "shinyWidgets", "shinyalert",
  "shinydashboard", "bslib", "plotly", "ggplot2", "gridExtra", "viridis",
  "dplyr", "tidyr", "stringr", "data.table", "scales", "purrr",
  "mclust", "moments", "multimode", "umap", "tsne", "fastICA",
  "markdown", "rmarkdown", "readxl", "colourpicker", "seqinr",
  "philentropy", "PEkit", "sqldf", "tiff", "NMF"
))
```

### Launch

```r
shiny::runApp("app.R")
```

---

## Deploying to shinyapps.io

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "<your-account>",
                          token = "<your-token>",
                          secret = "<your-secret>")

rsconnect::deployApp(appDir = ".", appPrimaryDoc = "app.R",
                     appName = "Randkluft")
```

> **Note on R version:** shinyapps.io must support your local R version. Check available versions at [shinyapps.io](https://www.shinyapps.io) before deploying.

---

## Citation

If you use Randkluft in your research, please cite:

> Amiryousefi A. *et al.* **Randkluft: visually-guided detection of rare
high-intensity cells from CyCIF density skewness** *Bioinformatics* (under review).

Until the paper is formally published, please cite the application URL:
[https://irscope.shinyapps.io/Randkluft/](https://irscope.shinyapps.io/Randkluft/)

---

## Contact

For questions, bug reports, or feedback:

📧 ali_amiryousefi@hms.harvard.edu

---

## License

Randkluft is provided for academic and research use. See `LICENSE` for details.
