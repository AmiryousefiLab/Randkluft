# Randkluft — Vignette

**Automated unitary gating of CyCIF markers**

---

## Introduction

This vignette walks through a complete Randkluft analysis using the bundled example dataset (`exemplar-001--unmicst_cell.csv`), from file upload through gate estimation, visual inspection, manual refinement, phenotyping, and export.

Randkluft is built around a single core idea: in marker distributions from multiplexed imaging, true-positive cells form a small right-skewed tail attached to a large symmetric background population. The algorithm isolates this tail by iteratively trimming the upper bound of the distribution until the remaining sub-distribution reaches zero skewness — the point at which only the background remains. The trim boundary at that moment is the gate.

---

## 1. Launching the app

### Locally

```r
shiny::runApp("app.R")
```

### Online

[https://irscope.shinyapps.io/Randkluft/](https://irscope.shinyapps.io/Randkluft/)

Navigate to the **Randkluft** tab to begin.

---

## 2. Uploading data

### 2.1 File format

Randkluft expects a **cell-by-marker CSV file**. Each row represents one cell, and columns contain:

- `imageid` — sample or patient identifier
- `X_centroid`, `Y_centroid` — spatial coordinates
- One column per marker (e.g. `CD3`, `CD20`, `SOX10`, `ELANE`)

```
imageid,X_centroid,Y_centroid,CD3,CD20,SOX10,ELANE,...
LSP11350,312.4,185.2,4.1,3.8,12.7,5.1,...
LSP11350,314.1,188.0,3.9,4.2,8.4,4.8,...
...
```

Marker values may be raw fluorescence intensities or pre-log-transformed values. Randkluft automatically log-transforms raw data (values > 20) on upload.

DNA / Hoechst / DAPI channels are detected by column name pattern and masked automatically.

### 2.2 Uploading the example file

1. Open the **Upload file** tab in the sidebar.
2. Click **Browse** and select `exemplar-001--unmicst_cell.csv`.
3. Alternatively, click the **here** link to download this file from within the app.
4. Wait for the notification **"File Ready for Use"** to appear in the bottom-right corner.

The app will:
- Detect the sample identifier column (`imageid`)
- Log-transform intensities if needed
- Auto-populate the marker checklist and patient selector

---

## 3. Running automated gating

1. Navigate to the **Randkluft** sub-tab inside the sidebar (**Essential** section).
2. All detected markers are pre-selected. Deselect any markers you do not wish to gate.
3. All patients are pre-selected. Deselect specific patients if needed.
4. Click the **Randkluft** button.

A progress bar labelled **"Randkluft in Action…"** tracks completion. When finished, a confirmation dialog reads **"Randkluft Found"**.

### What the algorithm does

For each marker × patient pair:

1. Remove top and bottom 1% outliers.
2. Compute the skewness of the full distribution.
3. If skewness < 0 (negatively skewed), fall back to a two-component Gaussian Mixture Model gate.
4. Otherwise, perform binary search on the upper bound `b` of the distribution, shrinking it until `skewness(data[data < b]) ≈ 0` (tolerance α = 0.01, maximum 100 iterations).
5. The converged `b` is returned as the gate.

---

## 4. Visual inspection

Enable the **Show gates** toggle to activate histogram rendering.

Each view shows a 4-panel layout:

| Panel | Content |
|---|---|
| Top-left | Marker histogram with the estimated gate (red line) |
| Top-right | Digital representation — spatial map coloured by marker intensity |
| Bottom-left | Positive cells — red/grey spatial overlay based on gate |
| Bottom-right | Positive density — kernel density contour of gate-positive cells |

Use the **Next/Previous Patient** and **Next/Previous Marker** buttons to page through all combinations.

---

## 5. Manual gate refinement

If the automated gate is unsatisfactory for a particular marker:

1. With the histogram visible, type a new threshold into the **Type Gate Value** field.
2. Click **Update Gate**.
3. The histogram and all four panels update immediately with the new gate.
4. The revised gate is saved into the gate table, and is reflected in all downloads.

---

## 6. Bivariate and trivariate views

### Bivariate

Navigate to the **Bivariate** sub-tab. Select X and Y markers and a patient. A scatter plot appears showing:

- Per-cell density (colour-coded)
- Vertical and horizontal gate lines
- Quadrant proportions annotated in red (+/+), green (+/−), blue (−/+), and black (−/−)

Below the scatter plot, per-marker spatial overlays show the gate-positive cells for each axis separately.

Update either gate using the **Enter Gate** fields and **Update Gates** button.

### Trivariate

The **Trivariate** sub-tab provides an interactive 3-D scatter plot for three simultaneous markers.

---

## 7. Phenotyping

Navigate to **Extra → Phenotyping**.

### 7.1 Upload a phenotype workflow

Upload a two-column CSV defining cell types and their marker criteria:

```
phenotype,markers
"CD4+ T cells","CD3+, CD4+, CD45+"
"Cytotoxic T cells","CD3+, CD8+, CD45+"
"B cells","CD20+, CD45+"
"Macrophages","CD68+, CD45+"
```

An example file (`phenotype_table_help.csv`) is provided. Click **here** inside the app to download it.

### 7.2 Define phenotypes manually

Alternatively, use the checkbox interface:

1. Check the markers that define a phenotype.
2. Toggle the switch beside each marker to `+` (positive) or `−` (negative).
3. Enable **Any positive** if at least one positive marker is sufficient (OR logic instead of AND).
4. Type a phenotype name.
5. Click **Add phenotype definition**.

Repeat for each cell type.

### 7.3 Apply phenotypes

Click **Phenotype my data**. The app assigns phenotypes to every cell, annotates the original data table with a `phenotype` column, and displays:

- A summary table of phenotype definitions
- A stacked bar chart of phenotype proportions
- A **Partition Diversity Estimate** (MLE) summarising compositional diversity

### 7.4 Export

- **Download Phenotyped Data** — the original CSV with a `phenotype` column appended.
- **Download Workflow** — the phenotype definition table for reproducibility.

---

## 8. Exporting results

### Gate estimates

Click **Download Gate Estimates** to download a CSV table:

```
Patient,Marker,Gate
LSP11350,CD3,8.12
LSP11350,CD20,7.84
...
```

### Histogram plots

- **Download Current Plot** — downloads a PDF of the four-panel view currently displayed (histogram + spatial panels).
- **Download All Plots** — downloads a multi-page PDF with one four-panel page per marker × patient combination.

> **Note:** "Show gates" must be toggled ON and plots must be visible before "Download Current Plot" is available.

---

## 9. Worked example

The bundled file `exemplar-001--unmicst_cell.csv` contains a subset of the MCMICRO exemplar dataset (CyCIF, Harvard Medical School). It has a single patient (`Patient1`) and 27 marker channels.

**Suggested workflow:**

```
1. Upload exemplar-001--unmicst_cell.csv
2. Click Randkluft with all markers selected
3. Toggle Show gates ON
4. Page through markers using Next Marker
5. Identify any gate that looks misplaced (e.g. DNA channels if not masked)
6. Adjust manually using Type Gate Value
7. Download Gate Estimates CSV
8. Upload phenotype_table_help.csv in the Phenotyping tab
9. Click Phenotype my data
10. Download Phenotyped Data
```

---

## 10. Tips and troubleshooting

| Symptom | Likely cause | Solution |
|---|---|---|
| Proportions or cell counts are zero | Marker is negatively skewed — no positive tail | Inspect manually and set gate by hand |
| "Disconnected from server" on upload | File does not match expected format | Check that `imageid`, `X_centroid`, `Y_centroid` columns exist |
| Slow upload or gating | File has many unused columns | Remove spatial or DNA columns before upload |
| Download Current Plot produces no file | Show gates not yet activated | Toggle Show gates ON and wait for plots to render |
| Gates look identical across markers | All markers already log-transformed and near zero | Confirm data is correctly formatted |

---

## Session info

```r
sessionInfo()
# R version 4.4.x or later recommended
# Key packages: shiny, ggplot2, gridExtra, dplyr, moments, multimode, mclust
```

---

## Citation

> Amiryousefi A. *et al.* **Randkluft: automated unitary gating of CyCIF markers.** *Bioinformatics* (under review).

Contact: ali_amiryousefi@hms.harvard.edu
