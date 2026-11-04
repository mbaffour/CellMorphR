# CellMorphR <img src="https://img.shields.io/badge/R-Shiny-blue?logo=r" alt="R Shiny"/> <img src="https://img.shields.io/badge/Python-3.8+-green?logo=python" alt="Python 3.8+"/> <img src="https://img.shields.io/badge/License-MIT-yellow" alt="MIT License"/>

## Quick start - download and run (no command line)

**[Download the latest release](https://github.com/mbaffour/CellMorphR/releases/latest)** (or clone this repo), unzip it, then double-click:

| OS | Double-click |
|----|--------------|
| **Windows** | `run_cellmorphr.bat` |
| **macOS** | `Run CellMorphR.command` |
| **Linux** | `./run.sh` |

The **first run** installs the required R packages (a few minutes); after that the app opens in your web browser instantly. Keep the console window open while you work - closing it stops the app.

### Or run it from RStudio

Open **`CellMorphR.Rproj`** (or just **`app.R`**) and click the green **Run App** button. From any R console in the folder, `shiny::runApp()` does the same.

Prefer the terminal? `Rscript launcher.R` (Windows) or `Rscript run_app.R` (macOS/Linux) works too. Every route runs the same app and installs missing packages on first use.

### Versions

| Version | App file | Status |
|---------|----------|--------|
| **3.4.0** (current) | `CellMorphR_v3.R` | Maintained. 11 analysis tabs, SuperPlot options, fraction-of-cells-beyond-threshold analysis, replicate power planner, baseline-normalised trajectories, colourblind-safe default palette with colour-vision preview, per-measurement axis units, SD/SEM/95% CI error bars, plot-data CSV export, multi-file import, click-to-edit figure text, editable PDF/SVG/PowerPoint export, dose-response trend test, anti-pseudoreplication statistics throughout. **The global models were corrected in 3.2.0 — p-values and effect sizes from 3.1.2 and earlier are not comparable. See CHANGELOG.md.** |
| 1.0.0 | `archive/CellMorphR_v1.0.R` | Archived and unmaintained - see [archive/README.md](archive/README.md). |

Every release is downloadable from the [Releases page](https://github.com/mbaffour/CellMorphR/releases).

### Try it with the example data

Two example files ship with the app - load either from the Data tab:

| File | Design |
|------|--------|
| `example_data.csv` | Time course: 4 conditions x 4 timepoints x 3 biological replicates |
| `example_data_no_timepoints.csv` | Single snapshot: 4 samples x 3 replicates, no time column |

The app also has a built-in **Guide** tab: a quick start, a "Check My Data" report that tells you in plain language what it detected and whether your design supports the statistics you want, a statistics helper, and a **Script Builder** that writes a ready-to-run R script for your analysis.

**[DATA_FORMAT.md](DATA_FORMAT.md)** explains exactly how to lay out your own file (one row per cell; `Condition` and `Replicate` required, `Time_min` optional). See [METHODS.md](METHODS.md) for a citable description of how figures are plotted and statistics computed.

**A statistical analysis suite for time-resolved single-cell morphometry data from bacteriophage infection experiments.**

CellMorphR provides publication-ready visualization and statistically rigorous analysis of cell morphological changes over time, with built-in safeguards against pseudoreplication â€” the most common statistical error in microscopy-based studies.

---

## Table of Contents

- [Overview](#overview)
- [In-app Guide and Script Builder](#scripting-your-analysis)
- [The Problem CellMorphR Solves](#the-problem-cellmorphr-solves)
- [Features](#features)
- [Installation](#installation)
  - [R Shiny App](#r-shiny-app)
  - [Python Command-Line Tool](#python-command-line-tool)
- [Data Format](#data-format)
  - [Required Columns](#required-columns)
  - [Measurement Columns](#measurement-columns)
  - [Example Data](#example-data)
  - [Handling Technical Replicates](#handling-technical-replicates)
- [R Shiny App â€” User Guide](#r-shiny-app--user-guide)
  - [Data Tab](#1-data-tab)
  - [Distributions Tab](#2-distributions-tab)
  - [Trajectories Tab](#3-trajectories-tab)
  - [Statistics Tab](#4-statistics-tab)
  - [Effect Sizes Tab](#5-effect-sizes-tab)
  - [Publication Figure Tab](#6-publication-figure-tab)
- [Python Script â€” User Guide](#python-script--user-guide)
  - [Basic Usage](#basic-usage)
  - [Command-Line Options](#command-line-options)
  - [Examples](#examples)
- [Statistical Methodology](#statistical-methodology)
  - [Why Per-Replicate Analysis?](#why-per-replicate-analysis)
  - [The Pseudoreplication Problem](#the-pseudoreplication-problem)
  - [Repeated-Measures ANOVA on Replicate Summaries](#repeated-measures-anova-on-replicate-summaries)
  - [Linear Mixed Models (Optional)](#linear-mixed-models-optional)
  - [Pairwise Comparisons](#pairwise-comparisons)
  - [Effect Size Metrics](#effect-size-metrics)
- [Visualization Philosophy](#visualization-philosophy)
  - [Violin Plots with Replicate Overlays](#violin-plots-with-replicate-overlays)
  - [Temporal Trajectory Plots](#temporal-trajectory-plots)
  - [Empirical Cumulative Distribution Functions (ECDFs)](#empirical-cumulative-distribution-functions-ecdfs)
  - [Effect Size Trajectories](#effect-size-trajectories)
  - [Composite Publication Figures](#composite-publication-figures)
- [Figure Export Specifications](#figure-export-specifications)
- [Experimental Design Considerations](#experimental-design-considerations)
- [Demo Data](#demo-data)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Overview

CellMorphR was developed to analyze single-cell morphometry data from phase-contrast microscopy experiments studying bacteriophage N4 lysis inhibition. When bacteriophage N4 infects *Escherichia coli* at high multiplicity of infection (MOI), the lysis inhibition (LIN) phenotype causes infected cells to continue growing without lysing, leading to progressive cell enlargement (filamentation) over time. Quantifying this morphological change requires tracking cell size distributions across multiple time points, conditions, and biological replicates.

The tool addresses two critical needs:

1. **Correct statistical handling of hierarchical microscopy data** â€” where hundreds of cells are measured per biological replicate, but the true sample size is the number of independent biological replicates (typically N=3).

2. **Publication-ready visualization** â€” generating figures that simultaneously show single-cell distributions (revealing biological heterogeneity) and replicate-level summaries (showing statistical power honestly).

CellMorphR is available as both an interactive R Shiny application (for
exploratory analysis and figure customization) and a Python command-line tool
(for batch processing and pipeline integration).

> **The Python command-line tool is not distributed in this repository.** This
> repository ships the R Shiny app (`CellMorphR_v3.R`) only; there is no `.py`
> file here. The Python sections below, the R/Python feature matrix and the
> cross-implementation reproducibility section describe it for reference and are
> **not** a description of anything you can run from this checkout. They also
> predate the 3.2.0 statistical corrections, so the two implementations no longer
> agree even where both are available — see
> [`AUDIT-2026-09.md`](AUDIT-2026-09.md) (CM-16).

---

## The Problem CellMorphR Solves

In a typical microscopy experiment:

```
Experiment
â”œâ”€â”€ Uninfected
â”‚   â”œâ”€â”€ 30 min
â”‚   â”‚   â”œâ”€â”€ Replicate 1 (Dish A) â†’ 3 images â†’ ~200 cells
â”‚   â”‚   â”œâ”€â”€ Replicate 2 (Dish B) â†’ 3 images â†’ ~180 cells
â”‚   â”‚   â””â”€â”€ Replicate 3 (Dish C) â†’ 3 images â†’ ~210 cells
â”‚   â”œâ”€â”€ 60 min
â”‚   â”‚   â””â”€â”€ ...
â”‚   â””â”€â”€ ...
â”œâ”€â”€ +N4 (MOI=5)
â”‚   â”œâ”€â”€ 30 min
â”‚   â”‚   â””â”€â”€ ...
â”‚   â””â”€â”€ ...
```

Each biological replicate yields ~200 measured cells. With 2 conditions Ã— 4 time points Ã— 3 replicates, you have ~4,800 total cell measurements. The critical question is: **what is your sample size?**

**The wrong answer:** N = ~4,800 cells (treating every cell as an independent observation).

**The right answer:** N = 3 biological replicates per condition Ã— time group (24 total replicate-level observations).

Cells from the same dish share the same culture conditions, the same infection event, and the same imaging session. They are **correlated**, not independent. Treating them as independent â€” a statistical error called **pseudoreplication** â€” inflates your sample size ~200-fold and produces artificially small p-values. A difference of 0.01 ÂµmÂ² can appear "significant" at p < 0.0001 simply because of the massive pseudo-N.

CellMorphR handles this correctly:

1. **Displays** all cells (violins, ECDFs) â€” to show biological heterogeneity
2. **Summarizes** at the replicate level (medians) â€” to define independent observations
3. **Tests** on replicate summaries â€” to produce honest p-values

---

## Features

| Feature | R Shiny App | Python Script |
|---------|:-----------:|:-------------:|
| Interactive data exploration | âœ“ | â€” |
| Violin plots with replicate overlays | âœ“ | âœ“ |
| Ridge plots | âœ“ | â€” |
| Box + strip plots | âœ“ | â€” |
| ECDF plots | âœ“ | âœ“ |
| Temporal trajectory plots | âœ“ | âœ“ |
| Baseline-normalised trajectories (fold change / difference) | âœ“ | â€” |
| Fraction of cells beyond a threshold | âœ“ | â€” |
| Replicate power planner | âœ“ | â€” |
| Colour-vision preview | âœ“ | â€” |
| Multi-file import | âœ“ | â€” |
| Repeated-measures ANOVA | âœ“ | âœ“ |
| Linear mixed models | âœ“ | â€” |
| Pairwise comparisons (Holm-corrected) | âœ“ | âœ“ |
| Effect size trajectories | âœ“ | âœ“ |
| Composite multi-panel figures | âœ“ | âœ“ |
| Multiple measurement variables | âœ“ (dropdown) | âœ“ (`--measure` flag) |
| Custom DPI/dimensions for export | âœ“ | âœ“ |
| PDF vector export | âœ“ | âœ“ |
| PNG raster export | âœ“ | âœ“ |
| SVG vector export | âœ“ | âœ“ |
| TIFF export (LZW compressed) | âœ“ | â€” |
| EPS export (cairo) | âœ“ | â€” |
| Batch processing | â€” | âœ“ |
| Demo data generator | âœ“ | âœ“ |

---

## Installation

### R Shiny App

**Requirements:** R â‰¥ 4.1.0, RStudio recommended

**Step 1: Install R packages**

```r
install.packages(c(
  "shiny",       # Web application framework
  "bslib",       # Bootstrap 5 theming
  "ggplot2",     # Grammar of graphics plotting
  "dplyr",       # Data manipulation
  "tidyr",       # Data tidying
  "readxl",      # Excel file import
  "DT",          # Interactive data tables
  "scales",      # Axis scaling utilities
  "lme4",        # Linear mixed models (optional, for LMM analysis)
  "lmerTest",    # Satterthwaite p-values for LMM ANOVA tables (optional, recommended with lme4)
  "emmeans",     # Estimated marginal means (optional, for LMM pairwise tests)
  "ggridges",    # Ridge plots (optional)
  "ggforce",     # Sina plots (optional, recommended)
  "patchwork",   # Multi-panel figure composition (required for Publication Figure tab)
  "writexl",     # Excel export (optional)
  "ggbeeswarm",  # Beeswarm plots (optional)
  "colourpicker",# Color selection widget (optional)
  "svglite"      # High-quality SVG export (optional, recommended)
))
```

**Step 2: Launch the app**

```r
# From RStudio: open CellMorphR_v3.R and click "Run App"
# Or from the R console:
shiny::runApp(shiny::shinyAppFile("CellMorphR_v3.R"), launch.browser = TRUE)
```

The app opens in your default web browser at `http://127.0.0.1:xxxx`.

### Python Command-Line Tool

**Requirements:** Python â‰¥ 3.8

**Step 1: Install dependencies**

```bash
pip install pandas numpy matplotlib scipy openpyxl xlrd
```

**Step 2: Verify installation**

```bash
python cellmorphr.py --demo --measure Area --all-plots
```

This generates all plots using simulated data and saves them to the current directory.

---

## Data Format

CellMorphR expects a **single-sheet** Excel file (`.xlsx`, `.xls`) or CSV/TSV file with one row per cell measurement. All cells from all conditions, time points, and replicates should be in the same file.

### Required Columns

| Column | Type | Description | Example Values |
|--------|------|-------------|----------------|
| `Condition` | String | Treatment group identifier | `Uninfected`, `+N4` |
| `Time_min` | Numeric | Time post-infection in minutes | `30`, `60`, `90`, `120` |
| `Replicate` | Integer | Biological replicate identifier | `1`, `2`, `3` |

### Measurement Columns

Any additional numeric columns are automatically detected as measurement variables. Common morphometric parameters include:

| Column | Unit | Description |
|--------|------|-------------|
| `Area` | ÂµmÂ² | Cell area from phase-contrast segmentation |
| `Perimeter` | Âµm | Cell boundary perimeter |
| `FeretMax` | Âµm | Maximum Feret diameter (longest axis) |
| `FeretMin` | Âµm | Minimum Feret diameter (shortest axis) |
| `Circularity` | Dimensionless (0â€“1) | 4Ï€ Ã— Area / PerimeterÂ²; 1.0 = perfect circle |
| `AspectRatio` | Dimensionless | FeretMax / FeretMin; higher = more elongated |

You may include as many measurement columns as needed. In the R app, a dropdown allows switching between variables. In the Python script, the `--measure` flag selects the variable.

### Example Data

```
Condition   Time_min  Replicate  Area     Perimeter  FeretMax  FeretMin  Circularity
Uninfected  30        1          2.89     6.12       2.41      1.53      0.82
Uninfected  30        1          3.21     6.55       2.58      1.62      0.79
Uninfected  30        1          2.95     6.28       2.45      1.55      0.85
...
Uninfected  30        2          3.05     6.33       2.49      1.58      0.81
...
+N4         120       3          6.44     11.23      5.82      1.41      0.64
+N4         120       3          7.12     12.01      6.31      1.38      0.58
```

### Handling Technical Replicates

A common question: *"I took 3 images per dish. How do I handle them?"*

Multiple images from the same dish are **technical replicates** â€” different fields of view of the same biological sample. They are **not** independent biological observations. All cells from all images within one biological replicate receive the **same** `Replicate` value.

```
Dish A (Replicate 1):
  Image 1 â†’ 45 cells  â”€â”
  Image 2 â†’ 62 cells  â”€â”¤â”€â”€ All labeled Replicate = 1
  Image 3 â†’ 38 cells  â”€â”˜
  Total: 145 cells, all Replicate = 1

Dish B (Replicate 2):
  Image 1 â†’ 52 cells  â”€â”
  Image 2 â†’ 71 cells  â”€â”¤â”€â”€ All labeled Replicate = 2
  Image 3 â†’ 44 cells  â”€â”˜
  Total: 167 cells, all Replicate = 2
```

The technical replicates serve their purpose by providing more cells to estimate that replicate's distribution. When CellMorphR computes the per-replicate median, all 145 cells from Dish A collapse into a single number â€” the median area for Replicate 1. That single median is your one independent observation from that dish.

**Optional column:** You may include an `Image` column to track which image each cell came from, but CellMorphR does not use it for analysis. It pools all cells within a replicate automatically.

---

## R Shiny App - User Guide

The app has **10 tabs**. You do not need all of them: load data, pick a
measurement, and the Distributions tab already gives you a publication figure.

### 1. Data Tab

Upload `.csv`, `.tsv`, `.xlsx`, or `.xls` (one row per cell - see
[DATA_FORMAT.md](DATA_FORMAT.md)), or click **Load Demo Data** to explore
immediately. Column roles are auto-detected; if detection fails you are asked to
map them, and you may choose *"(no timepoints - compare samples directly)"*.

The sidebar controls everything downstream:

| Control | What it does |
|---------|--------------|
| **Measurement Variable** | Which column all the plots and statistics use |
| **Measurement unit** | Optional, remembered per measurement (`µm²`, `µm`, `a.u.`); appended to the axis titles on the Distributions, Trajectories, Effect Sizes and Pub Figure tabs, e.g. `Area (µm²)`, or `log(Area [µm²])` under a transform |
| **Reference (control) condition** | The baseline for effect sizes and "vs control" comparisons |
| **Data transform** | None / log / sqrt, applied to every measurement |
| **Colour palette** | Okabe-Ito (colourblind-safe) by default, 8 other presets including Tol and the pre-3.3 Blue/Red, plus a per-condition colour picker |
| **Preview colours as** | Simulates deuteranopia, protanopia, tritanopia or greyscale on every figure (caption says so) so you can check a palette before export |
| **Several files at once** | Select multiple CSV/Excel files (one per image, dish or condition); they are stacked with a `Source_File` column you can map to Condition or Replicate |
| **Rename samples** | Edit condition names and timepoint labels (see Editing figures below) |
| **Split Condition into factors** | Turn `pSJExD_gp75_100nM` into a `Condition_Base` plus a numeric `Dose_nM` column, enabling factorial grouping and the trend test |
| **Save / Resume session** | Write an `.rds` you can reload later; it now carries every tab's figure text, annotations, typography and export size |
| **Save / Load figure style** | A `.json` of the figure settings alone, to reuse one look across datasets |

### 2. Distributions Tab

The main figure tab. 11 plot types: violin, strip, sina, box+strip, replicate
dot plot, ridge, ECDF, histogram, stacked quartile bar, pie, donut.

**SuperPlot options** (the recommended default presentation):

- **Distinguish replicates by shape** - each biological replicate gets its own
  point shape, so reviewers see sample-to-sample agreement directly
- **Colour individual cells by replicate** - the classic Lord et al. (2020)
  SuperPlot
- **Show n** - adds `n=cells, N=replicates` under each group
- **Choose the shapes** - pick the point shape for each biological replicate
  (filled, solid or open circles, squares, triangles, diamonds, plus, cross,
  asterisk) and a separate shape for the individual cell points. The legend
  follows your choice.

**Significance brackets** are optional and off by default. They compare
conditions within each timepoint on the replicate summaries, using the test and
correction chosen on the Statistics tab, corrected across every
condition-within-timepoint comparison as one family. The Statistics tab ships
with the matching comparison mode ("Pairwise within each time point") selected,
so at the shipped defaults the brackets and the Pairwise table are the same
family and the same numbers. Move the tab off that configuration and the two
become different families: the bracket caption then says so, and names what to
set the tab back to, instead of leaving two p-values for one contrast
unexplained. Brackets are capped per timepoint so a crowded experiment stays
readable; the caption states how many of the total comparisons are shown.

Other controls: subset which samples/timepoints appear, swap which variable sits
on the X axis, log axis, manual Y limits, label angle and wrapping.

### 3. Trajectories Tab

Replicate summaries over time, with loess or linear trend and an error band or
bars showing your choice of **SD, SEM or 95% CI** across biological replicates
(never across cells). The subtitle always names the one drawn and the number
of replicates per group. Each point is one biological replicate; dashed lines
connect replicates. **Normalise to baseline** divides (fold change) or
subtracts (difference) each replicate's own value at the earliest timepoint;
fold changes are summarised by their geometric mean with log-scale error bars,
and the subtitle names any dish dropped for lacking a baseline. **Download
Plot Data (CSV)** exports the replicate summaries together with the group
centre, SD, SEM, 95% CI and N behind the figure. With no time column this tab
explains that there is nothing to plot and points you back to Distributions.

### 3b. Thresholds Tab

The filamentation readout: for every biological replicate, the percentage of
cells whose raw measurement lies above (or below) a threshold, either a
multiple of the reference condition's pooled baseline median (default 2×) or
an absolute value. The per-dish percentages are plotted as points with the
replicate mean and SD/SEM/95% CI across dishes, optionally with significance
brackets that use the Statistics tab's test and correction on those per-dish
values. Cells are never pooled into one proportion. Table and CSV export.

### 4. Morphospace (PCA) Tab

PCA over all measurement columns. Sub-tabs: the cell-level morphospace scatter
(optionally with centroids, density contours and loading vectors), a loadings
bar chart showing which features drive PC1/PC2, centroid **trajectories** over
time, morphological **diversity** metrics (convex hull area, PC variance, mean
distance to centroid), and a scree plot. Features can be excluded from the PCA
in the sidebar.

### 5. Phenotypes Tab

Unsupervised clustering (k-means or Ward hierarchical), optionally on principal
components, with automatic k selection by silhouette. Shows the optimal-k curve,
morphospace coloured by phenotype, phenotype composition over time, and a
per-cluster characterisation table.

### 6. Statistics Tab

Build the comparison you actually want:

1. Pick the measurement and the **metadata columns that define your groups**
   (Condition alone, Condition + Time, Vector + Dose, and so on).
2. Choose the **statistical unit** - replicate medians or means.
3. Choose the **comparison mode** - each group vs a reference, all pairwise,
   pairwise within each timepoint, or pairwise within each condition.
4. Choose the **test** (Welch's t-test or Wilcoxon rank-sum) and the
   **correction** (Holm, Benjamini-Hochberg, Bonferroni, none).

**Suggest Test & Settings** inspects your design and explains what it recommends
and why, rather than silently choosing for you.

Sub-tabs: Group Preview, Comparison Preview, Test Suggestion, Descriptive
statistics, Results, **Pairwise** (with 95% CI, degrees of freedom, effect
direction and per-comparison warnings), **Trend**, EMM Plot, and About.

An optional global model is available: a **repeated-measures ANOVA** on the
replicate summaries, which puts the dish in its own error stratum so Condition
is tested against dish-to-dish variation and Time and Condition x Time against
the dish x timepoint variation, or a **linear mixed model** on cells with random
intercepts for dish *and* for dish x timepoint, with emmeans contrasts corrected
across the whole table as one family. Both name the unit each term was tested
against and report partial eta-squared or the variance components. Assumption
diagnostics (Shapiro-Wilk, Levene's test, variance ratio) accompany them and are
advisory only -- see
[Repeated-Measures ANOVA on Replicate Summaries](#repeated-measures-anova-on-replicate-summaries).

**Small-n guardrails.** With two replicates per group the app warns you that a
Wilcoxon test can never reach p < 0.05 (minimum attainable p = 0.33) and that a
t-test has about 1 degree of freedom, and points you at the mixed model or the
trend test instead.

**Trend sub-tab (dose-response).** Fits replicate summaries against a numeric
predictor (for example `Dose_nM`) within each stratum and reports slope, 95% CI,
p, R-squared and a Spearman rank check, Holm-corrected across strata. For
gradient designs this is usually both more powerful and more appropriate than
dozens of pairwise tests.

### 7. Effect Sizes Tab

Absolute difference, percent change, or **Hedges' g** against the reference
condition over time, as lines, points, or magnitude-scaled bubbles. Hedges' g
carries a 95% noncentral-t confidence interval, drawn as an error bar on every
point; with three dishes per group that interval is usually several units wide,
which is the honest answer and the reason it is shown. Timepoints where a
condition is absent, or has a single replicate, report "not estimable" with the
reason instead of stopping the tab. P-values tell you whether an effect exists;
this tab tells you how big it is and how well determined that size is.

### 8. Multi-Feature Tab

Runs the analysis across **every** measurement at once: a significance heatmap
(repeated-measures ANOVA per feature, stars or p-values in the cells), a raw-value heatmap
(optionally z-scored), a spider/radar comparison, a partial eta-squared feature
importance chart, and a divergence table giving the earliest timepoint at which
each feature separates from control.

### 9. Cell Metrics Tab

Population-level summaries (mean, median, SD, IQR, CV%, skewness, quartiles, n)
for any measurement, drawn as bar, line, violin, boxplot, ridgeline, ECDF,
scatter, hexbin, heatmap, bubble, or radar charts - with the same
replicate-level significance brackets available.

### 10. Pub Figure Tab

Assembles up to nine panels (distributions, trajectories, effect sizes,
morphospace, phenotypes, cell metrics, radar, morphological trajectory,
diversity) into one lettered composite with a single shared legend. Export it as
a figure, or as a **PowerPoint deck with one editable slide per panel**.

---

### Editing figures

**Double-click any figure** to open the editor and change:

- title, subtitle, X and Y axis labels, legend title, caption
- Y-axis range (every tab except Distributions, which has its own axis controls)
- free-text annotations, one per line as `x | y | text` (x can be a number or a
  category label such as `30 min`)
- **condition names** and **timepoint labels** (for example "0 h", "Uninduced")
- a **Reset text to automatic** button

Renaming applies everywhere at once - every figure, every statistics table,
every export. Timepoint labels are display-only: the underlying numeric time is
untouched, so ordering, trajectories and the trend test keep working. There is a
one-click "replace underscores with spaces" button for cleaning up
instrument-generated names.

Every plot tab also has sidebar controls for fonts, sizes, per-element
typography, legend position, gridlines, panel borders, point size and export
dimensions - plus **Apply to all tabs** to propagate typography and export size
everywhere.

**Journal size presets** (Nature, Science, Cell, PLOS, 16:9 slide, poster
panel) fill in width, height, 300 dpi and a journal-style base font; they are
common guidance, so check the journal's current instructions. **Preview at
export size** renders the on-screen figure at exactly the export dimensions,
through the same code path as the download, so what you see is what the PDF
will be. A **transparent background** option covers PNG, PDF, SVG and EPS.

All of these settings are stored in a saved session, and **Save Figure Style
(.json)** on the Data tab writes them alone so one look can be loaded onto
another dataset.

### Exporting

| Format | Use it for |
|--------|-----------|
| **PDF**, **SVG**, **EPS** | Vector figures for journals and for editing in Illustrator or Inkscape |
| **PPTX (editable)** | PowerPoint where text and shapes stay editable (ungroup the graphic) |
| **PNG**, **TIFF** | Raster; TIFF is LZW-compressed for journal submission |
| **CSV / TXT** | Statistics tables, cell metrics, multi-feature results |
| **RDS** | Session file to resume the analysis later |

Set width, height and DPI per tab before exporting; what you see on screen is
exactly what is written to the file.

### Scripting your analysis

The **Script** tab writes a ready-to-run R script that reproduces your current
setup outside the app - same helper functions, same statistics, same theme.
Choose your measurement, grouping, test and plot type, click **Generate script**,
then copy it or download it as a `.R` file and run it with
`Rscript my_analysis.R`. Use it to re-run an analysis when data is updated, to
batch over many measurements, or to attach a reproducible script to a paper.

`analysis_template.R` in this folder is the same starting point as a plain file.
## Python Script â€” User Guide

### Basic Usage

```bash
# Analyze a dataset with default settings (violin + trajectory plots)
python cellmorphr.py data.xlsx --measure Area

# Generate all plots and save statistics
python cellmorphr.py data.xlsx --measure Area --all-plots --save-stats

# Test with demo data
python cellmorphr.py --demo --measure Area --all-plots
```

### Command-Line Options

| Flag | Default | Description |
|------|---------|-------------|
| `file` | â€” | Input file path (Excel or CSV) |
| `--demo` | `False` | Use simulated demo data instead of a file |
| `--measure`, `-m` | `Area` | Measurement column to analyze |
| `--ref` | Auto-detect | Reference/control condition name |
| `--stat` | `median` | Summary statistic: `median` or `mean` |
| `--format`, `-f` | `pdf` | Output format: `pdf`, `png`, or `svg` |
| `--dpi` | `300` | Output resolution |
| `--font-size` | `11` | Base font size for all text |
| `--outdir`, `-o` | `.` | Output directory |
| `--all-plots` | `False` | Generate all plot types |
| `--violin` | `False` | Generate distribution plot |
| `--trajectory` | `False` | Generate temporal trajectory plot |
| `--ecdf` | `False` | Generate ECDF plot |
| `--effect` | `False` | Generate effect size plot |
| `--composite` | `False` | Generate composite publication figure |
| `--save-stats` | `False` | Save statistical output to text/CSV files |
| `--plot-style` | `violin` | Distribution plot style: `violin`, `strip`, `sina`, `box`, `dotplot` |
| `--show-sig` | `False` | Show significance brackets on distribution/composite plots |
| `--sig-style` | `stars` | Significance label: `stars` (*, **), `stars_p` (stars + p), `pval` (p only) |
| `--hide-ns` | `False` | Suppress non-significant brackets (show only \*, \*\*, \*\*\*) |
| `--seed` | `42` | Random seed for jitter reproducibility; set to 0 for non-deterministic |

### Examples

```bash
# Analyze multiple measurement variables
python cellmorphr.py data.xlsx --measure Area --all-plots --outdir results/
python cellmorphr.py data.xlsx --measure FeretMax --all-plots --outdir results/
python cellmorphr.py data.xlsx --measure Circularity --all-plots --outdir results/

# High-resolution composite figure for journal submission
python cellmorphr.py data.xlsx -m Area --composite -f pdf --dpi 600 --font-size 9

# Strip plot (jittered raw data points) â€” preferred by eLife, Nature Methods
python cellmorphr.py data.xlsx -m Area --violin --plot-style strip

# Sina plot (density-proportional jitter) â€” best of violin + strip
python cellmorphr.py data.xlsx -m Area --violin --plot-style sina

# Replicate-only dot plot â€” shows true N, follows SuperPlots framework
python cellmorphr.py data.xlsx -m Area --violin --plot-style dotplot

# Composite figure with strip plot instead of violin
python cellmorphr.py data.xlsx -m Area --composite --plot-style strip -f pdf

# Add significance brackets with stars
python cellmorphr.py data.xlsx -m Area --violin --show-sig -f pdf

# Significance brackets with p-values on composite figure
python cellmorphr.py data.xlsx -m Area --composite --show-sig --sig-style stars_p -f pdf

# Only show significant brackets (suppress ns), Nature-style
python cellmorphr.py data.xlsx -m Area --violin --show-sig --hide-ns --plot-style sina -f pdf

# PNG for presentations
python cellmorphr.py data.xlsx -m Area --composite -f png --dpi 150 --font-size 14

# Specify control condition explicitly
python cellmorphr.py data.xlsx -m Area --ref "Mock" --all-plots

# Use mean instead of median
python cellmorphr.py data.xlsx -m Area --stat mean --trajectory --save-stats
```

**Output files** are named by plot type and measurement variable:
- `distributions_Area.pdf`
- `trajectory_Area.pdf`
- `ecdf_Area.pdf`
- `effect_size_Area.pdf`
- `figure_Area.pdf` (composite)
- `statistics_Area.txt` (ANOVA and pairwise results)
- `pairwise_Area.csv` (machine-readable pairwise comparisons)

---

## Statistical Methodology

### Why Per-Replicate Analysis?

The statistical unit in a microscopy experiment is the **biological replicate** (independent culture/dish), not the individual cell. Each replicate constitutes one independent realization of the experimental conditions. The hundreds of cells measured per replicate are repeated measurements of that single biological observation â€” they improve the precision of your estimate of that replicate's cell size distribution, but they do not constitute additional independent observations.

This is analogous to measuring body weight: if you weigh yourself 10 times on the same scale, you have 1 measurement with 10 technical replicates, not 10 independent measurements.

### The Pseudoreplication Problem

Pseudoreplication occurs when non-independent observations are analyzed as though they are independent, inflating the apparent sample size and the power to detect differences. In microscopy data, this happens when individual cells from the same dish are treated as independent data points.

**Consequences of pseudoreplication:**
- **Inflated sample size**: True N=3 becomes apparent N=600+ per group.
- **Artificially small p-values**: Standard errors shrink proportionally to âˆšN; a 200-fold increase in N reduces SE by ~14-fold.
- **False positives**: Trivially small, biologically meaningless differences appear "highly significant" (p < 0.0001).
- **Misleading error bars**: SEM computed over all cells is ~14Ã— smaller than SEM computed over replicate summaries.
- **Reviewer rejection**: High-impact journals (Nature, Science, Cell) routinely flag and reject papers with pseudoreplicated microscopy analyses.

**CellMorphR's approach:**
1. Compute the summary statistic (median or mean) for each biological replicate â€” collapsing ~200 cells to 1 value.
2. Perform all statistical tests on these replicate-level summaries (N=3 per group).
3. Display both the cell-level data (for biological insight) and replicate-level summaries (for statistical honesty) in every figure.

### Repeated-Measures ANOVA on Replicate Summaries

The default global model is a repeated-measures ANOVA with factors **Condition**
and **Time** on the per-replicate medians, with the biological replicate — the
dish — as its own error stratum.

> **Changed in 3.2.0.** Up to v3.1.2 this was a plain two-way ANOVA,
> `aov(y ~ Condition * Time)`, on the same replicate summaries. In a time course
> the *same* dishes are imaged at every timepoint, so 3 dishes measured at 5
> timepoints contributed 15 rows per condition while only 3 independent units
> existed, and the Condition F-test was built on a residual divided by roughly
> five times too many degrees of freedom. On null data that produced a **41%**
> false-positive rate for Condition, where 5% is nominal. See
> [`AUDIT-2026-09.md`](AUDIT-2026-09.md) (CM-1) and [`CHANGELOG.md`](CHANGELOG.md).

**The model:**

```
median_val ~ Condition * Time + Error(Rep_unique)
```

`Rep_unique` is the dish, identified uniquely per condition (see
[unique replicate identifiers](#linear-mixed-models-optional) below). The error
stratum splits the analysis in two, because the two kinds of term have two
different error terms:

| Term | Varies | Tested against | DenDF for 2 conditions × 3 dishes × 4 times |
|---|---|---|---|
| Condition | between dishes | dish-to-dish variation | 4 |
| Time | within a dish | dish × timepoint variation | 12 |
| Condition × Time | within a dish | dish × timepoint variation | 12 |

Earlier versions of this document gave a single residual of `N - ab = 24 - 8 =
16` for that design. There are 6 dishes in it, not 24 independent units, and the
three terms do not share one error term. The Statistics tab now prints the error
stratum and the denominator df beside every term, so the unit each test actually
used is visible rather than inferred.

`Time` is deliberately not nested further (`Error(Rep_unique/Time_fac)`): there
is exactly one summary value per dish per timepoint, so that stratum would be
left with zero residual degrees of freedom and no p-values at all.

**Degenerate designs.** With a single timepoint the design collapses to a one-way
ANOVA over dishes — which already has the dish as its unit — and the output says
so. Where there are no more replicate groups than conditions there is no
dish-to-dish variation left to test Condition against, and the app reports that
instead of fitting a model.

**Effect sizes:** Partial eta-squared (η²_p) is reported for each ANOVA term,
relative to that term's own error stratum, quantifying the proportion of variance
explained by each factor relative to its own effect plus its residual variance.
Benchmarks: small = 0.01, medium = 0.06, large = 0.14 (Cohen, 1988).

**Assumption diagnostics:** CellMorphR reports three diagnostics, computed on the
same replicate-level unit the tests use, in the Statistics tab and in the
downloaded report:

1. **Normality of residuals:** Shapiro-Wilk test (Shapiro & Wilk, 1965). With
   N=3 per cell this test has low power, but ANOVA is robust to moderate
   non-normality with balanced designs (Maxwell & Delaney, 2004).
2. **Homogeneity of variance:** Levene's test with median centering (Brown &
   Forsythe, 1974) — the same test used by GraphPad Prism. Implemented in base R
   and checked against `car::leveneTest` in the test suite.
3. **Variance ratio:** the largest group SD over the smallest, which is usually
   the more readable of the two.

These are **advisory only**. CellMorphR does not switch test on the strength of
them — the test you selected is the test that runs — and the report says so
rather than implying otherwise. If an assumption looks severely violated,
consider a transform (log, for positively skewed cell size data) or the
rank-based alternative; that decision stays with you.

### Linear Mixed Models (Optional)

The R Shiny app also supports a linear mixed model (LMM) that analyzes all
cell-level data while correctly partitioning variance:

```
Measurement ~ Condition * Time + (1|Rep_unique) + (1|Rep_unique:Time_fac)
```

> **Changed in 3.2.0.** The model was previously
> `Measurement ~ Condition * Time + (1|Rep_unique)`. A random *intercept* per
> dish absorbs only what is constant across the whole time course, so the
> Condition test was correctly protected (DenDF 4 for 6 dishes) while `Time` and
> `Condition × Time` — which vary *within* a dish — fell back on the cell-level
> residual: DenDF in the thousands, i.e. the cell count. Once dishes drift over
> time, as they do in any time-lapse, that gave measured false-positive rates of
> **99%** and **100%** on null data. See
> [`AUDIT-2026-09.md`](AUDIT-2026-09.md) (CM-2).

The two random terms do different jobs:

- `(1|Rep_unique)` — cells within the same dish are correlated. This protects the
  `Condition` test.
- `(1|Rep_unique:Time_fac)` — each dish has its own offset at each timepoint
  (focus drift, medium evaporation, per-session segmentation quality). This is
  the error term `Time` and `Condition × Time` belong to, and it is what protects
  them.

The dish × timepoint term is only estimable when dishes are measured at more than
one timepoint with more than one cell per dish per timepoint. Where it is not, it
is dropped and the report says plainly what that costs, and points at the
repeated-measures ANOVA instead.

Used within its assumptions the LMM makes fuller use of the data than the ANOVA
on summaries, because it uses all cell-level information and weights replicates
with unequal cell counts. It does assume normally distributed residuals and
random effects, which may not hold for skewed cell size data — consider a log
transform.

**Singular fits.** When a variance component is estimated at exactly zero the fit
is reported, not silently replaced: the report names which term collapsed, notes
that the resulting tests are conservative rather than wrong, and points at the
repeated-measures ANOVA as the model to trust. (Up to v3.1.2 this branch returned
nothing while announcing "Falling back to Two-way ANOVA" — and no ANOVA was ever
computed, so the report promised a global model and then showed none. Audit
CM-7.)

**Important implementation detail — unique replicate identifiers:** In many
datasets, biological replicates are numbered 1, 2, 3 within each condition.
Replicate "1" in Uninfected is a *different dish* from Replicate "1" in +N4. If
replicate IDs are passed directly to the random effect term `(1|Replicate)`, the
model incorrectly treats cross-condition replicates sharing the same number as
the same random-effect group — causing model convergence failures or incorrect
variance estimates. CellMorphR automatically creates unique replicate identifiers
by combining condition and replicate number (e.g., `Uninfected.1`,
`Uninfected.2`, `+N4.1`, `+N4.2`) before fitting the model. The same identifier
is the error stratum in the repeated-measures ANOVA.

**P-values in the ANOVA table:** Base `lme4::lmer()` does not produce p-values for
fixed effects by default because the appropriate denominator degrees of freedom
for F-tests in mixed models are debated. CellMorphR uses `lmerTest::lmer()` when
available, which provides Satterthwaite-approximated degrees of freedom and
p-values — the most widely accepted approach in biological sciences. Without
`lmerTest` the table has no p-value column at all, and the report says so.
Install it alongside `lme4` for complete output:

```r
install.packages(c("lme4", "lmerTest", "emmeans"))
```

**emmeans contrasts.** Condition contrasts at each timepoint are taken from the
fitted model on Kenward-Roger degrees of freedom, pulled *unadjusted*, and then
corrected here across the whole table as one family, with the family size stated
in the output. Previously `emmeans(model, ~ Condition | Time_fac)` made the
timepoint a `by` variable, so `adjust = "holm"` corrected within each timepoint
separately — one contrast per stratum with two conditions, i.e. no correction at
all, printed under a "holm" label (audit CM-8).

### Pairwise Comparisons

When the interaction is significant, pairwise comparisons at each time point identify **when** the divergence becomes statistically significant:

- At each time point, a **Welch's two-sample t-test** (unequal variance assumed) compares the per-replicate medians between conditions (N=3 vs. N=3). Welch's test is used rather than Student's t-test because:
  - **Unequal variances are biologically expected.** Infected cells exhibit progressively higher variance as some cells filament while others do not. The equal-variance assumption of Student's t-test is violated.
  - **With N=3, equal variance cannot be verified.** Levene's test or Bartlett's test have negligible power at N=3, so assuming equal variance is unjustified.
  - **Welch's is the recommended default.** It reduces to Student's t-test when variances are equal, but does not produce inflated Type I error rates when they are not. Major journals and statistical guidelines (e.g., Ruxton 2006, *Behavioral Ecology*) recommend Welch's as the default for two-sample comparisons.
- P-values are adjusted using the **Holm-Bonferroni** method (or Benjamini-Hochberg, Bonferroni, or none, at your choice) to control the family-wise error rate across **every comparison in the table, as one family**. The family size is stated in the output. Both the Pairwise table and the emmeans contrasts define the family this way; so do the figure brackets at the shipped defaults.
- Significance is reported as: ns (p â‰¥ 0.05), * (p < 0.05), ** (p < 0.01), *** (p < 0.001).

**Expected pattern for lysis inhibition:**
- 30 min: No significant difference (infection just initiated; cells haven't yet diverged).
- 60 min: Approaching or reaching significance (infected cells beginning to enlarge).
- 90â€“120 min: Highly significant difference (progressive filamentation has created a large size gap).

### Effect Size Metrics

CellMorphR computes three effect size metrics at each time point, all on
replicate-level summaries:

1. **Absolute difference**: `mean(Treatment medians) - mean(Control medians)`. In
   the original measurement units (e.g., µm²).

2. **Percent change**: `(Treatment - Control) / Control × 100`. Scale-free and
   intuitive for reporting.

3. **Hedges' g**: `J × (Treatment - Control) / pooled_SD`, where
   `pooled_SD = √[((n₁-1)·SD_ctrl² + (n₂-1)·SD_trt²) / (n₁+n₂-2)]` and
   `J = 1 - 3/(4(n₁+n₂-2)-1)` is Hedges' small-sample correction. Reported with a
   95% **noncentral-t confidence interval**, drawn as an error bar on every point
   and carried into the table and the CSV export.

> **Changed in 3.2.0.** This was an uncorrected Cohen's d over
> `pooled_SD = √[(SD_ctrl² + SD_trt²)/2]` — the unweighted average of the
> variances, which is the pooled SD only when the two groups have equal n — with
> no interval, while the plot drew dotted "large effect" guides at ±0.8 across
> the top of it. At n = 3 vs 3, J = 0.80, so every d reported was 25% larger than
> the unbiased estimate. The ±0.8 guides remain, but only alongside the interval
> and with a caption saying they are Cohen's conventional benchmark rather than a
> threshold for your experiment. Audit CM-3.

Where a condition is absent at a timepoint, or has a single replicate, g is
reported as not estimable with the reason, and the count of such timepoints
appears in the subtitle. Previously that case reached `if (NA > 0)` and replaced
the plot, the table and the download button with an R error block (audit CM-4).

---

## Visualization Philosophy

### Distribution Plots with Replicate Overlays

CellMorphR offers 7 distribution plot styles (R Shiny) or 5 styles (Python) to suit different journals and reviewer preferences. All styles except ECDF and Dot Plot overlay per-replicate medians as individual points to convey true statistical power.

**Choosing a plot style:**

| Style | Best For | Limitations |
|-------|----------|-------------|
| **Violin** | Revealing distribution shape, multimodality, skewness | Kernel smoothing creates artificial tails; width is density estimate, not raw data |
| **Strip (Jitter)** | Maximum transparency â€” nothing hidden or smoothed | Can look noisy with >500 cells per group; hard to see density |
| **Sina** | Best of both: raw data with density information | Requires `ggforce` (R) or `scipy` (Python) |
| **Dot Plot** | Emphasizing true N; follows SuperPlots (Lord et al., 2020) | Hides within-replicate variability |
| **Box + Strip** | Familiarity with traditional reviewers | Hides distribution shape beyond quartiles |
| **Ridge** | Many groups; space-efficient | Harder to compare specific groups |
| **ECDF** | Assumption-free distribution comparison | Less intuitive for non-statistical audiences |

**Recommendation for Nature-tier journals:** Use **Sina** or **Strip** for the main figure (shows every cell without smoothing), with **Dot Plot** in supplementary to demonstrate that the true N=3 replicates are consistently separated. This follows current best practices for transparent data visualization (Weissgerber et al., 2015, *PLOS Biology*; Lord et al., 2020, *Journal of Cell Biology*).

### Temporal Trajectory Plots

**What they show:** Per-replicate summary statistics plotted over time, with condition-level trend lines and uncertainty bands.

**Why this design:**
- Each point is one biological replicate â€” the unit of statistical inference. Connecting points from the same replicate across time (dashed lines) shows individual trajectories.
- The smooth trend line (loess) captures the overall temporal pattern per condition.
- The error band or bars (SD, SEM or 95% CI, chosen in the sidebar and named in the subtitle) reflect between-replicate variability at the correct level (N = number of dishes), not the misleadingly narrow SEM that would result from computing over all cells.
- **Diverging trajectories are the visual signature of the condition Ã— time interaction** â€” the statistical result you're reporting.

### Empirical Cumulative Distribution Functions (ECDFs)

**What they show:** The cumulative proportion of cells at or below each measurement value, plotted for each condition within each time-point panel.

**Why this design:**
- ECDFs are completely non-parametric â€” no binning (unlike histograms), no bandwidth selection (unlike KDE/violins). Every cell contributes exactly one step.
- A rightward shift of the +N4 curve = larger cells. The magnitude of the horizontal gap equals the difference in that quantile of the distribution.
- The vertical gap at any x-value = the proportion of the population that differs.
- Particularly useful for detecting subpopulation effects: if only a fraction of infected cells become filamentous, the ECDF curves may overlap in the lower range but diverge in the upper tail.

### Effect Size Trajectories

**What they show:** The magnitude of the infection effect (percent change, absolute difference, or Cohen's d) plotted over time.

**Why this design:**
- Communicates the biological story directly: "Cell area increased by 30% at 60 min and by 125% at 120 min."
- P-values alone don't convey magnitude. A p = 0.01 from a 5% change and a p = 0.01 from a 125% change are very different biological results.
- The upward trend in effect size over time is the quantitative representation of progressive cell enlargement.

### Composite Publication Figures

**What they show:** Multi-panel figures (A, B, C) combining distributions, trajectories, and effect sizes.

**Why this design:**
- High-impact journals expect consolidated, multi-panel figures that tell a complete story.
- Panel A shows the raw data (distributions) â€” what the biology looks like.
- Panel B shows the statistical structure (replicate trajectories) â€” the evidence for your claim.
- Panel C shows the quantitative conclusion (effect sizes) â€” the magnitude of the effect.
- Panel labels (A, B, C) follow journal conventions for cross-referencing in the text.

---

## Figure Export Specifications

### Supported Export Formats

| Format | Type | Best For | Notes |
|--------|------|----------|-------|
| **PDF** | Vector | Journal submission (default) | Text remains searchable and editable. Scales to any size without quality loss. Required by Nature, Science, Cell. |
| **SVG** | Vector | Post-editing in Illustrator/Inkscape | Ideal for adding annotations, adjusting colors, or combining with schematics. Uses `svglite` if installed for cleaner output. |
| **EPS** | Vector | Legacy journal systems | Uses `cairo_ps` device for reliable font embedding. Some older submission systems require EPS. |
| **PNG** | Raster | Presentations, web, quick sharing | White background. Set DPI â‰¥ 300 for print quality. |
| **TIFF** | Raster | Journals requiring TIFF specifically | LZW lossless compression to keep file sizes reasonable. Some journals (e.g., PLOS, Wiley) specifically request TIFF. |

### Journal Requirements

For submission to Nature, Science, Cell, and related journals:

| Parameter | Single Column | Double Column | Full Page |
|-----------|:------------:|:------------:|:---------:|
| Width | 3.5 in (89 mm) | 7.0 in (178 mm) | 10.0 in (254 mm) |
| Min DPI (photos) | 300 | 300 | 300 |
| Min DPI (line art) | 600 | 600 | 600 |
| Preferred format | PDF or EPS | PDF or EPS | PDF or EPS |
| Font | Sans-serif (Arial) | Sans-serif (Arial) | Sans-serif (Arial) |
| Min font size | 6 pt | 6 pt | 6 pt |
| Max font size | 8 pt | 8 pt | 12 pt |

CellMorphR defaults (300 DPI, PDF output, sans-serif fonts, 11 pt base) meet these requirements. Adjust font size downward (8â€“9 pt) for single-column figures.

---

## Experimental Design Considerations

CellMorphR is designed for the following experimental structure:

| Component | Expected Values |
|-----------|----------------|
| Conditions | 2+ (e.g., Uninfected, +N4, +N4-mutant) |
| Time points | 2+ numeric values (e.g., 30, 60, 90, 120 min) |
| Biological replicates | â‰¥3 per condition Ã— time combination |
| Cells per replicate | ~50â€“500 (no minimum; unequal N is handled correctly) |
| Technical replicates (images) | Any number; pooled within biological replicates |

**Minimum requirements:**
- At least 2 conditions and 2 time points (otherwise, the interaction test is not applicable).
- At least 3 biological replicates per group (for meaningful variance estimation).
- At least ~30 cells per replicate (for stable median/mean estimation).

**Limitations:**
- The repeated-measures ANOVA assumes equal variance across groups (homoscedasticity), and sphericity across timepoints, for which no correction is applied. Cell size distributions in phage infection experiments often become more variable over time in infected conditions. Log-transformation may help if this is severe.
- With N=3 replicates per group the global model has 4 denominator degrees of freedom for Condition, not 16, and power is correspondingly limited. This is the honest number. Prioritize effect sizes with their intervals, and biological consistency, over p-values.
- The app does not currently support more complex designs (e.g., nested random effects for imaging plates, split-plot designs, or longitudinal within-replicate tracking).

---

## Cross-Implementation Reproducibility

> **Out of date as of 3.2.0, and unverifiable from this repository.** The
> Python script is not present in this checkout (see the note in
> [Overview](#overview)), and the R implementation's global models, effect sizes
> and multiplicity families all changed in 3.2.0. Any Python implementation
> written against the description below will **not** reproduce current R output.
> Treat this section as a record of the pre-3.2.0 intent until the Python tool is
> shipped and revalidated against the corrected R models.

The intent was that both the R Shiny app and the Python command-line tool
implement identical statistical methods and produce identical numerical results
when given the same input data, so that the choice of tool does not affect
reported findings.

### Verified Equivalences

| Component | R Implementation | Python Implementation |
|-----------|-----------------|----------------------|
| ANOVA | `aov(Y ~ Condition * Time)` (Type I SS) | Manual Type I SS computation |
| Pairwise t-tests | `t.test(trt, ctrl)` (Welch's, `var.equal=FALSE`) | `scipy.stats.ttest_ind(trt, ctrl, equal_var=False)` |
| Multiple comparison correction | `p.adjust(method="holm")` | Manual Holm-Bonferroni with monotonicity enforcement |
| Summary statistic | `median()` per replicate | `numpy.median()` per replicate |
| P-value computation | R's built-in F-distribution | `scipy.stats.f.sf()` (survival function) |

### Why Both Produce Identical Results

- **ANOVA:** With balanced designs (equal N per cell), Type I, Type II, and Type III sums of squares are all mathematically identical. The Python script automatically verifies design balance and prints a confirmation or warning. If your design is unbalanced (e.g., a missing replicate), the Python output will flag this and recommend switching to the R Shiny app's LMM option or using `car::Anova(model, type='III')` for Type III SS.
- **Welch's t-test:** Both implementations use Welch's t-test with Satterthwaite-approximated degrees of freedom. The argument order (treatment, control) is matched to produce consistent t-statistic signs.
- **Holm correction:** Both apply the same step-up procedure with monotonicity enforcement, producing identical adjusted p-values.

To verify reproducibility on your own data, run the analysis in both tools and compare the ANOVA F-statistics, pairwise t-statistics, and Holm-adjusted p-values. They should match to at least 3 decimal places.

---

## Demo Data

Both the R and Python implementations include a built-in demo data generator that produces simulated bacteriophage infection morphometry data. The simulation mimics the expected biology:

- **Uninfected cells**: Log-normally distributed cell area centered around ~3 ÂµmÂ² with slow, minimal growth over time. Between-replicate variability is ~3%.
- **+N4 infected cells**: Cell area increases progressively from ~3 ÂµmÂ² at 30 min to ~7 ÂµmÂ² median at 120 min, with increasing variance over time (some cells filament more than others). This models the lysis inhibition phenotype where cells continue growing without lysing.
- **Sample sizes**: ~200 cells per replicate, 3 replicates per condition Ã— time point, for a total of ~4,800 cells.

The demo data produces clear, significant results (Condition Ã— Time interaction p < 0.0001, pairwise significance emerging at 60 min) and serves as a reference for validating your analysis pipeline.

**R:**
```r
# Click "Load Demo Data" button in the Data tab sidebar
```

**Python:**
```bash
python cellmorphr.py --demo --measure Area --all-plots --save-stats
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `Error: could not find function "%\|\|%"` | Fixed in current version. If using an older version, add `library(rlang)` to the top of `CellMorphR_v3.R`, or replace the `%\|\|%` operator with base R (see commit history). |
| LMM crashes / app disconnects when running Linear Mixed Model | Fixed in current version. Two issues were resolved: (1) Replicate IDs like 1, 2, 3 repeated across conditions caused `lmer` to treat cross-condition replicates as the same random-effect group â€” the app now creates unique IDs automatically. (2) Missing `lmerTest` package meant ANOVA tables lacked p-values, causing downstream errors. Install with `install.packages("lmerTest")`. |
| LMM ANOVA table shows no p-values | Install `lmerTest`: `install.packages("lmerTest")`. Without it, `lme4::lmer` returns F-statistics but no p-values by design. `lmerTest` adds Satterthwaite-approximated p-values. |
| `Error in read_excel: path does not exist` | Ensure the file extension is `.xlsx` or `.xls`. If your file is `.csv`, the app detects this automatically. |
| Measurement dropdown is empty | Your file may lack numeric columns beyond the metadata. Check that `Area`, `Perimeter`, etc. are numeric (not stored as text). |
| Publication Figure tab is blank | Install the `patchwork` package: `install.packages("patchwork")`. This package is required for multi-panel figure composition. |
| Linear Mixed Model option grayed out | Install `lme4`, `lmerTest`, and `emmeans`: `install.packages(c("lme4", "lmerTest", "emmeans"))`. |
| Ridge plot not available | Install `ggridges`: `install.packages("ggridges")`. |
| Sina plot falls back to regular jitter | Install `ggforce`: `install.packages("ggforce")`. Without it, the sina plot uses regular jitter as a fallback. |
| Figures look different across platforms | PDF output is vector-based and renders identically everywhere. PNG rendering may vary slightly by system. |
| Python: `ModuleNotFoundError: No module named 'openpyxl'` | Install with `pip install openpyxl`. Required for reading `.xlsx` files. |
| Very slow with large datasets (>50,000 cells) | The violin plots are the bottleneck. Consider subsampling for exploration and using the full dataset for final figures. |
| Unequal cell counts across replicates | This is handled correctly. The per-replicate median/mean is computed regardless of cell count. Replicates with 100 cells and replicates with 300 cells contribute equally to the ANOVA (one data point each). |
| TIFF files are very large | TIFF uses LZW lossless compression by default. For smaller files, consider PNG (lossy) or PDF (vector). Some journals require TIFF specifically; in that case, reduce DPI to 300 for photographic content. |
| SVG text renders incorrectly | Install `svglite` for better SVG output: `install.packages("svglite")`. Without it, the base R SVG device is used, which may have font issues on some systems. |
| EPS fonts look wrong | CellMorphR uses `cairo_ps` for EPS export, which embeds fonts properly. If fonts still look wrong, ensure Cairo is installed on your system (`capabilities("cairo")` in R should return `TRUE`). |
| Loess warning with few time points | With fewer than 4 time points, loess smoothing is automatically replaced with a linear trend line. This is expected behavior, not an error. |

---

## Citation

If you use CellMorphR in your research, please cite:

```bibtex
@software{cellmorphr2026,
  author  = {Awuah, Michael},
  title   = {CellMorphR: Single-Cell Morphometry Analysis Suite for Bacteriophage Infection Experiments},
  year    = {2026},
  url     = {https://github.com/mbaffour/CellMorphR},
  version = {3.4.0}
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contributing

Contributions are welcome. Please open an issue or submit a pull request on [GitHub](https://github.com/mbaffour/CellMorphR).

**Areas for future development:**
- Support for longitudinal cell tracking (within-replicate time series)
- Additional measurement variables from fluorescence microscopy
- Automated outlier detection and quality control
- Integration with CellProfiler and ImageJ/FIJI output formats
- Bayesian hierarchical models as an alternative to frequentist ANOVA
