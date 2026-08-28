# How to format your data for CellMorphR

CellMorphR reads **CSV, TSV, or Excel (.xlsx/.xls)** files in *long format*:
**one row per cell**.

Two ready-to-use examples ship with the app:

| File | Design |
|------|--------|
| `example_data.csv` | Time course: 4 conditions x 4 timepoints x 3 replicates |
| `example_data_no_timepoints.csv` | Single snapshot: 4 samples x 3 replicates, no time column |

---

## Required columns

| Column | What it is | Example values |
|--------|------------|----------------|
| **Condition** (or `Sample`, `Group`, `Treatment`) | The experimental group the cell belongs to | `WT`, `Phage_MOI_1`, `pSJExD_gp75_100nM` |
| **Replicate** (or `Rep`, `Experiment`, `Bio_Rep`) | The **biological** replicate. This is your true *n* - every p-value in the app is computed across these, never across cells | `1`, `2`, `3` |

Plus **at least one numeric measurement column** (see below).

## Optional columns

| Column | What it is | Notes |
|--------|------------|-------|
| **Time_min** (or `Time`, `Timepoint`, `Hour`, `Frame`) | Timepoint of the image | **Optional.** Leave it out entirely for single-snapshot experiments - the app then simply compares your samples side by side. |
| `Cell_ID`, `ID`, `Image`, `Label`, `Index` | Identifiers | Ignored automatically, never treated as measurements |

Column names are matched case-insensitively, and near-misses are detected
(`sample`, `treatment`, `bio_rep`, ...). If a role cannot be detected, the app
asks you to map the columns yourself, where you may also choose
*"(no timepoints - compare samples directly)"*.

## Measurement columns

Every remaining **numeric** column is treated as a measurement you can plot and
analyse: `Length`, `Width`, `Area`, `Perimeter`, `Circularity`, `AspectRatio`,
`Solidity`, `Roundness`, and so on. Names are yours to choose - anything numeric
that varies is offered in the Measurement dropdown. Columns that are constant or
look like identifiers are excluded automatically.

---

## Example with a time course

```csv
Condition,Time_min,Replicate,Cell_ID,Length,Width,Area,Perimeter,Circularity,AspectRatio
Uninfected,0,1,1,2.511,1.072,2.028,7.083,0.527,2.342
Uninfected,0,1,2,2.744,0.913,1.912,7.166,0.468,3.005
Phage_MOI_1,60,2,884,6.902,1.221,6.371,16.482,0.295,5.653
```

## Example without timepoints

```csv
Sample,Replicate,Cell_ID,Length,Width,Area,Perimeter,Circularity,AspectRatio
WT,1,1,2.511,1.072,2.028,7.083,0.527,2.342
delta_gp75,1,352,3.284,1.104,2.884,8.955,0.452,2.975
Phage_infected,3,1201,5.118,1.243,4.921,13.166,0.357,4.117
```

---

## Design notes that matter for the statistics

- **Biological replicates are your n.** Three independent experiments with 100
  cells each is far stronger than one experiment with 300 cells. All hypothesis
  tests summarise each replicate to one value first (SuperPlots; Lord et al.
  2020), so cells cannot inflate significance.
- **Two replicates is thin.** With n=2 per group, a Wilcoxon test can never reach
  p<0.05 and a t-test has ~1 degree of freedom. The app warns you and suggests
  the mixed model or the dose-response trend test instead. Three or more is
  much better.
- **Technical replicates (multiple images of the same sample) are not biological
  replicates.** Pool them under the same `Replicate` value.
- **Encode factors in the condition name** if you have them, e.g.
  `pSJExD_gp75_100nM`. The Data tab can split such labels into separate
  vector and numeric dose columns for factorial and dose-response analysis.
- **Timepoints stay numeric** (`0`, `30`, `60`). You can relabel how they are
  *displayed* on figures (e.g. "0 h", "Uninduced") without affecting ordering,
  trajectories, or the trend test.

See [METHODS.md](METHODS.md) for how the figures and statistics are computed.
