# Changelog

All notable changes to CellMorphR are recorded here. Versions before 3.2.0 were
not changelogged; for those, see the version table in the
[README](README.md#overview) and the commit history.

## 3.2.0 — 2026-09-12

**This release changes reported p-values and effect sizes.** Both global models
were pseudoreplicated and are corrected here. Numbers produced by v3.1.2 and
earlier are **not comparable** with 3.2.0 output, and any figure or table whose
p-values came from the Statistics tab's global model, or whose effect sizes came
from the Effect Sizes tab, should be regenerated before publication. The pairwise
path — the one most users run — was correct throughout and is unchanged in
substance.

The release closes the findings of the 2026-09 statistical and software audit;
the full disposition of all 22 findings is in [`AUDIT-2026-09.md`](AUDIT-2026-09.md).

### Corrected — statistics

- **Repeated-measures ANOVA replaces the two-way ANOVA on replicate summaries**
  (CM-1). `aov(y ~ Condition * Time)` counted the same dish once per timepoint,
  so 3 dishes measured at 5 timepoints contributed 15 rows per condition where
  only 3 independent units existed. `fit_rm_anova()` puts the dish in its own
  error stratum, `Error(Rep_unique)`: Condition is tested against dish-to-dish
  variation, Time and Condition × Time against the dish × timepoint variation.
  Measured type-I error for Condition on null data: **41.0% → ~5%**.
- **The mixed model gains a dish × timepoint random term** (CM-2).
  `(1|Rep_unique)` absorbs only what is constant across the whole time course, so
  `Time` and `Condition × Time` were tested against the cell-level residual
  (DenDF in the thousands). The model is now
  `Measurement ~ Condition * Time_fac + (1|Rep_unique) + (1|Rep_unique:Time_fac)`.
  Measured type-I error on null data with dish drift: Time **99.2% → ~5%**,
  Condition × Time **100% → ~5%**. The term is dropped, with an explicit
  statement of what that costs, where it is not estimable.
- **Hedges' g replaces Cohen's d, with a confidence interval** (CM-3). The old
  estimate used `√[(var₁ + var₂)/2]`, which is the pooled SD only at equal n, had
  no small-sample correction (25% upward bias at n = 3 vs 3), and carried no
  interval — while the plot drew "large effect" guides at ±0.8 over the top of
  it. `hedges_g()` uses the (n−1)-weighted pooled SD, Hedges' J, and a 95%
  noncentral-t interval, drawn as an error bar on every point.
- **emmeans contrasts are corrected across the right family** (CM-8).
  `emmeans(model, ~ Condition | Time_fac)` made the timepoint a `by` variable, so
  `adjust = "holm"` corrected within each timepoint — one contrast per stratum
  with two conditions, i.e. no correction at all, printed under a "holm" label.
  Contrasts are now taken unadjusted and corrected over the whole table, with the
  family size stated.
- **One multiplicity family, stated** (CM-6, CM-21). The figure brackets and the
  Statistics tab were two pipelines producing two different p-values for the same
  contrast; at the shipped defaults the table did not even contain the contrast
  the figure drew. The Statistics tab now defaults to "Pairwise within each time
  point", which over the default grouping is exactly the set of contrasts the
  brackets draw. `detect_earliest_divergence()` no longer uses a third,
  different family. Every table names its family and its size.
- **Levene's test implemented** (CM-12). The README has documented it since v3.0;
  no such test existed. `levene_test()` implements Brown–Forsythe (median
  centering) in base R, verified against `car::leveneTest` in the test suite,
  alongside `variance_ratio()` and the existing Shapiro-Wilk.
- **Diagnostics are labelled advisory** (CM-18). Nothing in CellMorphR switches
  test on the strength of a diagnostic. That is a defensible design, and it is
  now said rather than implied.

### Corrected — reporting

- The singular-fit branch announced "Falling back to Two-way ANOVA" and then
  computed no ANOVA at all, so the report promised a global model and showed none
  (CM-7). Singular fits are now reported: which variance component collapsed,
  that the tests are conservative rather than wrong, and which model to trust.
  Designs that cannot support a global model say so and point at the pairwise
  results.
- The comparison mode printed in the report is the mode that was **built**, not
  the one requested (CM-9). "Within each time point" silently degraded to
  "vs reference" whenever its stratum column was not among the grouping columns.
- The transform, including the log shift constant, appears in **every** statistics
  report rather than only on the mixed-model path (CM-13). The log transform's
  data-dependent shift (`log(x − min(x) + 1)`) now travels with the numbers as a
  printable label, because two datasets get different shifts and their effect
  sizes are therefore not comparable.
- The mixed model and the pairwise tests read the same transformed data, so one
  report can no longer hold pairwise tests on log data beside a model on raw data
  under a single "Transform:" line (CM-14). The duplicate LMM-only transform
  selector is gone.
- The on-screen report and the downloaded `.txt` are one function and cannot
  drift apart.
- The Data tab's "N per group" readout reports per-group counts and says
  "unbalanced" when they differ, instead of counting distinct replicate labels
  across the whole file (CM-19).
- `run_anova_for_feature()`'s eta-squared columns are named `partial_eta2`, which
  is what they compute (CM-20).

### Fixed — input handling

- `validate_data()` reports "no numeric measurement columns found" instead of
  dying on `sapply(character(0), …)` in exactly the case that message was written
  for; columns that look numeric but parsed as text are named (CM-5).
- A time column already named `Time_min` is coerced to numeric like any other, so
  `"30 min"` or `"0.5 h"` no longer stays text and sorts lexically (0, 120, 30,
  60). A column that cannot be read as minutes goes to the mapping panel with the
  offending values quoted (CM-17).
- The Effect Sizes tab survives a condition that is absent at a timepoint, or has
  a single replicate, instead of replacing the plot, the table and the download
  button with an R error block (CM-4).
- Dead `vapply(meta_cols, is.null, …)` filter removed (CM-22).

### Added

- **A test suite** (CM-10). The repository previously had no tests of any kind,
  which is how two global models shipped with 41% and ~100% type-I error rates.
  `Rscript tests/run.R` drives the statistical core standalone — everything above
  the `END OF STATISTICAL CORE` marker in `CellMorphR_v3.R` is plain R with no
  Shiny reactivity — across eight files covering replication-unit invariants,
  type-I error simulations for every global-model term, effect sizes, transforms,
  diagnostics, comparison modes, input validation and the Statistics tab.
  `CELLMORPHR_FAST=1` shortens the simulations for a quick local check.

### Documentation

- `METHODS.md` describes the corrected models, states the unit of replication for
  each, and carries a version note. Its copy-paste Methods paragraph is rewritten
  and now cites 3.2.0 and R 4.6.0 (it was stuck at v3.0.0 / R 4.5.0 — CM-15). The
  claim that "all hypothesis tests are computed on per-replicate summary values"
  is qualified: the mixed model is cell-level by construction, and the sentence
  was the strongest claim in the repository and not true as written.
- `README.md` rewrites the statistical methodology sections, removes the "the
  figure and the tables can never disagree" claim, replaces the
  `N − ab = 24 − 8 = 16` residual arithmetic with the per-stratum denominator df,
  and marks the Python sections as documenting a tool not distributed in this
  repository and not revalidated against the corrected models (CM-16).
- `PUBLICATION_QUALITY_NOTES.md` is marked as recommendations rather than applied
  code: none of the items it lists as "Fixed" are present in `CellMorphR_v3.R`,
  including the colourblind-safe palette (CM-11).
