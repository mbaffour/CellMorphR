# Changelog

All notable changes to CellMorphR are recorded here. Versions before 3.2.0 were
not changelogged; for those, see the version table in the
[README](README.md#overview) and the commit history.

## 3.11.0 — 2026-11-10

**Statistics downloads: the whole analysis as one workbook.**

- **Download Analysis Workbook (.xlsx)** on the Statistics tab: every table
  the analysis produced, one sheet each -- Settings, Descriptive, Pairwise,
  Global_model (repeated-measures ANOVA table or mixed-model ANOVA),
  Effect_sizes (partial eta-squared), EMMs, Pairwise_EMM, Random_effects,
  Diagnostics, Test_suggestion, Replicate_summaries, Comparison_preview and
  the Report text. Without the writexl package it is a ZIP of CSVs.
- **Descriptive Stats (CSV)** is a download of its own, next to **Pairwise
  Results (CSV)** and **Full Report (.txt)** (the former Download Results and
  Download Pairwise CSV buttons).
- The same four buttons sit in the sidebar and under the Descriptive, Results
  and Pairwise panels, so they are where the result is being read.
- The **report workbook** (Data, Cell Metrics and Legend tabs) gains
  Descriptive_stats, Global_model, Effect_sizes and EMMs sheets.

## 3.10.0 — 2026-11-09

**Legend key on every figure; Trend, Power and EMM downloads in any format.**

- **Legend key for the secondary figures.** The twelve figures that had none
  now carry the same legend key as the main figures: PCA loadings, morphospace
  trajectory, morphospace diversity and scree; optimal k and morphospace by
  phenotype; raw heatmap, spider/radar and feature importance; dose-response
  trend, estimated marginal means and power curves. Each has a collapsed
  **Legend key** panel under its plot with the on-figure switch, the legend
  draft and **Download Legend Key** (in the tab's export format, size and DPI).
  The key is cut from the live figure, so it follows renamed conditions,
  colours and selection. Figures with nothing mapped to colour, shape or line
  type (scree, optimal k, trend) say so and still get a draft.
- **Legend tab** lists all 21 figures (main first), and **Download All Legend
  Drafts** and the report workbook's Legend_drafts sheet include them.
- **Drafts name the right measurement.** The trend and EMM figures are
  snapshots from when their analysis was run; their drafts name the
  measurement they were run on, not the one currently selected.
- **Trend, Power and EMM figures** now have the format dropdown beside their
  Download button (PDF, PNG, SVG, TIFF, EPS, editable PPTX). The trend figure
  was PDF-only and the power figure had no download.
- **Two-panel figures share one key.** The diversity figure's two panels now
  share one collected legend instead of drawing it twice; with the embed
  switch on, the key sits inside the first panel. The loadings figure's bar
  colour (loading sign) is keyed and legended the same way.
- **Colourblind-safe phenotypes and effects.** Phenotype clusters (composition
  and morphospace-by-phenotype) and the importance figure's ANOVA terms use
  the Okabe-Ito palette instead of ggplot's default hues.
- **Optimal k** says plainly when k was set manually instead of showing nothing.

## 3.9.0 — 2026-11-09

**Downloads: format beside the button, batch export, EPS fix.**

- **File format beside Download Figure.** The format dropdown (PDF, PNG, SVG,
  TIFF, EPS, editable PPTX) now sits right next to each tab's Download button
  instead of at the bottom of the long Export Settings list, where it was
  easy to miss. Same input, so the legend key export follows it as before.
- **Batch export** on the Distributions, Trajectories, Thresholds, Effect
  Sizes and Pub Figure tabs: tick any measurements and any formats and get
  one ZIP with that figure, every current edit applied (titles, labels,
  annotations, typography, size, palette, filters), for each measurement in
  each format. Implemented by a measurement override that every figure
  reactive reads, so the batch uses exactly the on-screen code path; the
  figures return to the selected measurement when it finishes.
- **EPS export wrote nothing on some Macs.** R reported cairo as available
  while its library could not load (no XQuartz), so `ggsave` failed silently.
  The EPS branch now checks that the file exists and falls back to the plain
  PostScript device with a notice.

## 3.8.0 — 2026-11-09

**Legend tab, fuller legend drafts, report workbook, unit picker.**

- **Legend tab** (after Pub Figure): one card per main figure with its legend
  key drawn live, its complete legend draft, and a Download Key button in a
  chosen format; **Download All Legend Drafts (.txt)** writes every draft to
  one file. Key and draft come from the same reactives as the figure tabs, so
  the two can never differ. Figures not yet generated on their own tab say so.
- **Legend drafts now also state** what the marks are (one sentence per plot
  type: violins are kernel-density outlines of single-cell values, each point
  is one biological replicate, hexagons are shaded by cell count, and so on)
  and what data is shown: cells, conditions shown out of the total (named when
  filtered), timepoints shown out of the total, biological replicates per
  group, the source file name, the app version and the date. On Distributions
  the counts follow the sample and timepoint filters, so a filtered figure
  says "2 of 4 conditions shown (Uninfected, Phage_MOI_10)".
- **Report workbook (.xlsx)**, from the Data tab, Cell Metrics tab and Legend
  tab: README (version, date, source, cells, conditions, timepoints,
  replicates per group, transform, units per measurement, reference, test,
  correction, comparison mode, statistical unit), Cell_metrics for **every**
  measurement as one long table, Replicate_summaries, Condition_summaries,
  Pairwise_stats, Thresholds, Power_plan, Colour_map and Legend_drafts.
  Tables not yet computed say so in their sheet. Falls back to a CSV of the
  cell metrics when `writexl` is missing.
- **Unit picker.** The measurement unit is now a dropdown of common
  microscopy units (µm, µm², µm³, nm, mm, px, px², a.u., %, ratio, degrees,
  fold) that also accepts anything typed; still remembered per measurement
  and restored with sessions. Pixel-to-micron conversion of the data itself
  is **not** included: it would rewrite values upstream of every analysis and
  deserves its own, deliberate feature.
- The source file name is now recorded (file, demo, or session) for legends
  and the report, and saved with the session.

## 3.7.0 — 2026-11-09

**The legend text is now a whole legend draft.** The box under each figure's
Legend Key used to describe only the key (which sample is which colour and
marker). It now reads, in order: the measurement with its unit and transform
("Area (µm²), log-transformed"; "raw values" on the Thresholds tab); the key
clauses; then the figure's own subtitle and caption sentences. Those two are
generated from the same inputs as the plot, so the draft states exactly what
is drawn -- replicate medians or means, N per group, SD / SEM / 95% CI,
baseline normalisation, the threshold value, the test and its correction --
and cannot disagree with the figure. Long drafts wrap in the box. Tabs that
plot several measurements at once (Morphospace, Phenotypes, Multi-Feature)
omit the measurement sentence. The composite, whose panels carry no
subtitles, gets one sentence per selected panel ("A, violins = single-cell
distributions; dots = replicate medians; B, ..."), and a trajectory
normalised to baseline says "raw values (normalised to baseline)".

## 3.6.0 — 2026-11-09

**Sample colours you can read, export and carry over.** The Data tab gains a
*Sample Colours (hex codes)* card: every sample with its exact hex code (the
chosen palette, never the colour-vision preview), a coloured swatch and a
plain colour name; a copy-ready `Sample<TAB>Hex` box; **Download Colour Map
(CSV)**; and **Load colour map**, which applies a CSV (`Sample, Hex`, or
`Condition, Colour`; hex codes or R colour names) by matching sample names
first and then the remaining rows in order to the samples that did not match,
so a palette moves to a differently named experiment too. Loading switches the
quick palette to Custom and reports how many samples matched by name and by
order.

- **Saved sessions now restore the per-sample colours.** They were the one
  setting a session did not keep.
- Helpers `normalise_hex()`, `parse_colour_map()` and `assign_colour_map()`
  are in the statistical core and tested (`tests/testthat/test-colour-map.R`).

## 3.5.0 — 2026-11-08

**Legend key on every main figure.** Each figure tab's controls gain a
"Legend Key" block with three things and nothing else: a one-click **Download
Legend Key** (in the tab's chosen format: PDF, PNG, SVG, TIFF, EPS or editable
PPTX), a **legend text starter**, and a **Place legend key on the figure**
switch.

- **The key is the figure's own legend.** It is lifted out of the live plot
  object, so colours, fills, marker shapes, line types, names and font are
  exactly those drawn, and it follows the current sample names, colours and
  selection: rename a condition or hide a sample and the key changes with
  the figure. Figures that hide their legend on purpose (sample names already
  on the x axis) still get a key. The standalone download and the embedded
  version are produced by the same function, so they match exactly.
- **Legend text starter.** One clause per legend, in words:
  "Condition: Uninfected (blue); +N4 (vermillion). Replicate: 1 (filled
  circles); 2 (filled triangles); 3 (filled squares)." Colour words come from
  the nearest of a small named set (Okabe-Ito names included), shape words
  from the plotting-symbol table; continuous scales read as "colour scale from
  white (0) to vermillion (10)". A variable encoded by more than one aesthetic
  (Replicate by cell colour and marker shape) gets one clause with every word
  per entry: "1 (reddish purple, filled circles)". When the colour-vision
  preview is on, the sentence says the words describe simulated colours.
- **Place on the figure.** Moves the legend inside the panel, top-left, on a
  translucent white box. The composite figure keeps its shared legend along
  the bottom instead (a collected legend cannot be placed inside), and its
  control says so.
- **Bundles.** The PPTX panel deck gets one more slide, "Legend key", as
  editable vector graphics. `save_figure()` now draws legend grobs as well as
  plots, including in the PPTX branch.
- Covered: Distributions, Trajectories, Thresholds, Effect Sizes, Morphospace
  scatter, Phenotype composition, Multi-Feature heatmap, Cell Metrics and the
  composite. Not yet: the secondary figures (PCA loadings, scree and
  diversity, raw heatmap, radar, trend, EMM, power).
- The embed switch is part of the saved figure settings and of "Apply to all
  tabs". Tests in `tests/testthat/test-legend-key.R`.

## 3.4.0 — 2026-11-04

**Analysis and data release.** Three new analyses, all built on per-replicate
summaries like everything else in the app, one new way to look at trajectories,
a colour-vision preview, and multi-file import. Nothing in the existing
Statistics, Effect Sizes or pairwise pipelines changes; figures from 3.3.0 are
reproduced exactly unless the new options are switched on.

### Added — analysis

- **Thresholds tab: fraction of cells beyond a cut-off.** The filamentation
  readout. For every dish, the percentage of cells whose raw measurement is
  above (or below) a threshold, defined either as a multiple of the reference
  condition's pooled baseline median (default 2×) or as an absolute value. The
  per-dish percentages are then treated as replicate summaries: mean across
  dishes, SD / SEM / 95% CI across dishes, Welch (or Wilcoxon) tests across
  dishes with the Statistics tab's correction and one stated family. Pooling
  every cell into one proportion is deliberately not offered. Intervals are
  drawn clipped to 0–100% and the subtitle says when clipping happened; the
  About text warns that a t-test on three percentages at 0% or 100% is crude.
  Works with and without a time column. Table and CSV export.
- **Statistics tab > Power: replicate planner.** Plans the *next* experiment.
  For each treatment group the between-replicate SD is pooled with the
  reference group's and the table gives the dishes per group needed to detect
  a chosen difference (default: the observed one) at the target power with a
  two-sample t-test, plus a power-versus-N curve. Because an SD from three
  dishes is itself uncertain, a second, conservative N is planned on the SD's
  upper 95% chi-square bound. Deliberately **not** "observed power" of the
  comparison already run, which is a function of its p-value. Counts above
  1000 are shown as ">1000" (the CSV keeps the number). The panel states its
  scale (the Data-tab transform applies, as on the Statistics tab), that it
  plans for the pairwise within-timepoint test, and that it pools SDs as equal.
- `compute_sig_table()` takes a `unit_label` so its caption can name what was
  tested ("per-replicate percentages of cells above the threshold").

### Added — plotting

- **Trajectories: normalise each replicate to its own baseline.** Fold change
  (÷ its own value at the earliest timepoint) or difference (−), per
  condition × replicate; dishes without a baseline are dropped and named in the
  subtitle. Fold changes are ratios, so their centre is the **geometric mean**
  and the SD / SEM / 95% CI are computed on log2 fold changes and
  back-transformed, giving the asymmetric interval a ratio deserves; the fold
  axis carries no unit; optional log-scale Y. Differences and raw values keep
  the arithmetic summary. The normalisation is always computed on
  **untransformed** values even when a Data-tab transform is active (a ratio of
  logs is not a fold change), and the subtitle says so; non-positive fold
  changes are plotted but excluded from the geometric summary and counted in
  the subtitle. Presentation only: the Statistics and Effect Sizes tabs, and
  panel B of the composite figure, keep showing the raw time course. The CSV
  export records the normalisation, the scale and which spread is shown.
- **Colour-vision preview** on the Data tab: deuteranopia, protanopia,
  tritanopia or greyscale, simulated with the `colorspace` package on the
  condition and replicate palettes. Every figure drawn through the shared
  label pipeline, the composite, and **every exported file** carry a caption
  saying so, and the export shows a warning. Viridis time fills and fixed greys
  are not simulated.

### Added — figure editing and export

- **Journal size presets.** Each tab's Export Settings offer Nature (89 /
  183 mm), Science (2.25 / 4.75 / 7.25 in), Cell (85 / 114 / 174 mm), PLOS full
  width, a 16:9 slide and a poster panel. A preset sets width, height, 300 dpi
  and a journal-style base font (7–8 pt). The base-font slider now goes down
  to 5 pt. The values are common guidance, labelled as such; check the
  journal's current instructions. "Apply to all tabs" now also copies the
  export size, dpi, format and preset.
- **Preview at export size.** A checkbox renders the on-screen figure through
  the same `ggsave` path as the download, at exactly the export width and
  height (96 dpi × a 100–300% zoom), so text size, spacing and aspect on
  screen are what the PDF or PNG will be. Off by default; the normal
  fill-the-panel view is unchanged.
- **Editor extensions.** The double-click editor and the sidebar "Labels &
  caption" panel gain a Y-axis range (every tab except Distributions, which
  keeps its own axis controls and bracket headroom), free-text annotations
  typed as `x | y | text` lines (x may be a number or a category label), and a
  "Reset text to automatic" button. All of it goes through the one shared
  label pipeline, so downloads match the screen.
- **Figure settings survive Save Session**, and can be saved on their own as a
  `.json` figure style and loaded onto another dataset: every tab's text,
  annotations, axis range, typography, grid/border, export size, format and
  transparency. Analysis choices (plot type, filters, thresholds,
  normalisation), colour pickers and the size-preset dropdown are not part of
  a style; after a restore the preset reads "Custom", which is accurate.
- **Transparent background** checkbox for PNG, PDF, SVG and EPS export.
- **Titles, subtitles and captions wrap to the export width.** They used to
  fold at character counts tuned for an 8-inch figure, so at journal column
  width the subtitle ran off the right edge; the fold now follows the export
  width and the font size.

### Added — data

- **Multi-file upload.** Select several CSV/TSV/Excel files (one per image,
  dish or condition); they are stacked with a `Source_File` column (file name
  without extension) that the usual column-mapping dialog lets you assign to
  Condition or Replicate. Nothing is inferred from file names.

### Fixed

- **Empty panels until the Data sidebar was opened.** The measurement and
  reference-condition selectors are rendered inside the Data sidebar, and Shiny
  does not render outputs it considers hidden, so a user who loaded data and
  went straight to another tab saw blank figures until they expanded that
  sidebar once. Both selectors now render regardless.

### Not done

- The Script Builder does not yet emit the Thresholds analysis, the planner or
  the baseline normalisation.
- No zipped "download every figure" and no theme presets; both are easy
  follow-ups if wanted.

## 3.3.0 — 2026-11-04

**Presentation release: figures change, numbers do not.** No statistic, model,
test or correction is touched. Figures regenerated with 3.3.0 will differ from
3.2.0 figures in their default colours and in what the Trajectories band/bars
show, both of which are stated below. Every change lives in the plotting code
and the Data-tab sidebar; the statistical core gains only the extra spread
columns listed under *Added*.

### Changed — figures

- **Conditions default to the Okabe-Ito colourblind-safe palette.** Up to 3.2.0
  a two-condition experiment was drawn in a fixed blue/red pair and anything
  larger in ggplot's rainbow `hue_pal()`, neither of which is distinguishable
  under deuteranopia or in greyscale print. The default is now Okabe-Ito
  (blue, vermillion, green, orange, ...), one distinct colour for up to eight
  conditions and viridis beyond. The old default is still offered as the
  *Blue/Red (pre-3.3 default)* preset, and the Paul Tol set as *Tol bright*.
  Per-condition colour pickers are unchanged.
- **Replicate colours no longer collide with condition colours.** When
  *Color individual cells by replicate* is on, replicates take the Okabe-Ito
  colours the conditions are **not** using (then dark greys), chosen against
  the full condition palette so a replicate keeps its colour whichever
  conditions are shown.
- **Hiding a condition no longer recolours the others.** The Distributions tab
  rebuilt its colours positionally from the filtered condition list, so
  unticking one sample shifted the colours of every remaining one and the tab
  stopped matching the Pub Figure. Colours are now looked up by condition
  name in the one shared palette.
- **The Trajectories band/bars are a stated choice: SD, SEM or 95% CI.** The
  old *Show SEM ribbon* box always drew SEM, except that with the trend line on
  it silently drew loess's own confidence band instead, and the subtitle never
  said which. The sidebar now offers SD / SEM / 95% CI, all computed across
  biological replicates (never cells), the subtitle names the one drawn
  together with the real replicate count per group, and the trend line is
  drawn without a band of its own so one figure never carries two
  differently-defined intervals. Panel B of the composite figure follows the
  same choice. The default is still SEM.
- **Titles no longer repeat the unit.** Axis titles carry the unit (below);
  plot titles name the measurement only.

### Added

- **Measurement unit on the Data tab.** An optional unit (`µm²`, `µm`, `a.u.`)
  is remembered per measurement, saved with the session, and appended to the
  axis titles of the Distributions, Trajectories, Effect Sizes and Pub Figure
  tabs: `Area (µm²)`, `Median Area (µm²) per Replicate`. (The Multi-Feature and
  Cell Metrics tabs plot several measurements at once and stay unit-free.) Under a
  transform the unit moves inside it, `log(Area [µm²])`, because a log axis is
  no longer in square microns. Blank leaves every title exactly as before.
- **Download Plot Data (CSV)** on the Distributions and Trajectories tabs: the
  per-replicate medians and means behind the figure (after the sample and
  timepoint filters), and for Trajectories the group mean/median, SD, SEM,
  95% CI and N as well, for replotting in Prism, Excel or another R session.
- `compute_condition_summaries()` gains `sd_median`, `sd_mean`, `ci95_median`
  and `ci95_mean` (half-widths of the t-based 95% interval on the replicate
  values; `NA` with a single replicate), alongside the existing SEM columns.
  Covered by `tests/testthat/test-plot-helpers.R`, as are the palette helpers.

### Not done (still recommendations in `PUBLICATION_QUALITY_NOTES.md`)

- A colourblind-simulation preview, significance annotations on trajectory
  plots, and per-replicate connectors on SuperPlots.
- The Script Builder and `analysis_template.R` now emit the Okabe-Ito palette,
  but their generated scripts do not carry the unit or the error-bar choice.

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
