# Methodology — How CellMorphR Figures Are Plotted

## Plain-English version

**1. Input data.** Each row of the input file is one segmented cell, with its
morphometric measurements (e.g., Length, Width, Area, Circularity) and three
labels: the experimental condition, the imaging timepoint (minutes), and the
biological replicate it came from.

**2. Two levels of data are drawn on every distribution figure (SuperPlot
style).**
- *Cell level:* the full single-cell distribution for each condition ×
  timepoint group is drawn as a violin (kernel density), jitter/sina points, or
  boxplot. This shows biological spread but is **never** used for p-values.
- *Replicate level:* each biological replicate is collapsed to a single summary
  value (its median, or mean, of the measurement). These per-replicate
  summaries are overlaid as large dots — optionally shape-coded per replicate —
  and the horizontal bar marks the mean of the replicate summaries. These dots
  are the true experimental n.

This two-layer presentation follows the "SuperPlots" recommendation of Lord et
al. (2020): the reader sees both cell-to-cell variability and sample-to-sample
reproducibility in one figure, without any statistics required.

**3. Statistics (when shown).** The unit of replication is the biological
replicate — the dish — everywhere, including in the two global models. In a time
course the *same* dish is imaged at every timepoint, so a design with 3 dishes
per condition and 5 timepoints has 3 independent units per condition, not 15.

- *Pairwise comparisons* are computed on the per-replicate summary values, never
  on individual cells. Welch's t-test (or Wilcoxon rank-sum) with Holm
  correction applied across **every comparison in the table as one family**; the
  output states the family size. Significance brackets on figures use the same
  test, correction and family as the Statistics tab when that tab is left on its
  default within-timepoint configuration, and the figure caption says so when it
  is not.
- *Repeated-measures ANOVA* (the default global model) is fitted on the
  replicate summaries as `y ~ Condition * Time + Error(Replicate)`, with
  `Replicate` uniquely identified per condition. Condition varies between dishes
  and is tested against dish-to-dish variation; Time and Condition × Time vary
  within a dish and are tested against the dish × timepoint variation. The
  output names the error stratum and the denominator degrees of freedom for
  every term. With a single timepoint the design collapses to a one-way ANOVA
  over dishes.
- *Linear mixed model* (optional) is fitted on cells as
  `Measurement ~ Condition * Time + (1|Replicate) + (1|Replicate:Time)`.
  The dish intercept protects the Condition test; the dish × timepoint term is
  what protects Time and Condition × Time, which otherwise fall back on the
  cell-level residual. Denominator degrees of freedom are Satterthwaite
  (`lmerTest`); `emmeans` contrasts use Kenward-Roger and are corrected across
  all Condition-within-timepoint contrasts as one family. Singular fits are
  reported, with the variance component that collapsed named; they are not
  silently replaced by another model.
- *Standardised effect sizes* are Hedges' g on the replicate summaries: the
  difference in means over the (n−1)-weighted pooled SD, multiplied by Hedges'
  small-sample correction J, with a 95% noncentral-t confidence interval. At
  n = 3 per group J = 0.80 and the interval is several units wide; report the
  interval, not the point estimate.
- *Dose-gradient designs* use a linear trend test: the per-replicate summaries
  are regressed against the numeric dose within each stratum (e.g. construct ×
  timepoint), with Holm correction across strata and 95% CIs on the slope.
- *Assumption diagnostics* — Shapiro-Wilk on the residuals, Levene's test with
  median centering (Brown & Forsythe, 1974), and the largest/smallest group SD —
  are computed and reported, and are **advisory only**: CellMorphR does not
  switch test on the strength of them. The test selected is the test that runs.
- *Transforms* are applied once, in the Data tab, and feed every analysis. The
  log transform applies a data-dependent shift when a column contains
  non-positive values; the shift constant is stated in every report, because two
  datasets get different shifts and their effect sizes are therefore not
  comparable.
- With very few replicates (n = 2) the app flags that rank-based tests cannot
  reach significance and that pairwise tests are underpowered; the trend test or
  the mixed model is preferred.

**3b. Trajectory figures.** Each point is one biological replicate's summary
value; the band or error bars show the SD, the SEM, or the 95% t-interval of
those replicate values (the user's choice, named in the figure subtitle
together with the number of replicates per group), never a spread computed
across cells. The trend line is a loess (or, with fewer than four timepoints,
linear) fit through the replicate values and carries no interval of its own.
Conditions are coloured with the Okabe-Ito colourblind-safe palette unless the
user chooses otherwise.

When *normalise to baseline* is on, each replicate's summary is divided by
(fold change) or has subtracted (difference) that same replicate's value at the
earliest timepoint; replicates with no baseline measurement are excluded and
named. The normalisation is always computed on untransformed values, whatever
transform is selected on the Data tab, because a ratio of log values is not a
fold change; the subtitle states this when it applies. Fold changes are summarised by the geometric mean (mean of log2 fold
changes, back-transformed) with SD, SEM or 95% CI computed on the log2 scale
and back-transformed. This normalisation affects the figure only.

**3c. Fraction of cells beyond a threshold.** For each biological replicate
the percentage of cells whose untransformed measurement exceeds (or falls
below) a cut-off is computed. The cut-off is either an absolute value or a
multiple (default 2) of the pooled median of the reference condition's cells at
its earliest timepoint; the exact value is printed on the figure. The
per-replicate percentages are then the unit of analysis: the figure shows each
replicate, the mean across replicates and the SD, SEM or 95% t-interval across
replicates (drawn clipped to 0–100%), and any significance brackets are
Welch's t-tests (or Wilcoxon) on those percentages with the Statistics tab's
correction over one stated family. Cells are never pooled into a single
proportion. With three replicates near 0% or 100% such tests are crude and
should be read descriptively.

**3d. Replicate planning.** For each treatment group at each timepoint the
between-replicate SD of the replicate summaries is pooled with the reference
group's (weighted by degrees of freedom). The number of replicates per group
needed to detect a stated difference at the chosen power with a two-sample
t-test is obtained from `power.t.test`, both for the SD point estimate and for
its upper 95% chi-square bound (conservative). The target difference defaults
to the observed one but should be set to the smallest biologically meaningful
difference. The calculation plans for the pairwise within-timepoint comparison
(the test behind the figure brackets), not for the repeated-measures model, and
it pools the two groups' SDs as if equal although the pairwise test is Welch's;
this is the usual simplification for planning. Values are on the scale the
Statistics tab uses, i.e. after any Data-tab transform. This is prospective
planning; "observed power" of a completed comparison is not reported.

**4. Export.** Figures are rendered with ggplot2 and exported at user-set
dimensions and resolution as PDF, SVG, TIFF, EPS, PNG, or as editable vector
graphics in PowerPoint (officer/rvg), so text and shapes remain editable.

---

## Copy-paste Methods paragraph

> Single-cell morphometric data were visualized and analyzed in CellMorphR
> v3.4.0 (Awuah, M.B.; https://github.com/mbaffour/CellMorphR), an R (v4.6.0)
> Shiny application built on ggplot2. Distributions of single-cell
> measurements were plotted per condition and timepoint as violin (kernel
> density) plots with individual biological-replicate summary values (median
> per replicate) overlaid as points, following the SuperPlot format of Lord et
> al. (2020, J. Cell Biol. 219:e202001064); horizontal bars indicate the mean
> of replicate summaries. To avoid pseudoreplication, all statistical
> comparisons treated the biological replicate (dish) as the experimental unit
> (n = replicates), not the individual cell, including in the global models.
> Where indicated, pairwise comparisons used Welch's two-sample t-tests on
> per-replicate summary values, with Holm correction applied across the whole
> family of comparisons reported (* p < 0.05, ** p < 0.01, *** p < 0.001).
> Global tests of condition and time used a repeated-measures ANOVA on the
> per-replicate summaries with the dish as the error stratum, and/or a linear
> mixed model on cells with random intercepts for dish and for dish × timepoint
> and Satterthwaite degrees of freedom (lme4/lmerTest). Standardised effect
> sizes are Hedges' g on replicate summaries with 95% noncentral-t confidence
> intervals. Dose–response
> relationships were assessed by linear regression of replicate summary values
> against inducer concentration within each construct and timepoint, with
> Holm correction across strata; slope estimates are reported with 95%
> confidence intervals. Figures were exported as vector graphics.

*(If you used the mean instead of the median as the replicate summary, or a
different test/correction in the app, swap those words accordingly — the figure
captions and the downloaded statistics report generated by the app state exactly
which settings, which model and which transform were used. Cite the version you
actually ran: the in-app footer shows it, and the statistics report carries it.)*

> **Version note.** CellMorphR 3.2.0 corrected both global models; before it, the
> two-way ANOVA and the mixed model's within-dish terms were pseudoreplicated,
> and standardised effect sizes were uncorrected Cohen's d. P-values and effect
> sizes produced by v3.1.2 and earlier are not comparable with v3.2.0 output.
> See `CHANGELOG.md` and `AUDIT-2026-09.md`.

## References

- Brown, M.B., Forsythe, A.B. (1974). Robust tests for the equality of
  variances. *Journal of the American Statistical Association* 69(346):364-367.
- Hedges, L.V. (1981). Distribution theory for Glass's estimator of effect size
  and related estimators. *Journal of Educational Statistics* 6(2):107-128.
- Kuznetsova, A., Brockhoff, P.B., Christensen, R.H.B. (2017). lmerTest package:
  tests in linear mixed effects models. *Journal of Statistical Software*
  82(13):1-26.
- Lord, S.J., Velle, K.B., Mullins, R.D., Fritz-Laylin, L.K. (2020).
  SuperPlots: Communicating reproducibility and variability in cell biology.
  *Journal of Cell Biology* 219(6):e202001064. doi:10.1083/jcb.202001064
- Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis.*
  Springer-Verlag New York.
- R Core Team (2025). R: A language and environment for statistical computing.
  R Foundation for Statistical Computing, Vienna, Austria.
