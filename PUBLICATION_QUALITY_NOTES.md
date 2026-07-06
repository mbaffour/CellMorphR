# Publication-Quality Notes

Audit of CellMorphR's plotting/presentation against a figure-publication rubric,
plus the safe presentation improvements applied on the `claude/publication-quality`
branch and the deeper changes deferred as recommendations.

> **Scope guard.** These changes touch **presentation only**. The statistical
> methodology — per-replicate summaries, two-way ANOVA on replicate summaries,
> linear mixed model with replicate as a random effect, Welch pairwise tests,
> Holm correction, effect sizes, assumption diagnostics — was **not altered**.
> The tool's anti-pseudoreplication design is its core strength and was left
> untouched. Statistics were **not re-run/verified** in this pass (no R/Python
> stats validation was performed here); only plotting was smoke-tested.

---

## Audit against the publication rubric

| Rubric criterion | Before | Status |
|---|---|---|
| Colourblind-safe palette (Okabe-Ito / viridis) | Custom blue/red for n=2; `scales::hue_pal()` for n>2 (R) — **not colourblind-safe**; fixed blue/red only, no n>2 (Py) | **Fixed** |
| Axis titles **with units** (e.g. µm²) | Bare measure name ("Area"); no unit mechanism | **Fixed (opt-in)** |
| Publication theme (theme_classic / clean mpl) | `theme_publication` (R); `set_pub_style` (Py) | Already good |
| Legible, configurable fonts | Configurable base size both | Already good |
| Error bars explicitly labelled (SD/SEM/CI) | "SEM ribbon" text present | Improved wording |
| n shown | Present but "N=3" **hard-coded** in several subtitles | **Fixed (computed)** |
| Superplot (per-replicate means over cell cloud) | Implemented: replicate medians overlaid on violin/strip/sina/box + dedicated replicate-only dot plot | Already present |
| Vector (PDF/SVG) + ≥300 dpi + configurable size | R: PDF/PNG/SVG/TIFF/EPS, dpi, W/H. Py: pdf/png/svg, dpi, no figsize CLI | Improved (Py) |
| Significance annotation names the test | "Welch's t-test on per-replicate medians…Holm-corrected" caption | Already good |

**Overall:** the tool was already strong on theme, export, superplots, and
test-named significance captions. The real gaps were the **non-colourblind-safe
palette**, the **absence of unit-aware axis titles**, and **hard-coded "N=3"**
text that would be wrong for any experiment without exactly three replicates.

---

## Changes applied (safe, presentation-only)

### `R-code` (Shiny app)

- **Colourblind-safe palettes.** Added an `okabe_ito` vector and a
  `palette_colors(n, palette)` helper (Okabe-Ito default, viridis option); both
  scale to any number of conditions. Replaced the old `cond_colors()` reactive
  (which used non-colourblind-safe `hue_pal()` for n>2 and a fixed blue/red for
  n=2) with a call to `palette_colors()`. Added a **"Colour palette"** selector
  to the Data-tab sidebar.
- **Unit-aware axis titles.** Added an `axis_label(measure, unit, prefix)` helper
  and a **"Measurement unit (optional)"** text input on the Data tab. When set
  (e.g. `µm²`), axis titles across the Distributions, Trajectories, Effect Sizes,
  and Publication-Figure tabs render as `Area (µm²)`, `Median Area (µm²)`, etc.
  Blank leaves titles exactly as before.
- **De-hard-coded replicate count.** Subtitles/annotations now compute the actual
  number of biological replicates per group (`N = k`, or `N = lo–hi` if
  unbalanced) instead of literally printing "N=3".
- **Clearer error-bar wording.** Trajectory subtitle now states
  "band/bars = SEM of replicate summaries" when the SEM layer is shown, and the
  box-plot subtitle now names boxes = IQR.

### `python_script` (CLI)

- **Colourblind-safe palettes.** Added `OKABE_ITO`, `VIRIDIS_ANCHORS`, and a
  `palette_colors(n, palette)` helper; the default `COLORS` constant now maps to
  Okabe-Ito blue/vermillion. Condition colours are assigned via `palette_colors()`
  (one distinct colour per condition, any count). New `--palette {okabe,viridis}`
  flag (default `okabe`).
- **Unit-aware axis titles.** Added `axis_label()` and a `--unit` flag threaded
  through every plotting function (`plot_distributions`, `plot_trajectories`,
  `plot_effect_sizes`, `plot_ecdf`, `plot_composite`). Default (no unit) is
  unchanged.
- **Configurable figure size.** New `--figsize "W,H"` flag for the
  distribution/trajectory plots (defaults preserved when omitted).
- **De-hard-coded replicate count.** Added `_n_rep_text()`; distribution style
  captions and the trajectory caption now report the real replicate N.
- **Clearer error-bar wording.** Trajectory caption now reads
  "ribbon = SEM of replicate summaries".

**Verification of the presentation changes.** `python_script` was
compile-checked and smoke-run on the built-in demo data with the new flags
(`--unit "µm²" --palette viridis --all-plots --show-sig`, plus a dotplot
composite in PNG): all five plot types and the composite figure rendered
without error and the unit appears in the exported axis titles. The R app has no
R runtime available in this environment, so it was validated by static review
and a delimiter/quote-balance check only — **not executed**.

---

## Deferred recommendations (not implemented — need review / runtime testing)

1. **Superplot polish (gold standard for this data).** The superplot pattern is
   already implemented (per-replicate medians overlaid on the cell cloud, plus a
   replicate-only dot plot). Recommended refinements, deferred because they are
   more than cosmetic:
   - Colour/shape each replicate consistently *by replicate ID* across all
     panels (Lord et al., 2020, *J Cell Biol* recommend one colour per replicate
     so the reader can track a replicate's cloud and its mean). The current code
     colours by **condition** and only varies **shape** by replicate in the
     trajectory/dotplot layers.
     - **Reason deferred:** the cell-cloud layers are coloured by condition to
       keep two-condition comparisons legible; switching to per-replicate colour
       changes the visual grammar of every distribution plot and interacts with
       the palette selector — worth a design decision + visual QA, not a blind edit.
   - Draw a connector between each replicate's mean marker and its cloud, and add
     an explicit legend entry "large markers = biological-replicate means".

2. **Units metadata per measurement.** Currently one global unit input applies to
   the selected measure. A per-column unit map (e.g. Area→µm², Perimeter→µm,
   Circularity→a.u.) carried through the Excel/CSV loader would be more robust for
   multi-measure files. Deferred: touches the data model and UI.

3. **Consistent vector-export polish.** `savefig(..., bbox_inches="tight")` is
   applied to some Python plots but not all (trajectory/ecdf/effect omit it);
   standardising would tidy margins. Fonts embedded in PDF/EPS
   (`pdf.fonttype=42`) would guarantee editable text in Illustrator/Inkscape —
   deferred as it changes output file internals and should be checked on a real
   figure.

4. **Colourblind simulation preview.** A one-click "preview under
   deuteranopia/protanopia" toggle in the Shiny app would let users confirm
   accessibility. Non-trivial (needs `colorspace`/`dichromat`); deferred.

5. **Explicit error-bar choice.** Trajectories always use SEM. Offering
   SD / SEM / 95% CI as a user choice (with the caption updating to match) is a
   presentation feature but changes what the ribbon means numerically, so it is
   left as a recommendation rather than a silent default change.

6. **Significance annotations beyond distribution plots.** Brackets/tests are
   shown on distribution plots only. Annotating the trajectory divergence points
   would require deciding *which* test label to attach and is a design/stats
   presentation decision — deferred.
