# =============================================================================
# Presentation helpers added in 3.3.0: the spread columns behind the
# Trajectories band / bars, and the colourblind-safe condition / replicate
# palettes. Pure functions above the END OF STATISTICAL CORE marker.
# =============================================================================

load_cellmorphr_core()

test_that("condition summaries carry SD, SEM and 95% CI across replicates", {
  set.seed(11)
  df <- make_cellmorph_data(cond_effect = 2, times = c(0, 60), ncell = 50)
  rs <- compute_replicate_summaries(df, "Area")
  cs <- compute_condition_summaries(rs)

  expect_true(all(c("sd_median", "sd_mean", "sem_median", "sem_mean",
                    "ci95_median", "ci95_mean", "n_reps") %in% names(cs)))

  one <- rs[rs$Condition == "Trt" & rs$Time_min == 60, ]
  row <- cs[cs$Condition == "Trt" & cs$Time_min == 60, ]
  expect_equal(row$n_reps, 3)
  expect_equal(row$sd_median,   sd(one$median_val))
  expect_equal(row$sem_median,  sd(one$median_val) / sqrt(3))
  expect_equal(row$ci95_median, qt(0.975, 2) * sd(one$median_val) / sqrt(3))
  expect_equal(row$sd_mean,     sd(one$mean_val))
  expect_equal(row$ci95_mean,   qt(0.975, 2) * sd(one$mean_val) / sqrt(3))
  # The CI is wider than the SEM, which is narrower than the SD at N = 3.
  expect_gt(row$ci95_median, row$sem_median)
  expect_lt(row$sem_median, row$sd_median)
})

test_that("a single replicate yields an undefined, not a zero, interval", {
  df <- make_cellmorph_data(times = c(0, 60), reps = "R1", ncell = 30)
  cs <- compute_condition_summaries(compute_replicate_summaries(df, "Area"))
  expect_true(all(cs$n_reps == 1))
  expect_true(all(is.na(cs$sd_median)))
  expect_true(all(is.na(cs$ci95_median)))
  expect_true(all(is.na(cs$ci95_mean)))
})

test_that("error_measure_label names each spread measure", {
  expect_equal(error_measure_label("sem"),  "SEM")
  expect_equal(error_measure_label("sd"),   "SD")
  expect_equal(error_measure_label("ci95"), "95% CI")
})

test_that("conditions default to Okabe-Ito and stay distinct at any count", {
  expect_equal(default_cond_colors(2), c("#0072B2", "#D55E00"))
  expect_equal(default_cond_colors(8), OKABE_ITO)
  for (n in c(1, 3, 8, 12, 20)) {
    cols <- default_cond_colors(n)
    expect_length(cols, n)
    expect_equal(length(unique(cols)), n)
  }
  expect_equal(names(PALETTE_PRESETS)[1], "Okabe-Ito (colourblind-safe)")
  expect_true(all(vapply(PALETTE_PRESETS, length, integer(1)) >= 8))
})

test_that("replicate colours never coincide with the condition colours on the figure", {
  cond_cols <- setNames(default_cond_colors(3), c("A", "B", "C"))
  rep_cols  <- rep_colors_for(c("R1", "R2", "R3"), cond_cols)
  expect_named(rep_cols, c("R1", "R2", "R3"))
  expect_length(intersect(toupper(rep_cols), toupper(cond_cols)), 0)
  expect_equal(length(unique(rep_cols)), 3)
  # Yellow is dropped: it is invisible at cell-point size on a white panel.
  expect_false("#F0E442" %in% toupper(rep_cols))
  # Still works when the conditions have used up the whole palette.
  rep8 <- rep_colors_for(paste0("R", 1:4), OKABE_ITO)
  expect_length(intersect(toupper(rep8), toupper(OKABE_ITO)), 0)
  expect_equal(length(unique(rep8)), 4)
})

test_that("figure size presets are sane and start with Custom", {
  expect_equal(names(FIGURE_SIZE_PRESETS)[1], "Custom")
  expect_null(FIGURE_SIZE_PRESETS[["Custom"]])
  for (nm in names(FIGURE_SIZE_PRESETS)[-1]) {
    ps <- FIGURE_SIZE_PRESETS[[nm]]
    expect_true(all(c("width", "height", "base") %in% names(ps)), info = nm)
    expect_gt(ps$width, 1); expect_lt(ps$width, 20)
    expect_gt(ps$height, 1); expect_lt(ps$height, 12)
    expect_gte(ps$base, 5); expect_lte(ps$base, 26)
  }
  expect_equal(FIGURE_SIZE_PRESETS[["Nature single column (89 mm)"]]$width, 89 / 25.4)
  expect_equal(FIGURE_SIZE_PRESETS[["Science 2 column (7.25 in)"]]$width, 7.25)
})

test_that("parse_annotations reads x | y | text lines and keeps x numeric only when all are", {
  a <- parse_annotations("60 | 12.5 | lysis onset\n90 | 3 | a | b")
  expect_equal(nrow(a), 2)
  expect_type(a$x, "double"); expect_equal(a$x, c(60, 90))
  expect_equal(a$label, c("lysis onset", "a | b"))
  b <- parse_annotations("30 min | 8 | *\n60 | 9 | **")
  expect_type(b$x, "character"); expect_equal(b$x, c("30 min", "60"))
  # malformed lines are skipped; blank input is an empty frame
  # malformed or empty-label lines are skipped
  c <- parse_annotations("no pipes here\n | 5 | \n1 | notanumber | x\n\n2 | 4 | ok")
  expect_equal(nrow(c), 1)
  expect_equal(c$label, "ok"); expect_equal(c$x, 2)
  expect_equal(nrow(parse_annotations("")), 0); expect_equal(nrow(parse_annotations(NULL)), 0)
})
