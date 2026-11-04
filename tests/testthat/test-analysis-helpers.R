# =============================================================================
# Analysis helpers added in 3.4.0: baseline normalisation, threshold fractions,
# the replicate planner, the colour-vision simulator and the significance
# caption's unit label. All pure functions above the END OF STATISTICAL CORE
# marker.
# =============================================================================

load_cellmorphr_core()

test_that("normalise_to_baseline divides or subtracts each dish's own first-timepoint value", {
  rs <- data.frame(Condition = rep(c("A", "B"), each = 4),
                   Replicate = rep(c("R1", "R2"), times = 4),
                   Time_min  = rep(c(0, 0, 60, 60), 2),
                   median_val = c(2, 4, 6, 8, 1, 2, 5, 0),
                   mean_val   = c(2, 4, 6, 8, 1, 2, 5, 0), stringsAsFactors = FALSE)
  f <- normalise_to_baseline(rs, "median_val", "fold")
  expect_equal(f$t0, 0)
  expect_equal(f$rs$median_val[f$rs$Condition == "A" & f$rs$Time_min == 60], c(3, 2))   # 6/2, 8/4
  expect_equal(f$rs$median_val[f$rs$Condition == "B" & f$rs$Time_min == 60], c(5, 0))   # 5/1, 0/2
  expect_true(all(f$rs$median_val[f$rs$Time_min == 0] == 1))
  expect_length(f$dropped, 0)
  d <- normalise_to_baseline(rs, "median_val", "diff")
  expect_equal(d$rs$median_val[d$rs$Condition == "A" & d$rs$Time_min == 60], c(4, 4))
  expect_true(all(d$rs$median_val[d$rs$Time_min == 0] == 0))
  # mean_val is left alone: only the plotted statistic is normalised
  expect_equal(f$rs$mean_val, rs$mean_val)
  none <- normalise_to_baseline(rs, "median_val", "none")
  expect_equal(none$rs$median_val, rs$median_val)
})

test_that("normalise_to_baseline drops and names dishes without a usable baseline", {
  rs <- data.frame(Condition = "A", Replicate = c("R1", "R1", "R2", "R3", "R3"),
                   Time_min  = c(0, 60, 60, 0, 60),
                   median_val = c(2, 4, 9, 0, 3), mean_val = 1, stringsAsFactors = FALSE)
  f <- normalise_to_baseline(rs, "median_val", "fold")
  expect_setequal(f$dropped, c("A / R2", "A / R3"))      # R2: no t0 row; R3: zero baseline
  expect_equal(unique(f$rs$Replicate), "R1")
  d <- normalise_to_baseline(rs, "median_val", "diff")
  expect_equal(d$dropped, "A / R2")                        # a zero baseline is fine for a difference
})

test_that("threshold_value uses the reference condition's baseline median, or an absolute value", {
  df <- data.frame(Condition = rep(c("Ctrl", "Trt"), each = 6),
                   Time_min  = rep(c(0, 0, 0, 60, 60, 60), 2),
                   Area = c(1, 2, 3, 10, 10, 10,  5, 5, 5, 20, 20, 20))
  t <- threshold_value(df, "Area", "multiple", 2, "Ctrl")
  expect_equal(t$base, 2); expect_equal(t$threshold, 4); expect_equal(t$t0, 0); expect_equal(t$ref, "Ctrl")
  t2 <- threshold_value(df, "Area", "multiple", 1.5, "Trt")
  expect_equal(t2$threshold, 7.5)
  # unknown reference falls back to the first condition present
  t3 <- threshold_value(df, "Area", "multiple", 2, "Nope")
  expect_equal(t3$ref, "Ctrl")
  a <- threshold_value(df, "Area", "absolute", 12)
  expect_equal(a$threshold, 12); expect_equal(a$mode, "absolute")
})

test_that("compute_fraction_summaries gives one percentage per dish, never pooled", {
  df <- data.frame(Condition = "A", Time_min = 0,
                   Replicate = rep(c("R1", "R2"), c(4, 2)),
                   Area = c(1, 1, 5, 5,  5, NA))
  f <- compute_fraction_summaries(df, "Area", threshold = 3, direction = "above")
  expect_equal(nrow(f), 2)
  expect_equal(f$frac_val[f$Replicate == "R1"], 50)
  expect_equal(f$n_cells[f$Replicate == "R2"], 1)      # NA is not a cell
  expect_equal(f$frac_val[f$Replicate == "R2"], 100)
  b <- compute_fraction_summaries(df, "Area", threshold = 3, direction = "below")
  expect_equal(b$frac_val[b$Replicate == "R1"], 50)
  expect_equal(b$frac_val[b$Replicate == "R2"], 0)
  # the per-dish fractions feed the ordinary across-replicate summaries
  f$median_val <- f$frac_val; f$mean_val <- f$frac_val
  cs <- compute_condition_summaries(f)
  expect_equal(cs$grand_mean, 75); expect_equal(cs$n_reps, 2)
})

test_that("plan_replicates matches power.t.test and guards degenerate groups", {
  set.seed(3)
  df <- make_cellmorph_data(cond_effect = 2, rep_sd = 1, times = c(0, 60), ncell = 40)
  rs <- compute_replicate_summaries(df, "Area")
  plan <- plan_replicates(rs, "Ctrl", "median_val", power = 0.8, sig_level = 0.05)
  expect_equal(nrow(plan), 2)                              # Trt at 0 and 60
  expect_equal(plan$Condition, c("Trt", "Trt"))
  row <- plan[plan$Time_min == 60, ]
  v_t <- rs$median_val[rs$Condition == "Trt"  & rs$Time_min == 60]
  v_c <- rs$median_val[rs$Condition == "Ctrl" & rs$Time_min == 60]
  expect_equal(row$Observed_Difference, mean(v_t) - mean(v_c))
  expect_equal(row$Between_Replicate_SD, sqrt((var(v_t) * 2 + var(v_c) * 2) / 4))
  expect_equal(row$Target_Difference, abs(row$Observed_Difference))
  expect_equal(row$Replicates_Needed_Per_Group,
               ceiling(power.t.test(delta = row$Target_Difference, sd = row$Between_Replicate_SD, power = 0.8)$n))
  # a user-chosen delta overrides the observed one
  plan2 <- plan_replicates(rs, "Ctrl", "median_val", delta = 5)
  expect_true(all(plan2$Target_Difference == 5))
  # the curve agrees with power.t.test at every N
  cv <- plan_power_curve(plan, n_range = 2:5)
  expect_equal(nrow(cv), 8)
  expect_equal(cv$Power[cv$Time_min == 60 & cv$N == 3],
               power.t.test(n = 3, delta = row$Target_Difference, sd = row$Between_Replicate_SD)$power)
  # degenerate: one dish per group -> no SD, no N; unknown reference -> empty
  one <- rs[rs$Replicate == "R1", ]
  p1 <- plan_replicates(one, "Ctrl")
  expect_true(all(is.na(p1$Between_Replicate_SD))); expect_true(all(is.na(p1$Replicates_Needed_Per_Group)))
  expect_equal(nrow(plan_replicates(rs, "Nope")), 0)
  # zero SD (identical dishes) must not error
  flat <- rs; flat$median_val <- ifelse(flat$Condition == "Trt", 5, 1)
  expect_true(all(is.na(plan_replicates(flat, "Ctrl")$Replicates_Needed_Per_Group)))
})

test_that("simulate_cvd changes the colours, keeps their names, and is a no-op for normal vision", {
  cols <- setNames(default_cond_colors(3), c("A", "B", "C"))
  expect_identical(simulate_cvd(cols, "none"), cols)
  for (mode in c("deutan", "protan", "tritan", "grey")) {
    out <- simulate_cvd(cols, mode)
    expect_named(out, c("A", "B", "C"))
    expect_false(identical(toupper(out), toupper(cols)))
  }
  expect_equal(cvd_label("deutan"), "deuteranopia")
  expect_null(cvd_label("none")); expect_null(cvd_label(NULL))
})

test_that("compute_sig_table names what was tested in its caption", {
  rs <- data.frame(Condition = rep(c("A", "B"), each = 3), Time_min = 0, Replicate = rep(1:3, 2),
                   frac_val = c(10, 12, 11, 40, 45, 42), median_val = 1, stringsAsFactors = FALSE)
  s <- compute_sig_table(rs, list(c("A", "B")), times = 0, sig_yvar = "frac_val",
                         unit_label = "per-replicate percentages of cells above the threshold")
  expect_match(s$caption, "on per-replicate percentages of cells above the threshold")
  expect_equal(nrow(s$annot), 1)
  s2 <- compute_sig_table(rs, list(c("A", "B")), times = 0, sig_yvar = "median_val")
  expect_match(s2$caption, "on replicate medians")
})

test_that("summarise_trajectory uses a geometric centre and log-scale error for fold changes", {
  rs <- data.frame(Condition = "A", Time_min = 60, Replicate = 1:3,
                   median_val = c(2, 4, 8), mean_val = 1, stringsAsFactors = FALSE)
  g <- summarise_trajectory(rs, "median_val", "fold", "sd")
  expect_equal(g$centre, 4)                               # geometric mean of 2, 4, 8
  expect_equal(g$sd_val, sd(log2(c(2, 4, 8))))            # = 1 on the log2 scale
  expect_equal(g$lo, 2); expect_equal(g$hi, 8)            # 2^(2 -/+ 1): asymmetric around 4
  expect_match(g$scale, "geometric")
  a <- summarise_trajectory(rs, "median_val", "none", "sd")
  expect_equal(a$centre, mean(c(2, 4, 8)))
  expect_equal(a$lo, a$centre - sd(c(2, 4, 8))); expect_equal(a$hi, a$centre + sd(c(2, 4, 8)))
  ci <- summarise_trajectory(rs, "median_val", "none", "ci95")
  expect_equal(ci$hi - ci$centre, qt(0.975, 2) * sd(c(2, 4, 8)) / sqrt(3))
  none <- summarise_trajectory(rs, "median_val", "none", NULL)
  expect_true(is.na(none$lo) && is.na(none$hi)); expect_equal(none$centre, a$centre)
  # one dish: centre but no interval
  one <- summarise_trajectory(rs[1, ], "median_val", "fold", "ci95")
  expect_equal(one$centre, 2); expect_true(is.na(one$ci95))
})

test_that("plan_replicates also plans on the SD's upper 95% bound", {
  set.seed(5)
  df <- make_cellmorph_data(cond_effect = 2, rep_sd = 1, times = c(0, 60), ncell = 40)
  rs <- compute_replicate_summaries(df, "Area")
  plan <- plan_replicates(rs, "Ctrl", "median_val")
  row <- plan[plan$Time_min == 60, ]
  expect_equal(row$SD_df, 4)
  expect_equal(row$Between_Replicate_SD_Upper95, row$Between_Replicate_SD * sqrt(4 / qchisq(0.025, 4)))
  expect_gt(row$Between_Replicate_SD_Upper95, row$Between_Replicate_SD)
  expect_gte(row$Replicates_Needed_Conservative, row$Replicates_Needed_Per_Group)
  expect_equal(row$Replicates_Needed_Conservative,
               ceiling(power.t.test(delta = row$Target_Difference, sd = row$Between_Replicate_SD_Upper95, power = 0.8)$n))
})
