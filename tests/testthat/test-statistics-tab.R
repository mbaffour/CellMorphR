# =============================================================================
# End-to-end check of the Statistics tab, driven through shiny::testServer.
#
# raw_data() is a reactiveVal, so a synthetic dataset can be pushed straight in
# and the real observeEvent(input$run_stats) exercised -- including the global
# model branches and the emmeans family (CM-8), which live in the server and
# cannot be reached from the statistical core alone.
# =============================================================================

app_env <- new.env(parent = globalenv())
ok <- tryCatch({
  suppressWarnings(suppressMessages(
    sys.source(file.path(cellmorphr_root(), "CellMorphR_v3.R"), envir = app_env)))
  TRUE
}, error = function(e) { message("could not source the app: ", conditionMessage(e)); FALSE })

test_that("the whole app file sources without error", {
  expect_true(ok)
  expect_true(is.function(app_env$server))
})

run_stats_once <- function(inputs = list(), data = NULL) {
  if (is.null(data)) {
    set.seed(55)
    data <- make_cellmorph_data(ncell = 40, times = c(0, 30, 60))
  }
  data$Condition <- factor(data$Condition, levels = c("Ctrl", "Trt"))
  out <- NULL
  shiny::testServer(app_env$server, {
    raw_data(data)
    measure_cols("Area")
    session$setInputs(
      measure = "Area", stat_measure = "Area", stat_summary = "median_val",
      stat_group_cols = c("Condition", "Time_min"),
      stat_compare_mode = "within_time", stat_test = "welch",
      stat_p_adjust = "holm", data_transform = "none",
      ref_condition = "Ctrl", stat_method = "anova")
    for (nm in names(inputs)) do.call(session$setInputs, inputs[nm])
    session$setInputs(run_stats = 1)
    out <<- stat_results()
  })
  out
}

test_that("the default global model is the repeated-measures ANOVA", {
  skip_if_not(ok)
  r <- run_stats_once()
  expect_equal(r$global_model_kind, "rm_anova")
  expect_true(r$global_model$ok)
  expect_match(r$anova_text, "Repeated-measures ANOVA")
  # Condition tested on dish-to-dish variation: 6 dishes, 2 conditions -> 4 df.
  cond <- r$global_model$table[r$global_model$table$Term == "Condition", ]
  expect_equal(cond$DenDF, 4)
  expect_equal(cond$Error_Stratum, "Rep_unique")
})

test_that("the report states the family size and the transform (CM-8, CM-13)", {
  skip_if_not(ok)
  r <- run_stats_once()
  expect_equal(r$pairwise_family_size, nrow(r$pairwise))
  expect_match(r$pairwise_family, "one family")
  expect_match(r$pairwise_family, as.character(nrow(r$pairwise)))
  expect_equal(r$transform_label, "none (raw values)")

  rlog <- run_stats_once(list(data_transform = "log"))
  expect_equal(rlog$transform_used, "log")
  expect_match(rlog$transform_label, "log(x", fixed = TRUE)
})

test_that("the report names the realised comparison mode, not the request (CM-9)", {
  skip_if_not(ok)
  r <- run_stats_once(list(stat_group_cols = "Condition"))
  expect_equal(r$comparison_mode_requested, "within_time")
  expect_equal(r$comparison_mode, "vs_reference")
  expect_match(r$comparison_mode_message, "needs Time_min")
})

test_that("advisory diagnostics, Levene included, reach the report (CM-12, CM-18)", {
  skip_if_not(ok)
  r <- run_stats_once()
  expect_false(is.null(r$diagnostics))
  txt <- format_assumption_diagnostics(r$diagnostics)
  expect_match(txt, "Brown & Forsythe", fixed = TRUE)
  expect_match(txt, "ADVISORY")
})

test_that("the mixed-model branch corrects over the whole contrast family (CM-8)", {
  skip_if_not(ok && has_lme4 && has_lmerTest && has_emmeans,
              "lme4 + lmerTest + emmeans needed")
  r <- run_stats_once(list(stat_method = "lmm"))
  expect_equal(r$global_model_kind, "lmm")
  expect_true(isTRUE(r$lmm$fitted))
  expect_match(r$lmm_formula, "Rep_unique:Time_fac", fixed = TRUE)

  cons <- as.data.frame(r$pairwise_emm)
  # 3 timepoints x 1 Ctrl-vs-Trt contrast, corrected together rather than one
  # stratum at a time (which with two conditions was no correction at all).
  expect_equal(nrow(cons), 3)
  expect_match(r$emm_family, "one family|as one family")
  expect_match(r$emm_family, "3")
})

test_that("Holm over the emmeans family actually moves the p-values", {
  skip_if_not(ok && has_lme4 && has_lmerTest && has_emmeans,
              "lme4 + lmerTest + emmeans needed")
  set.seed(56)
  d <- make_cellmorph_data(ncell = 60, times = c(0, 30, 60),
                           cond_effect = 0, interaction_effect = 1.5, rep_sd = 0.4)
  d$Condition <- factor(d$Condition, levels = c("Ctrl", "Trt"))
  r <- run_stats_once(list(stat_method = "lmm"), data = d)
  adj <- as.data.frame(r$pairwise_emm)$p_adj
  raw <- as.data.frame(r$pairwise_emm)$p_raw
  expect_equal(adj, p.adjust(raw, "holm"), tolerance = 1e-8)
  expect_true(any(adj > raw + 1e-12))   # the old `by` family did not correct at all
})

test_that("a design too small for a global model says so, and the pairwise survives", {
  skip_if_not(ok)
  set.seed(57)
  d <- make_cellmorph_data(ncell = 20, reps = "R1", times = c(0, 30))
  r <- run_stats_once(list(stat_method = "lmm"), data = d)
  expect_true(is.null(r$model))
  expect_match(r$global_model_message, "Pairwise results above are still valid")
  expect_gt(nrow(r$pairwise), 0)
})

test_that("the Effect Sizes tab survives a condition that is absent at a timepoint (CM-4)", {
  skip_if_not(ok)
  # The ordinary "drug/phage added after the first frame" design: the treatment
  # has no dishes at t = 0. This used to reach `if (NA > 0)` and replace the
  # plot, the table and the download button with a red error block.
  set.seed(58)
  d <- make_cellmorph_data(ncell = 30, times = c(0, 30, 60))
  d <- d[!(d$Condition == "Trt" & d$Time_min == 0), ]
  d$Condition <- factor(d$Condition, levels = c("Ctrl", "Trt"))
  esd <- NULL
  shiny::testServer(app_env$server, {
    raw_data(d); measure_cols("Area")
    session$setInputs(measure = "Area", data_transform = "none",
                      ref_condition = "Ctrl", es_summary = "median_val",
                      es_type = "hedges_g")
    esd <<- es_data()
  })
  expect_false(is.null(esd))
  t0 <- esd[!esd$is_control & esd$Time_min == 0, ]
  expect_equal(nrow(t0), 1)
  expect_true(is.na(t0$hedges_g))
  expect_match(t0$note, "at least 2 biological replicates")
  # Later timepoints are still estimated, with an interval.
  later <- esd[!esd$is_control & esd$Time_min == 60, ]
  expect_true(is.finite(later$hedges_g))
  expect_true(is.finite(later$hedges_g_low) && is.finite(later$hedges_g_high))
})
