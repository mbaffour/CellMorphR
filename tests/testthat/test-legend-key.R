# =============================================================================
# Legend key (3.5.0): colour / shape / line-type words, the text starter, and
# the legend-grob extraction, on plain ggplots. Pure functions above the core
# marker; the per-tab wiring is exercised headlessly through shiny::testServer.
# =============================================================================

load_cellmorphr_core()
library(ggplot2)

test_that("colour, shape and line-type words are human and cover the defaults", {
  expect_equal(colour_word(c("#0072B2", "#D55E00", "#009E73", "#E69F00")), c("blue", "vermillion", "bluish green", "orange"))
  expect_equal(colour_word("black"), "black")
  expect_equal(colour_word("#FFFFFF"), "white")
  expect_equal(colour_word("grey30"), "dark grey")
  expect_equal(shape_word(c(21, 22, 24, 16, 17)), c("filled circles", "filled squares", "filled triangles", "circles", "triangles"))
  expect_equal(shape_word(99), "shape 99")
  expect_equal(linetype_word(c("dashed", 1)), c("dashed lines", "solid lines"))
})

test_that("legend_text_starter writes one clause per legend from the plot's own scales", {
  d <- data.frame(x = rep(1:3, 2), y = 1:6, Condition = rep(c("Ctrl", "Trt"), each = 3), Replicate = factor(rep(1:3, 2)))
  p <- ggplot(d, aes(x, y, colour = Condition, shape = Replicate)) + geom_point() +
    scale_colour_manual(values = c(Ctrl = "#0072B2", Trt = "#D55E00")) +
    scale_shape_manual(values = c(21, 22, 24))
  txt <- legend_text_starter(p)
  expect_match(txt, "Condition: Ctrl \\(blue\\); Trt \\(vermillion\\)\\.", fixed = FALSE)
  expect_match(txt, "Replicate: 1 \\(filled circles\\); 2 \\(filled squares\\); 3 \\(filled triangles\\)\\.")
  expect_match(legend_text_starter(p, cvd = "deuteranopia"), "simulated for deuteranopia")
  # colour and fill of the same variable collapse to one clause
  p2 <- ggplot(d, aes(x, y, colour = Condition, fill = Condition)) + geom_point(shape = 21) +
    scale_colour_manual(values = c(Ctrl = "#0072B2", Trt = "#D55E00")) +
    scale_fill_manual(values = c(Ctrl = "#0072B2", Trt = "#D55E00"))
  expect_equal(length(gregexpr("Condition:", legend_text_starter(p2))[[1]]), 1)
  expect_equal(legend_text_starter(p2), "Condition: Ctrl (blue); Trt (vermillion).")
  # one variable encoded by colour AND shape: each entry gets both words, once
  p3 <- ggplot(d, aes(x, y, colour = Replicate, shape = Replicate)) + geom_point() +
    scale_colour_manual(values = c("#CC79A7", "#56B4E9", "#999999")) + scale_shape_manual(values = c(21, 22, 24))
  expect_equal(legend_text_starter(p3),
               "Replicate: 1 (reddish purple, filled circles); 2 (sky blue, filled squares); 3 (grey, filled triangles).")
  # continuous fill reads as a colour scale, with its limits
  h <- ggplot(data.frame(a = 1:3, b = 1:3, v = c(0, 5, 10)), aes(a, b, fill = v)) + geom_tile() +
    scale_fill_gradient(low = "white", high = "#D55E00")
  expect_match(legend_text_starter(h), "v: colour scale from white \\(0\\) to vermillion \\(10\\)\\.")
  # nothing mapped -> empty
  expect_equal(legend_text_starter(ggplot(d, aes(x, y)) + geom_point()), "")
})

test_that("legend_key_plot re-enables hidden discrete legends and leaves continuous ones alone", {
  d <- data.frame(x = 1:4, y = 1:4, g = c("a", "a", "b", "b"))
  hidden <- ggplot(d, aes(x, y, fill = g)) + geom_col() + scale_fill_manual(values = c(a = "red", b = "blue"), guide = "none") +
    theme(legend.position = "none")
  expect_null(extract_legend_grob(hidden))
  shown <- legend_key_plot(hidden)
  g <- extract_legend_grob(shown)
  expect_false(is.null(g))
  expect_setequal(intersect(legend_grob_labels(g), c("a", "b")), c("a", "b"))
  # the embedded version carries the same labels as the standalone key
  g_in <- extract_legend_grob(legend_key_plot(hidden, inside = TRUE))
  expect_equal(sort(legend_grob_labels(g_in)), sort(legend_grob_labels(g)))
  # continuous: still a colourbar, not binned
  h <- ggplot(data.frame(a = 1:3, b = 1:3, v = c(0, 5, 10)), aes(a, b, fill = v)) + geom_tile()
  gh <- extract_legend_grob(legend_key_plot(h))
  expect_false(is.null(gh))
  expect_true(any(grepl("colourbar|colorbar", vapply(gh$grobs, function(x) x$name, character(1)), ignore.case = TRUE)) ||
              length(legend_grob_labels(gh)) >= 3)
  sz <- legend_key_size(g)
  expect_true(sz$w > 0.5 && sz$w < 5 && sz$h > 0.3 && sz$h < 5)
})

test_that("legend_draft_text assembles measurement, key, subtitle and caption into sentences", {
  expect_equal(tidy_sentence("violins = cell distributions;\ndots = replicate means"), "Violins = cell distributions; dots = replicate means.")
  expect_equal(tidy_sentence("Already a sentence."), "Already a sentence.")
  expect_length(tidy_sentence(NULL), 0); expect_length(tidy_sentence("  "), 0)
  d <- data.frame(x = rep(1:3, 2), y = 1:6, Condition = rep(c("Ctrl", "Trt"), each = 3))
  p <- ggplot(d, aes(x, y, colour = Condition)) + geom_point() +
    scale_colour_manual(values = c(Ctrl = "#0072B2", Trt = "#D55E00")) +
    labs(subtitle = "each point = one replicate median (N = 3 per group); bars = SEM of replicate medians",
         caption = "Welch's t-test on replicate medians, Holm-corrected.")
  txt <- legend_draft_text(p, measurement = "Area (um2), log-transformed")
  expect_equal(txt, paste("Measurement: Area (um2), log-transformed.",
                          "Condition: Ctrl (blue); Trt (vermillion).",
                          "Each point = one replicate median (N = 3 per group); bars = SEM of replicate medians.",
                          "Welch's t-test on replicate medians, Holm-corrected."))
  # colour-vision note is added once, and not again when the caption already carries it
  expect_match(legend_draft_text(p, cvd = "deuteranopia"), "simulated for deuteranopia")
  p2 <- p + labs(caption = "Colours simulated for deuteranopia (preview, not the chosen palette).")
  expect_equal(length(gregexpr("simulated for", legend_draft_text(p2, cvd = "deuteranopia"))[[1]]), 1)
  # no measurement, no subtitle: just the key
  expect_equal(legend_draft_text(ggplot(d, aes(x, y, colour = Condition)) + geom_point() +
                                   scale_colour_manual(values = c(Ctrl = "#0072B2", Trt = "#D55E00"))),
               "Condition: Ctrl (blue); Trt (vermillion).")
})

test_that("plot descriptions and the data sentence read correctly", {
  expect_match(plot_description("dist", "violin"), "^Violins are kernel-density")
  expect_match(plot_description("cm", "hexbin"), "Hexagons")
  expect_match(plot_description("es", "hedges_g"), "Hedges' g")
  expect_match(plot_description("traj"), "one biological replicate")
  expect_equal(plot_description("traj", "no_such_type"), plot_description("traj"))   # falls back to the prefix entry
  expect_null(plot_description("dist", "no_such_type"))                              # no prefix-level entry for dist
  expect_null(plot_description("nope"))
  s <- legend_data_sentence(4558, conds_shown = c("A", "C"), conds_total = c("A", "B", "C"),
                            times_shown = c("0 min", "30 min"), times_total = c("0 min", "30 min", "60 min"),
                            reps_range = c(3, 3), source = "example_data.csv", version = "3.8.0", date = "2026-11-09")
  expect_equal(s, "Data: 4,558 cells; 2 of 3 conditions shown (A, C); 2 of 3 timepoints shown (0 min, 30 min); 3 biological replicates per group; source: example_data.csv; CellMorphR 3.8.0, 2026-11-09.")
  s2 <- legend_data_sentence(100, conds_shown = c("WT", "Mut"), conds_total = c("WT", "Mut"), reps_range = c(2, 3))
  expect_equal(s2, "Data: 100 cells; 2 conditions (WT, Mut); 2-3 biological replicates per group.")
  expect_true("µm²" %in% UNIT_CHOICES)
})

# 3.10.0: every figure in the Legend-tab registry has a description, and the
# secondary figures (loadings, scree, radar, trend, EMM, power ...) are in it.
test_that("every registered figure has a legend description, including the secondary figures", {
  expect_true(all(c("pca_loadings", "pca_traj", "pca_diversity", "pca_scree", "clust_k", "clust_pca",
                    "mf_raw", "mf_radar", "mf_importance", "trend", "emm", "power") %in% names(LEGEND_FIGURE_TITLES)))
  expect_equal(anyDuplicated(names(LEGEND_FIGURE_TITLES)), 0)
  # dist, es and cm are described per chart type, so pass one for them
  typed <- c(dist = "violin", es = "abs_diff", cm = "bar")
  for (pf in names(LEGEND_FIGURE_TITLES)) expect_false(is.null(plot_description(pf, typed[pf])), info = pf)
  expect_match(plot_description("power"), "replicates per group")
  # figures drawing several measurements never name one in the draft
  expect_true(all(c("pca_scree", "mf_radar", "clust_pca") %in% LEGEND_MULTI_MEASURE))
  expect_false(any(c("trend", "emm", "power") %in% LEGEND_MULTI_MEASURE))
})
