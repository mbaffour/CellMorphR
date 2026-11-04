# =============================================================================
# CellMorphR analysis template - starting point for scripting your analysis
#
# The Guide > Script Builder tab in the app writes a customised version of this
# for your own data, grouping, test and figure. Generated 2026-08-28
#
# Run it with:   Rscript my_analysis.R
# Keep this file next to CellMorphR_v3.R (it reuses the app's own functions,
# so the numbers and the figure match what you saw in the app exactly).
# =============================================================================

DATA_FILE <- "examples/example_data.csv"   # <- change to your own file
OUT_DIR   <- "cellmorphr_output"

suppressPackageStartupMessages(invisible(source("CellMorphR_v3.R")))
dir.create(OUT_DIR, showWarnings = FALSE)

# ---- 1. Load and validate ---------------------------------------------------
raw <- read.csv(DATA_FILE, check.names = FALSE)
v   <- validate_data(raw)
if (!isTRUE(v$ok)) stop("Data problem: ", v$msg)
d <- v$df
d$Replicate <- as.factor(d$Replicate)
cat("Loaded", nrow(d), "cells,", length(v$measures), "measurements\n")

MEASURES <- c("Length")
GROUP_COLS <- c("Condition", "Time_min")
REFERENCE  <- "Uninfected"   # control / baseline condition

for (MEASURE in MEASURES) {
  cat("\n=== ", MEASURE, " ===\n", sep = "")

  # ---- 2. Replicate-level summaries (this is your n) ------------------------
  rs <- compute_flexible_replicate_summary(d, MEASURE, GROUP_COLS,
                                          yvar = "median_val")

  # ---- 3. Pairwise statistics on replicate summaries ------------------------
  groups_present <- sort(unique(as.character(rs$Statistics_Group)))
  ref_group <- grep(REFERENCE, groups_present, fixed = TRUE, value = TRUE)[1]
  if (is.na(ref_group)) ref_group <- groups_present[1]
  cmp <- build_comparison_table(rs, mode = "vs_reference", reference_group = ref_group)
  res <- run_pairwise_tests(rs, cmp,
                            test = "welch", correction = "holm",
                            reference_group = ref_group, measurement = MEASURE)
  print(res[, c("Group_1", "Group_2", "N_BioRep_Group_1", "N_BioRep_Group_2",
                "Difference_Group2_minus_Group1", "CI95_low", "CI95_high",
                "p_adj", "Significance")])
  write.csv(res, file.path(OUT_DIR, paste0("stats_", MEASURE, ".csv")), row.names = FALSE)

  # ---- 4. SuperPlot figure --------------------------------------------------
  # reference condition first, so it sits at the left of each group
  d$Condition <- factor(d$Condition,
                        levels = c(REFERENCE, setdiff(sort(unique(as.character(d$Condition))), REFERENCE)))
  d$Time_fac <- factor(d$Time_min)
  rs_plot <- compute_replicate_summaries(d, MEASURE)
  rs_plot$Time_fac <- factor(rs_plot$Time_min)
  p <- ggplot(d, aes(x = Time_fac, y = .data[[MEASURE]], fill = Condition)) +
    geom_violin(position = position_dodge(0.8), alpha = 0.7, scale = "width",
                linewidth = 0.3, color = "grey30") +
    geom_point(data = rs_plot,
               aes(x = Time_fac, y = median_val, fill = Condition,
                   group = Condition, shape = factor(Replicate)),
               position = position_dodge(0.8), size = 2.6, colour = "black", stroke = 0.7) +
    scale_shape_manual(values = c(21, 22, 24, 23, 25), name = "Replicate") +
    # Okabe-Ito colourblind-safe palette, the app's default since 3.3.0
    scale_fill_manual(values = rep_len(c("#0072B2","#D55E00","#009E73","#E69F00","#CC79A7","#56B4E9","#F0E442","#999999"), length(unique(d$Condition)))) +
    labs(x = "Time", y = "Length", fill = "Condition",
         title = paste(MEASURE, "by condition"),
         subtitle = "Points = biological replicate medians",
         caption = "Welch's t-test on replicate medians, holm-corrected.") +
    theme_publication(12)

  ggsave(file.path(OUT_DIR, paste0("figure_", MEASURE, ".png")),
         p, width = 9, height = 5.5, dpi = 300)
}

cat("\nDone. Figures and statistics are in ", normalizePath(OUT_DIR), "\n", sep = "")
