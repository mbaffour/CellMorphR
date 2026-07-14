#!/usr/bin/env Rscript
# Auto-installs any missing packages on first run, then opens the Shiny app in your browser.
# Launched by the "Run CellMorphR" (.bat / .command) scripts, or run:  Rscript run_app.R
APP  <- "R-code"
PKGS <- c("shiny","bslib","DT","dplyr","tidyr","readxl","writexl","scales","svglite","emmeans","lme4","lmerTest","ggplot2","ggbeeswarm","ggforce","ggridges","patchwork","colourpicker")
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Use a writable personal library — R in "Program Files" usually isn't user-writable,
# so first-run installs would otherwise fail without admin.
ul <- Sys.getenv("R_LIBS_USER")
if (nzchar(ul)) {
  ul <- strsplit(ul, .Platform$path.sep, fixed = TRUE)[[1]][1]
} else {
  mm <- strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1]
  ul <- file.path(path.expand("~"), "R", "library", paste0(R.version$major, ".", mm))
}
if (!dir.exists(ul)) dir.create(ul, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(ul, .libPaths()))

missing <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing R packages into ", ul, " (first run only): ", paste(missing, collapse = ", "))
  install.packages(missing, lib = ul)
}
app <- source(APP)$value
message("Opening CellMorphR in your browser...")
shiny::runApp(app, launch.browser = TRUE)
