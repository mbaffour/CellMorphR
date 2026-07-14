#!/usr/bin/env Rscript
# Auto-installs any missing packages on first run, then opens the Shiny app in your browser.
# Launched by the "Run CellMorphR" (.bat / .command) scripts, or run:  Rscript run_app.R
APP  <- "R-code"
PKGS <- c("shiny","bslib","DT","dplyr","tidyr","readxl","writexl","scales","svglite","emmeans","lme4","lmerTest","ggplot2","ggbeeswarm","ggforce","ggridges","patchwork","colourpicker")
options(repos = c(CRAN = "https://cloud.r-project.org"))
missing <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing R packages (first run only): ", paste(missing, collapse = ", "))
  install.packages(missing)
}
app <- source(APP)$value
message("Opening ", "CellMorphR", " in your browser...")
shiny::runApp(app, launch.browser = TRUE)
