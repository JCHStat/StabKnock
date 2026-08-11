if (!file.exists("DESCRIPTION") || !dir.exists("analysis")) {
  stop("Run this script from the repository root.")
}
if (!file.exists("data/data.xlsx")) {
  stop("Private data file not found. Place it at data/data.xlsx; see data/README.md.")
}

scripts <- c(
  "analysis/01_cfps_analysis.R",
  "analysis/02_stage_specific_analysis.R",
  "analysis/03_sensitivity_analysis.R",
  "analysis/04_mediation_analysis.R",
  "analysis/05_additional_metrics.R"
)

for (script in scripts) {
  message("Running ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (!identical(status, 0L)) stop("Analysis failed: ", script)
}

message("All empirical analyses completed. Outputs are in results/.")
message("The computationally intensive simulation is separate: ",
        "Rscript analysis/06_simulation.R")
