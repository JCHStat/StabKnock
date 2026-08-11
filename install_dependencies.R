required_packages <- c(
  "ggplot2", "glmnet", "hdi", "knockoff", "RColorBrewer", "readxl", "Rfast"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
} else {
  message("All required packages are already installed.")
}
