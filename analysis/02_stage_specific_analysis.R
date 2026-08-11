## Stage-specific Stab-Knock selection and post-selection OLS SCE estimates.
## Stage-specific Stab-Knock selection and post-selection OLS SCE estimates.
## This is an exploratory robustness analysis; results should be reported as computed.

args <- commandArgs(trailingOnly = TRUE)
analysis_seed <- if (length(args) >= 1) as.integer(args[1]) else 20260605
set.seed(analysis_seed)

library(readxl)
library(glmnet)
library(knockoff)

source("R/stab_knock_functions.R")

outdir <- "results/stage_specific"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

data <- read_xlsx("data/data.xlsx")
colnames(data) <- c("AG", "GE", "HRS", "HEA", "EHEA",
                    "ES", "JS", "IRS", "SWB", "LS",
                    "FC", "SSES", "SSS", "HS", "DD",
                    "FA", "FMO", "FHEA", "FPP", "MA",
                    "MMO", "MHEA", "MPP", "FZ", "NFI",
                    "REL_GDP", "REL_GDP_Per_Cap")
data <- as.data.frame(data)

wave_counts <- c(`2010` = 2965, `2014` = 3653,
                 `2018` = 3712, `2020` = 2798)
data$wave <- rep(as.integer(names(wave_counts)),
                 times = as.integer(wave_counts))

candidate_names <- setdiff(names(data), c("SWB", "REL_GDP", "wave"))

estimate_stage <- function(dat, stage_label, q = 0.2, L = 50) {
  X <- as.matrix(dat[, candidate_names])
  y <- as.numeric(dat$SWB)
  treat_index <- which(colnames(X) == "HEA")

  fit_sel <- Knock_SPD(X = X, y = y, L = L, q = q, cv = 1)
  selected_idx <- fit_sel$S
  selected_names <- colnames(X)[selected_idx]

  controls <- setdiff(selected_names, "HEA")
  controls <- controls[sapply(dat[, controls, drop = FALSE],
                              function(z) length(unique(z)) > 1)]

  df <- data.frame(y = dat$SWB, D = dat$HEA)
  if (length(controls) > 0) {
    xx <- dat[, controls, drop = FALSE]
    names(xx) <- paste0("X", seq_along(controls))
    df <- cbind(df, xx)
  }

  if (length(unique(df$D)) < 2) {
    return(list(
      result = data.frame(
        stage = stage_label, n = nrow(df),
        estimate = NA_real_, se = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        p = NA_real_, selected_n = length(selected_idx),
        control_n = length(controls),
        selected_variables = paste(selected_names, collapse = " "),
        controls = paste(controls, collapse = " ")
      ),
      selection = fit_sel
    ))
  }

  fit_ols <- lm(y ~ ., data = df)
  co <- summary(fit_ols)$coefficients
  beta <- co["D", "Estimate"]
  se <- co["D", "Std. Error"]
  pval <- co["D", grep("Pr\\(", colnames(co))]

  list(
    result = data.frame(
      stage = stage_label, n = nrow(df),
      estimate = beta, se = se,
      ci_lower = beta - 1.96 * se,
      ci_upper = beta + 1.96 * se,
      p = pval, selected_n = length(selected_idx),
      control_n = length(controls),
      selected_variables = paste(selected_names, collapse = " "),
      controls = paste(controls, collapse = " ")
    ),
    selection = fit_sel,
    ols = fit_ols
  )
}

stages <- list(
  `Pre-2016` = data[data$wave %in% c(2010, 2014), , drop = FALSE],
  `Post-2016` = data[data$wave %in% c(2018, 2020), , drop = FALSE],
  `Pre-2020` = data[data$wave %in% c(2010, 2014, 2018), , drop = FALSE],
  `2020` = data[data$wave == 2020, , drop = FALSE]
)

stage_fits <- lapply(names(stages), function(nm) {
  message("Running stage-specific Stab-Knock for ", nm)
  estimate_stage(stages[[nm]], nm, q = 0.2, L = 50)
})
names(stage_fits) <- names(stages)

stage_results <- do.call(rbind, lapply(stage_fits, `[[`, "result"))
stage_results$ci_95 <- sprintf("[%.3f, %.3f]",
                               stage_results$ci_lower,
                               stage_results$ci_upper)

write.csv(stage_results,
          file.path(outdir, paste0("stage_specific_stab_knock_sce_seed_",
                                   analysis_seed, ".csv")),
          row.names = FALSE)

save(stage_results, stage_fits, candidate_names,
     file = file.path(outdir, paste0("stage_specific_stab_knock_sce_seed_",
                                     analysis_seed, ".RData")))

print(stage_results)
