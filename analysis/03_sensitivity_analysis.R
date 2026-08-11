## Statistical-analysis sensitivity checks.

set.seed(20260605)

library(readxl)

outdir <- "results/main"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

load(file.path(outdir, "fig3_and_results.RData"))

data <- read_xlsx("data/data.xlsx")
colnames(data) <- c("AG", "GE", "HRS", "HEA", "EHEA",
                    "ES", "JS", "IRS", "SWB", "LS",
                    "FC", "SSES", "SSS", "HS", "DD",
                    "FA", "FMO", "FHEA", "FPP", "MA",
                    "MMO", "MHEA", "MPP", "FZ", "NFI",
                    "REL_GDP", "REL_GDP_Per_Cap")
data <- as.data.frame(data)

X_all <- as.matrix(data[, candidate_names])
Y <- as.numeric(data$SWB)
D <- as.numeric(data$HEA)
treat_index <- which(colnames(X_all) == "HEA")
stab_controls <- setdiff(stab_full$S, treat_index)

extract_lm <- function(fit, term) {
  co <- summary(fit)$coefficients
  beta <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  data.frame(estimate = beta,
             se = se,
             ci_lower = beta - 1.96 * se,
             ci_upper = beta + 1.96 * se,
             p = co[term, grep("Pr\\(", colnames(co))])
}

effect_lm <- function(y, d, controls = NULL, weights = NULL, label) {
  dat <- data.frame(y = y, D = d)
  if (!is.null(controls) && ncol(as.matrix(controls)) > 0) {
    xx <- as.data.frame(controls)
    names(xx) <- paste0("X", seq_len(ncol(xx)))
    dat <- cbind(dat, xx)
  }
  fit <- if (is.null(weights)) lm(y ~ ., data = dat) else lm(y ~ ., data = dat, weights = weights)
  est <- extract_lm(fit, "D")
  data.frame(model = label, est, n = nrow(dat), controls = ncol(dat) - 2)
}

primary <- effect_lm(Y, D, X_all[, stab_controls, drop = FALSE], NULL,
                     "Stab-Knock + OLS")

ps_df <- data.frame(D = D, as.data.frame(X_all[, stab_controls, drop = FALSE]))
names(ps_df) <- c("D", paste0("X", seq_len(ncol(ps_df) - 1)))
ps_fit <- glm(D ~ ., data = ps_df, family = binomial)
ps_hat <- pmin(pmax(predict(ps_fit, type = "response"), 0.01), 0.99)
p_treat <- mean(D == 1)
ipw <- ifelse(D == 1, p_treat / ps_hat, (1 - p_treat) / (1 - ps_hat))
ipw <- pmin(pmax(ipw, quantile(ipw, 0.01)), quantile(ipw, 0.99))

ipw_only <- effect_lm(Y, D, NULL, ipw, "IPW")
ipw_adjusted <- effect_lm(Y, D, X_all[, stab_controls, drop = FALSE], ipw,
                          "IPW + Stab-Knock controls")

## Placebo using false 2014 timing in 2010/2014 waves.
wave_counts <- c(`2010` = 2965, `2014` = 3653, `2018` = 3712, `2020` = 2798)
wave_year <- rep(as.integer(names(wave_counts)), times = as.integer(wave_counts))
pre_idx <- which(wave_year %in% c(2010, 2014))
placebo_df <- data.frame(y = Y[pre_idx],
                         D = D[pre_idx],
                         placebo_post = as.integer(wave_year[pre_idx] == 2014),
                         as.data.frame(X_all[pre_idx, stab_controls, drop = FALSE]))
names(placebo_df) <- c("y", "D", "placebo_post",
                       paste0("X", seq_len(ncol(placebo_df) - 3)))
placebo_fit <- lm(y ~ D * placebo_post + ., data = placebo_df)
placebo_est <- extract_lm(placebo_fit, "D:placebo_post")
placebo <- data.frame(model = "Placebo test, false 2014 timing",
                      placebo_est, n = nrow(placebo_df), controls = length(stab_controls))

## E-value for the updated primary continuous estimate.
evalue_from_continuous <- function(beta, se, y_sd) {
  d <- abs(beta) / y_sd
  d_ci <- max(0, (abs(beta) - 1.96 * se) / y_sd)
  rr <- exp(0.91 * d)
  rr_ci <- exp(0.91 * d_ci)
  evalue <- rr + sqrt(rr * (rr - 1))
  evalue_ci <- ifelse(rr_ci <= 1, 1, rr_ci + sqrt(rr_ci * (rr_ci - 1)))
  data.frame(std_effect_abs = d,
             approximate_RR = rr,
             E_value_point = evalue,
             E_value_CI_limit = evalue_ci)
}
evalue <- evalue_from_continuous(primary$estimate, primary$se, sd(Y))

## Standard knockoff and RML/Lasso SCE are loaded from fig3_and_results.RData.
method_compare <- sce_table[, c("model", "estimate", "se", "ci_lower", "ci_upper", "p", "n", "selected_controls")]
names(method_compare)[names(method_compare) == "selected_controls"] <- "controls"

updated_sensitivity <- rbind(
  primary,
  ipw_only,
  ipw_adjusted,
  method_compare[method_compare$model == "Standard knockoff + OLS",
                 c("model", "estimate", "se", "ci_lower", "ci_upper", "p", "n", "controls")],
  method_compare[method_compare$model == "RML/Lasso + OLS",
                 c("model", "estimate", "se", "ci_lower", "ci_upper", "p", "n", "controls")],
  placebo
)

write.csv(updated_sensitivity,
          file.path(outdir, "updated_statistical_sensitivity.csv"),
          row.names = FALSE)
write.csv(evalue,
          file.path(outdir, "updated_evalue.csv"),
          row.names = FALSE)

save(updated_sensitivity, evalue, ps_fit, ps_hat, ipw, placebo_fit,
     file = file.path(outdir, "updated_statistical_analysis.RData"))

print(updated_sensitivity)
print(evalue)
