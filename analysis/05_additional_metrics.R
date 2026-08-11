## Comparable metrics for reporting analysis results.
## Outputs one table with SCE + 95% CI, standardized effect size,
## Cohen's d, OR, and standardized-difference summaries.

set.seed(20260605)

library(readxl)

outdir <- "results/additional_metrics"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cfps <- read_xlsx("data/data.xlsx")
colnames(cfps) <- c("AG", "GE", "HRS", "HEA", "EHEA",
                    "ES", "JS", "IRS", "SWB", "LS",
                    "FC", "SSES", "SSS", "HS", "DD",
                    "FA", "FMO", "FHEA", "FPP", "MA",
                    "MMO", "MHEA", "MPP", "FZ", "NFI",
                    "REL_GDP", "REL_GDP_Per_Cap")

Y <- as.numeric(cfps$SWB)
D <- as.numeric(cfps$HEA)
X_all <- as.data.frame(cfps[, setdiff(names(cfps), "SWB")])
X_cov <- X_all[, setdiff(names(X_all), "HEA"), drop = FALSE]

selected_path <- file.path("results/main", "StabKnock_selected_variables.csv")
selected_vars <- read.csv(selected_path, stringsAsFactors = FALSE)
selected_controls <- selected_vars$variable[selected_vars$role == "control"]
X_primary_controls <- X_all[, selected_controls, drop = FALSE]

fmt_ci <- function(lo, hi) sprintf("[%.3f, %.3f]", lo, hi)

extract_lm_effect <- function(fit, term = "D") {
  co <- summary(fit)$coefficients
  beta <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  data.frame(
    estimate = beta,
    se = se,
    ci_lower = beta - 1.96 * se,
    ci_upper = beta + 1.96 * se,
    p = co[term, "Pr(>|t|)"]
  )
}

primary_df <- data.frame(y = Y, D = D, X_primary_controls)
names(primary_df) <- c("y", "D", paste0("X", seq_len(ncol(primary_df) - 2)))
primary_fit <- lm(y ~ ., data = primary_df)
primary_sce <- extract_lm_effect(primary_fit, "D")

y_sd <- sd(Y)
std_sce <- primary_sce
std_sce$estimate <- std_sce$estimate / y_sd
std_sce$se <- std_sce$se / y_sd
std_sce$ci_lower <- std_sce$ci_lower / y_sd
std_sce$ci_upper <- std_sce$ci_upper / y_sd

## Unadjusted Cohen's d for raw SWB difference between HEA groups.
y1 <- Y[D == 1]
y0 <- Y[D == 0]
n1 <- length(y1)
n0 <- length(y0)
pooled_sd <- sqrt(((n1 - 1) * var(y1) + (n0 - 1) * var(y0)) / (n1 + n0 - 2))
cohen_d <- (mean(y1) - mean(y0)) / pooled_sd
cohen_d_se <- sqrt((n1 + n0) / (n1 * n0) + cohen_d^2 / (2 * (n1 + n0 - 2)))
cohen_d_ci <- c(cohen_d - 1.96 * cohen_d_se,
                cohen_d + 1.96 * cohen_d_se)

## OR for high SWB. The cutoff follows the primary analysis module:
## high SWB is defined as SWB >= sample median (= 8 in this dataset).
high_swb_cutoff <- median(Y)
high_swb <- as.integer(Y >= high_swb_cutoff)
or_unadj_fit <- glm(high_swb ~ D, family = binomial)
or_adj_df <- data.frame(high_swb = high_swb, D = D, X_primary_controls)
names(or_adj_df) <- c("high_swb", "D", paste0("X", seq_len(ncol(or_adj_df) - 2)))
or_adj_fit <- glm(high_swb ~ ., data = or_adj_df, family = binomial)

extract_or <- function(fit, term = "D") {
  co <- summary(fit)$coefficients
  beta <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  data.frame(
    estimate = exp(beta),
    se = NA_real_,
    ci_lower = exp(beta - 1.96 * se),
    ci_upper = exp(beta + 1.96 * se),
    p = co[term, "Pr(>|z|)"]
  )
}
or_unadj <- extract_or(or_unadj_fit, "D")
or_adj <- extract_or(or_adj_fit, "D")

smd_table <- function(d, controls, weights = NULL) {
  controls_df <- as.data.frame(controls)
  out <- data.frame(variable = names(controls_df), smd = NA_real_)
  for (j in seq_along(controls_df)) {
    x <- as.numeric(controls_df[[j]])
    if (is.null(weights)) {
      m1 <- mean(x[d == 1], na.rm = TRUE)
      m0 <- mean(x[d == 0], na.rm = TRUE)
      v1 <- var(x[d == 1], na.rm = TRUE)
      v0 <- var(x[d == 0], na.rm = TRUE)
    } else {
      w1 <- weights[d == 1]
      w0 <- weights[d == 0]
      x1 <- x[d == 1]
      x0 <- x[d == 0]
      m1 <- sum(w1 * x1) / sum(w1)
      m0 <- sum(w0 * x0) / sum(w0)
      v1 <- sum(w1 * (x1 - m1)^2) / sum(w1)
      v0 <- sum(w0 * (x0 - m0)^2) / sum(w0)
    }
    out$smd[j] <- (m1 - m0) / sqrt((v1 + v0) / 2)
  }
  out
}

ps_df <- data.frame(D = D, X_primary_controls)
names(ps_df) <- c("D", paste0("X", seq_len(ncol(ps_df) - 1)))
ps_fit <- glm(D ~ ., data = ps_df, family = binomial)
ps_hat <- pmin(pmax(predict(ps_fit, type = "response"), 0.01), 0.99)
p_treat <- mean(D == 1)
ipw <- ifelse(D == 1, p_treat / ps_hat, (1 - p_treat) / (1 - ps_hat))
ipw <- pmin(pmax(ipw, quantile(ipw, 0.01)), quantile(ipw, 0.99))

smd_unweighted <- smd_table(D, X_cov)
smd_ipw <- smd_table(D, X_cov, ipw)

smd_summary <- data.frame(
  specification = c("Unweighted HEA group difference",
                    "IPW-weighted HEA group difference"),
  mean_abs_smd = c(mean(abs(smd_unweighted$smd)),
                   mean(abs(smd_ipw$smd))),
  median_abs_smd = c(median(abs(smd_unweighted$smd)),
                     median(abs(smd_ipw$smd))),
  max_abs_smd = c(max(abs(smd_unweighted$smd)),
                  max(abs(smd_ipw$smd))),
  n_abs_smd_gt_0.10 = c(sum(abs(smd_unweighted$smd) > 0.10),
                        sum(abs(smd_ipw$smd) > 0.10)),
  n_abs_smd_gt_0.25 = c(sum(abs(smd_unweighted$smd) > 0.25),
                        sum(abs(smd_ipw$smd) > 0.25))
)

metrics_table <- rbind(
  data.frame(
    metric = "SCE from primary Stab-Knock + OLS",
    estimate = primary_sce$estimate,
    se = primary_sce$se,
    ci_95 = fmt_ci(primary_sce$ci_lower, primary_sce$ci_upper),
    p_value = primary_sce$p,
    details = paste0("n = ", length(Y), "; selected controls = ",
                     ncol(X_primary_controls))
  ),
  data.frame(
    metric = "Standardized SCE (SCE / SD of SWB)",
    estimate = std_sce$estimate,
    se = std_sce$se,
    ci_95 = fmt_ci(std_sce$ci_lower, std_sce$ci_upper),
    p_value = std_sce$p,
    details = sprintf("SD of SWB = %.3f", y_sd)
  ),
  data.frame(
    metric = "Unadjusted Cohen's d for SWB, HEA=1 vs HEA=0",
    estimate = cohen_d,
    se = cohen_d_se,
    ci_95 = fmt_ci(cohen_d_ci[1], cohen_d_ci[2]),
    p_value = NA_real_,
    details = sprintf("Mean SWB: HEA=1 %.3f; HEA=0 %.3f; pooled SD %.3f",
                      mean(y1), mean(y0), pooled_sd)
  ),
  data.frame(
    metric = "Unadjusted OR for high SWB",
    estimate = or_unadj$estimate,
    se = NA_real_,
    ci_95 = fmt_ci(or_unadj$ci_lower, or_unadj$ci_upper),
    p_value = or_unadj$p,
    details = paste0("High SWB defined as SWB >= ", high_swb_cutoff)
  ),
  data.frame(
    metric = "Adjusted OR for high SWB, Stab-Knock controls",
    estimate = or_adj$estimate,
    se = NA_real_,
    ci_95 = fmt_ci(or_adj$ci_lower, or_adj$ci_upper),
    p_value = or_adj$p,
    details = paste0("Adjusted for ", ncol(X_primary_controls),
                     " Stab-Knock-selected controls")
  ),
  data.frame(
    metric = "Standardized differences before IPW",
    estimate = smd_summary$mean_abs_smd[1],
    se = NA_real_,
    ci_95 = NA_character_,
    p_value = NA_real_,
    details = sprintf("Mean |SMD| %.3f; median |SMD| %.3f; max |SMD| %.3f; |SMD|>0.10: %d",
                      smd_summary$mean_abs_smd[1],
                      smd_summary$median_abs_smd[1],
                      smd_summary$max_abs_smd[1],
                      smd_summary$n_abs_smd_gt_0.10[1])
  ),
  data.frame(
    metric = "Standardized differences after IPW",
    estimate = smd_summary$mean_abs_smd[2],
    se = NA_real_,
    ci_95 = NA_character_,
    p_value = NA_real_,
    details = sprintf("Mean |SMD| %.3f; median |SMD| %.3f; max |SMD| %.3f; |SMD|>0.10: %d",
                      smd_summary$mean_abs_smd[2],
                      smd_summary$median_abs_smd[2],
                      smd_summary$max_abs_smd[2],
                      smd_summary$n_abs_smd_gt_0.10[2])
  )
)

write.csv(metrics_table,
          file.path(outdir, "results_comparable_metrics_table.csv"),
          row.names = FALSE)
write.csv(smd_unweighted,
          file.path(outdir, "results_smd_unweighted.csv"),
          row.names = FALSE)
write.csv(smd_ipw,
          file.path(outdir, "results_smd_ipw.csv"),
          row.names = FALSE)
write.csv(smd_summary,
          file.path(outdir, "results_smd_summary.csv"),
          row.names = FALSE)

save(metrics_table, primary_fit, primary_sce, std_sce, cohen_d,
     cohen_d_se, cohen_d_ci, or_unadj_fit, or_adj_fit, or_unadj,
     or_adj, smd_unweighted, smd_ipw, smd_summary, ps_fit, ps_hat,
     ipw, selected_vars, high_swb_cutoff,
     file = file.path(outdir, "results_comparable_metrics_table.RData"))

print(metrics_table)
