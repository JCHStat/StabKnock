## Main method comparisons and post-selection OLS analysis.

set.seed(20260605)

library(readxl)
library(glmnet)
library(knockoff)
library(ggplot2)
library(RColorBrewer)

source("R/stab_knock_functions.R")

outdir <- "results/main"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

data <- read_xlsx("data/data.xlsx")
colnames(data) <- c("AG", "GE", "HRS", "HEA", "EHEA",
                    "ES", "JS", "IRS", "SWB", "LS",
                    "FC", "SSES", "SSS", "HS", "DD",
                    "FA", "FMO", "FHEA", "FPP", "MA",
                    "MMO", "MHEA", "MPP", "FZ", "NFI",
                    "REL_GDP", "REL_GDP_Per_Cap")
data <- as.data.frame(data)

exclude_vars <- c("REL_GDP")
candidate_names <- setdiff(names(data), c("SWB", exclude_vars))
X_all <- as.matrix(data[, candidate_names])
Y <- as.numeric(data$SWB)
D <- as.numeric(data$HEA)
treat_index <- which(colnames(X_all) == "HEA")
n <- nrow(X_all)
p <- ncol(X_all)

methods <- c("Knockoff", "Lasso", "SIS", "Stab-EN", "Stab-AL",
             "Stab-Lasso", "Stab-Knock")
method_colors <- c("Knockoff" = "deepskyblue",
                   "Lasso" = "antiquewhite4",
                   "SIS" = "limegreen",
                   "Stab-EN" = "goldenrod1",
                   "Stab-AL" = "magenta",
                   "Stab-Lasso" = "#984EA3",
                   "Stab-Knock" = "#E41A1C")
line_shapes <- c("Knockoff" = 4, "Lasso" = 17, "SIS" = 8,
                 "Stab-EN" = 9, "Stab-AL" = 6,
                 "Stab-Lasso" = 3, "Stab-Knock" = 12)

pred_aic <- function(y, y_hat, n, d) {
  if (d <= 0 || any(is.na(y_hat))) return(NA_real_)
  RSS <- sum((y - y_hat)^2)
  sigma2_hat <- var(y - y_hat)
  (RSS + 2 * d * sigma2_hat) / (n * sigma2_hat)
}

fit_pred <- function(x_train, y_train, x_new) {
  if (ncol(as.matrix(x_train)) < 1) return(rep(NA_real_, nrow(as.matrix(x_new))))
  cv.fit <- cv.glmnet(x = as.matrix(x_train), y = y_train,
                      type.measure = "mse", nfolds = 10)
  mod <- glmnet(x = as.matrix(x_train), y = y_train,
                alpha = 1, lambda = cv.fit$lambda.min)
  as.numeric(predict(mod, newx = as.matrix(x_new)))
}

select_methods <- function(X_train, y_train, q = 0.2, stab_L = 50) {
  pp <- ncol(X_train)
  knockoffs <- function(X) create.fixed(X, method = "equi")
  knock_res <- knockoff.filter(X = X_train, y = y_train,
                               knockoffs = knockoffs,
                               statistic = stat.lasso_lambdasmax,
                               fdr = q)

  cv.fit <- cv.glmnet(x = X_train, y = y_train, type.measure = "mse", nfolds = 10)
  lambda <- cv.fit$lambda.min
  lasso_mod <- glmnet(x = X_train, y = y_train, alpha = 1, lambda = lambda)
  lasso_sel <- which(as.vector(lasso_mod$beta) != 0)

  sis_sel <- which(abs(cor(X_train, y_train)) > 0.1)

  sp_en <- rep(0, pp)
  sp_al <- rep(0, pp)
  sp_lasso <- rep(0, pp)
  L_stab <- 10
  ada_weight <- 1 / pmax(abs(as.vector(lasso_mod$beta)), 1e-6)
  for (m in seq_len(L_stab)) {
    idx_half <- sample(seq_len(nrow(X_train)), floor(nrow(X_train) / 2), replace = FALSE)
    en_fit <- glmnet(x = X_train[idx_half, , drop = FALSE], y = y_train[idx_half],
                     alpha = 0.5, lambda = lambda)
    en_idx <- which(as.vector(en_fit$beta) != 0)
    sp_en[en_idx] <- sp_en[en_idx] + 1

    al_fit <- glmnet(x = X_train[idx_half, , drop = FALSE], y = y_train[idx_half],
                     alpha = 1, lambda = lambda, penalty.factor = ada_weight)
    al_idx <- which(as.vector(al_fit$beta) != 0)
    sp_al[al_idx] <- sp_al[al_idx] + 1

    sl_fit <- glmnet(x = X_train[idx_half, , drop = FALSE], y = y_train[idx_half],
                     alpha = 1, lambda = lambda)
    sl_idx <- which(as.vector(sl_fit$beta) != 0)
    sp_lasso[sl_idx] <- sp_lasso[sl_idx] + 1
  }

  stab_res <- Knock_SPD(X = X_train, y = y_train, L = stab_L, q = q, cv = 1)

  list(Knockoff = knock_res$selected,
       Lasso = lasso_sel,
       SIS = sis_sel,
       `Stab-EN` = which(sp_en / L_stab > 0.9),
       `Stab-AL` = which(sp_al / L_stab > 0.8),
       `Stab-Lasso` = which(sp_lasso / L_stab > 0.8),
       `Stab-Knock` = stab_res$S)
}

## Fig. 3(b)/(c): model size and AIC over sample size.
N_grid <- c(500, 1000, 5000, 10000)
k_N <- 5
q_primary <- 0.2
ms_N_records <- data.frame()
aic_N_records <- data.frame()

for (rep_id in seq_len(k_N)) {
  for (N in N_grid) {
    idx <- sample(seq_len(n), N, replace = FALSE)
    X_train <- X_all[idx, , drop = FALSE]
    y_train <- Y[idx]
    sels <- select_methods(X_train, y_train, q = q_primary, stab_L = 10)
    for (method in names(sels)) {
      sel <- sels[[method]]
      ms_N_records <- rbind(ms_N_records,
                            data.frame(rep = rep_id, N = N, method = method,
                                       model_size = length(sel)))
      if (length(sel) > 0) {
        pred <- fit_pred(X_train[, sel, drop = FALSE], y_train,
                         X_train[, sel, drop = FALSE])
        aic_val <- pred_aic(y_train, pred, n = N, d = length(sel))
      } else {
        aic_val <- NA_real_
      }
      aic_N_records <- rbind(aic_N_records,
                             data.frame(rep = rep_id, N = N, method = method,
                                        AIC = aic_val))
    }
  }
}

## Fig. 3(b)/(c): model size and AIC over q.
q_grid_ms <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
q_grid_aic <- c(0.2, 0.4, 0.6, 0.8)
N_q <- 10000
k_q_aic <- 5
idx_q <- sample(seq_len(n), N_q, replace = FALSE)
X_q <- X_all[idx_q, , drop = FALSE]
y_q <- Y[idx_q]
ms_q_records <- data.frame()
for (qq in q_grid_ms) {
  sels <- select_methods(X_q, y_q, q = qq, stab_L = 10)
  for (method in names(sels)) {
    ms_q_records <- rbind(ms_q_records,
                          data.frame(q = qq, method = method,
                                     model_size = length(sels[[method]])))
  }
}

aic_q_records <- data.frame()
for (rep_id in seq_len(k_q_aic)) {
  idx <- sample(seq_len(n), N_q, replace = FALSE)
  X_train <- X_all[idx, , drop = FALSE]
  y_train <- Y[idx]
  for (qq in q_grid_aic) {
    sels <- select_methods(X_train, y_train, q = qq, stab_L = 10)
    for (method in names(sels)) {
      sel <- sels[[method]]
      if (length(sel) > 0) {
        pred <- fit_pred(X_train[, sel, drop = FALSE], y_train,
                         X_train[, sel, drop = FALSE])
        aic_val <- pred_aic(y_train, pred, n = N_q, d = length(sel))
      } else {
        aic_val <- NA_real_
      }
      aic_q_records <- rbind(aic_q_records,
                             data.frame(rep = rep_id, q = qq, method = method,
                                        AIC = aic_val))
    }
  }
}

for (nm in c("ms_N_records", "aic_N_records", "ms_q_records", "aic_q_records")) {
  tmp <- get(nm)
  tmp$method <- factor(tmp$method, levels = methods)
  assign(nm, tmp)
}

ms_N_summary <- aggregate(model_size ~ method + N, ms_N_records,
                          function(z) c(mean = mean(z), median = median(z), sd = sd(z)))
ms_N_summary <- do.call(data.frame, ms_N_summary)
names(ms_N_summary) <- c("method", "N", "mean_model_size", "median_model_size", "sd_model_size")
aic_N_summary <- aggregate(AIC ~ method + N, aic_N_records,
                           function(z) c(mean = mean(z, na.rm = TRUE),
                                         median = median(z, na.rm = TRUE),
                                         sd = sd(z, na.rm = TRUE)))
aic_N_summary <- do.call(data.frame, aic_N_summary)
names(aic_N_summary) <- c("method", "N", "mean_AIC", "median_AIC", "sd_AIC")
ms_q_summary <- aggregate(model_size ~ method + q, ms_q_records, mean)
aic_q_summary <- aggregate(AIC ~ method + q, aic_q_records,
                           function(z) c(mean = mean(z, na.rm = TRUE),
                                         median = median(z, na.rm = TRUE),
                                         sd = sd(z, na.rm = TRUE)))
aic_q_summary <- do.call(data.frame, aic_q_summary)
names(aic_q_summary) <- c("method", "q", "mean_AIC", "median_AIC", "sd_AIC")

write.csv(ms_N_records, file.path(outdir, "fig3_MS_N_data.csv"), row.names = FALSE)
write.csv(aic_N_records, file.path(outdir, "fig3_AIC_N_data.csv"), row.names = FALSE)
write.csv(ms_q_records, file.path(outdir, "fig3_MS_q_data.csv"), row.names = FALSE)
write.csv(aic_q_records, file.path(outdir, "fig3_AIC_q_data.csv"), row.names = FALSE)
write.csv(ms_N_summary, file.path(outdir, "fig3_MS_N_summary.csv"), row.names = FALSE)
write.csv(aic_N_summary, file.path(outdir, "fig3_AIC_N_summary.csv"), row.names = FALSE)
write.csv(ms_q_summary, file.path(outdir, "fig3_MS_q_summary.csv"), row.names = FALSE)
write.csv(aic_q_summary, file.path(outdir, "fig3_AIC_q_summary.csv"), row.names = FALSE)

fig_ms_N_line <- ggplot(ms_N_summary,
                        aes(x = N, y = mean_model_size, color = method, shape = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = line_shapes) +
  labs(x = "Sample size", y = "Average MS") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.35, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.text = element_text(size = 8.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, face = "bold"))

fig_ms_q_line <- ggplot(ms_q_summary,
                        aes(x = q, y = model_size, color = method, shape = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = line_shapes) +
  labs(x = "q", y = "MS") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.35, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.text = element_text(size = 8.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, face = "bold"))

fig_aic_N_box <- ggplot(aic_N_records, aes(x = factor(N), y = AIC, fill = method)) +
  geom_boxplot(outlier.size = 0.45, position = position_dodge(width = 0.85),
               linewidth = 0.25, na.rm = TRUE) +
  scale_fill_manual(values = method_colors) +
  labs(x = "Sample size", y = "AIC") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.35, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.text = element_text(size = 8.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, face = "bold"))

fig_aic_q_box <- ggplot(aic_q_records, aes(x = factor(q), y = AIC, fill = method)) +
  geom_boxplot(outlier.size = 0.45, position = position_dodge(width = 0.85),
               linewidth = 0.25, na.rm = TRUE) +
  scale_fill_manual(values = method_colors) +
  labs(x = "q", y = "AIC") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key.width = unit(0.35, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.text = element_text(size = 8.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, face = "bold"))

ggsave(file.path(outdir, "Fig3b_MS_N_line.png"), fig_ms_N_line,
       width = 6, height = 4, dpi = 1000)
ggsave(file.path(outdir, "Fig3b_MS_q_line.png"), fig_ms_q_line,
       width = 6, height = 4, dpi = 1000)
ggsave(file.path(outdir, "Fig3c_AIC_N_boxplot.png"), fig_aic_N_box,
       width = 9, height = 3.2, dpi = 1000)
ggsave(file.path(outdir, "Fig3c_AIC_q_boxplot.png"), fig_aic_q_box,
       width = 9, height = 3.2, dpi = 1000)

## Full-sample post-selection OLS SCE and comparable metrics.
sce_ols <- function(y, d, x, sel, label) {
  controls <- setdiff(sel, treat_index)
  dat <- data.frame(y = y, D = d)
  if (length(controls) > 0) {
    xx <- as.data.frame(x[, controls, drop = FALSE])
    names(xx) <- paste0("X", seq_along(controls))
    dat <- cbind(dat, xx)
  }
  fit <- lm(y ~ ., data = dat)
  co <- summary(fit)$coefficients
  est <- co["D", "Estimate"]
  se <- co["D", "Std. Error"]
  data.frame(model = label, estimate = est, se = se,
             ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se,
             p = co["D", "Pr(>|t|)"], n = nrow(dat),
             selected_controls = length(controls),
             selected_variables = paste(colnames(x)[sel], collapse = ", "))
}

set.seed(20260605)
stab_full <- Knock_SPD(X = X_all, y = Y, L = 50, q = 0.2, cv = 1)
knock_full <- knockoff.filter(X = X_all, y = Y,
                              knockoffs = function(X) create.fixed(X, method = "equi"),
                              statistic = stat.lasso_lambdasmax,
                              fdr = 0.2)
lasso_cv <- cv.glmnet(x = X_all, y = Y, family = "gaussian", alpha = 1, nfolds = 10)
lasso_fit <- glmnet(x = X_all, y = Y, family = "gaussian",
                    alpha = 1, lambda = lasso_cv$lambda.min)
lasso_sel <- which(as.vector(lasso_fit$beta) != 0)

sce_table <- rbind(
  sce_ols(Y, D, X_all, stab_full$S, "Stab-Knock + OLS"),
  sce_ols(Y, D, X_all, knock_full$selected, "Standard knockoff + OLS"),
  sce_ols(Y, D, X_all, lasso_sel, "RML/Lasso + OLS")
)
sce_table$std_effect <- sce_table$estimate / sd(Y)
write.csv(sce_table, file.path(outdir, "post_selection_OLS_SCE.csv"), row.names = FALSE)

selected_stab <- data.frame(
  variable = colnames(X_all)[stab_full$S],
  role = ifelse(colnames(X_all)[stab_full$S] == "HEA", "treatment", "control"),
  W = stab_full$W[stab_full$S],
  S_hat_original = stab_full$S_hat[stab_full$S],
  S_hat_knockoff = stab_full$S_hat[stab_full$S + ncol(X_all)]
)
write.csv(selected_stab, file.path(outdir, "StabKnock_selected_variables.csv"),
          row.names = FALSE)

## Comparable metrics for the Stab-Knock specification.
stab_controls <- setdiff(stab_full$S, treat_index)
stab_control_names <- colnames(X_all)[stab_controls]
stab_primary <- sce_table[sce_table$model == "Stab-Knock + OLS", ]

y_sd <- sd(Y)
y1 <- Y[D == 1]
y0 <- Y[D == 0]
n1 <- length(y1)
n0 <- length(y0)
pooled_sd <- sqrt(((n1 - 1) * var(y1) + (n0 - 1) * var(y0)) / (n1 + n0 - 2))
cohen_d <- (mean(y1) - mean(y0)) / pooled_sd
cohen_d_se <- sqrt((n1 + n0) / (n1 * n0) + cohen_d^2 / (2 * (n1 + n0 - 2)))
cohen_d_ci <- c(cohen_d - 1.96 * cohen_d_se, cohen_d + 1.96 * cohen_d_se)

high_swb_cutoff <- median(Y)
high_swb <- as.integer(Y >= high_swb_cutoff)
or_unadj_fit <- glm(high_swb ~ D, family = binomial)
or_adj_df <- data.frame(high_swb = high_swb, D = D,
                        as.data.frame(X_all[, stab_controls, drop = FALSE]))
names(or_adj_df) <- c("high_swb", "D", paste0("X", seq_len(ncol(or_adj_df) - 2)))
or_adj_fit <- glm(high_swb ~ ., data = or_adj_df, family = binomial)
or_extract <- function(fit) {
  co <- summary(fit)$coefficients
  beta <- co["D", "Estimate"]
  se <- co["D", "Std. Error"]
  c(estimate = exp(beta), ci_lower = exp(beta - 1.96 * se),
    ci_upper = exp(beta + 1.96 * se), p = co["D", "Pr(>|z|)"])
}
or_unadj <- or_extract(or_unadj_fit)
or_adj <- or_extract(or_adj_fit)

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

ps_df <- data.frame(D = D, as.data.frame(X_all[, stab_controls, drop = FALSE]))
names(ps_df) <- c("D", paste0("X", seq_len(ncol(ps_df) - 1)))
ps_fit <- glm(D ~ ., data = ps_df, family = binomial)
ps_hat <- pmin(pmax(predict(ps_fit, type = "response"), 0.01), 0.99)
p_treat <- mean(D == 1)
ipw <- ifelse(D == 1, p_treat / ps_hat, (1 - p_treat) / (1 - ps_hat))
ipw <- pmin(pmax(ipw, quantile(ipw, 0.01)), quantile(ipw, 0.99))
smd_unweighted <- smd_table(D, X_all[, setdiff(seq_len(ncol(X_all)), treat_index), drop = FALSE])
smd_ipw <- smd_table(D, X_all[, setdiff(seq_len(ncol(X_all)), treat_index), drop = FALSE], ipw)
smd_summary <- data.frame(
  specification = c("Before IPW", "After IPW"),
  mean_abs_smd = c(mean(abs(smd_unweighted$smd)), mean(abs(smd_ipw$smd))),
  median_abs_smd = c(median(abs(smd_unweighted$smd)), median(abs(smd_ipw$smd))),
  max_abs_smd = c(max(abs(smd_unweighted$smd)), max(abs(smd_ipw$smd))),
  n_abs_smd_gt_0.10 = c(sum(abs(smd_unweighted$smd) > 0.10),
                        sum(abs(smd_ipw$smd) > 0.10))
)
write.csv(smd_summary, file.path(outdir, "smd_summary.csv"), row.names = FALSE)

fmt_ci <- function(lo, hi) sprintf("[%.3f, %.3f]", lo, hi)
metrics_table <- rbind(
  data.frame(metric = "SCE, Stab-Knock + OLS",
             estimate = stab_primary$estimate, se = stab_primary$se,
             ci_95 = fmt_ci(stab_primary$ci_lower, stab_primary$ci_upper),
             p_value = stab_primary$p,
             details = paste0("Selected controls = ", stab_primary$selected_controls)),
  data.frame(metric = "Standardized SCE",
             estimate = stab_primary$estimate / y_sd, se = stab_primary$se / y_sd,
             ci_95 = fmt_ci(stab_primary$ci_lower / y_sd, stab_primary$ci_upper / y_sd),
             p_value = stab_primary$p,
             details = sprintf("SD(SWB) = %.3f", y_sd)),
  data.frame(metric = "Unadjusted Cohen's d, HEA=1 vs HEA=0",
             estimate = cohen_d, se = cohen_d_se,
             ci_95 = fmt_ci(cohen_d_ci[1], cohen_d_ci[2]),
             p_value = NA_real_,
             details = sprintf("Mean SWB: HEA=1 %.3f; HEA=0 %.3f", mean(y1), mean(y0))),
  data.frame(metric = "Unadjusted OR for high SWB",
             estimate = unname(or_unadj["estimate"]), se = NA_real_,
             ci_95 = fmt_ci(or_unadj["ci_lower"], or_unadj["ci_upper"]),
             p_value = unname(or_unadj["p"]),
             details = paste0("High SWB = SWB >= ", high_swb_cutoff)),
  data.frame(metric = "Adjusted OR for high SWB",
             estimate = unname(or_adj["estimate"]), se = NA_real_,
             ci_95 = fmt_ci(or_adj["ci_lower"], or_adj["ci_upper"]),
             p_value = unname(or_adj["p"]),
             details = paste0("Adjusted for ", length(stab_controls), " selected controls")),
  data.frame(metric = "Mean absolute SMD before IPW",
             estimate = smd_summary$mean_abs_smd[1], se = NA_real_,
             ci_95 = NA_character_, p_value = NA_real_,
             details = sprintf("Median %.3f; max %.3f; |SMD|>0.10: %d",
                               smd_summary$median_abs_smd[1],
                               smd_summary$max_abs_smd[1],
                               smd_summary$n_abs_smd_gt_0.10[1])),
  data.frame(metric = "Mean absolute SMD after IPW",
             estimate = smd_summary$mean_abs_smd[2], se = NA_real_,
             ci_95 = NA_character_, p_value = NA_real_,
             details = sprintf("Median %.3f; max %.3f; |SMD|>0.10: %d",
                               smd_summary$median_abs_smd[2],
                               smd_summary$max_abs_smd[2],
                               smd_summary$n_abs_smd_gt_0.10[2]))
)
write.csv(metrics_table, file.path(outdir, "comparable_metrics.csv"),
          row.names = FALSE)

latex_escape <- function(x) gsub("_", "\\\\_", x)
format_p <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p)))

sce_latex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Post-selection OLS SCE estimates across variable-selection methods}",
  "\\begin{tabular}{lccccc}",
  "\\hline",
  "Method & Estimate & SE & 95\\% CI & $p$ value & Selected controls \\\\",
  "\\hline"
)
for (i in seq_len(nrow(sce_table))) {
  sce_latex <- c(sce_latex, sprintf(
    "%s & %.3f & %.3f & %s & %s & %d \\\\",
    latex_escape(sce_table$model[i]), sce_table$estimate[i], sce_table$se[i],
    fmt_ci(sce_table$ci_lower[i], sce_table$ci_upper[i]),
    format_p(sce_table$p[i]), sce_table$selected_controls[i]))
}
sce_latex <- c(sce_latex, "\\hline", "\\end{tabular}", "\\end{table}")

metrics_latex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Comparable effect metrics for the Stab-Knock model}",
  "\\begin{tabular}{lcccc}",
  "\\hline",
  "Metric & Estimate & SE & 95\\% CI & $p$ value \\\\",
  "\\hline"
)
for (i in seq_len(nrow(metrics_table))) {
  metrics_latex <- c(metrics_latex, sprintf(
    "%s & %.3f & %s & %s & %s \\\\",
    latex_escape(metrics_table$metric[i]), metrics_table$estimate[i],
    ifelse(is.na(metrics_table$se[i]), "--", sprintf("%.3f", metrics_table$se[i])),
    ifelse(is.na(metrics_table$ci_95[i]), "--", metrics_table$ci_95[i]),
    format_p(metrics_table$p_value[i])))
}
metrics_latex <- c(metrics_latex, "\\hline", "\\end{tabular}", "\\end{table}")

ms_latex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Average model size by variable-selection method}",
  "\\begin{tabular}{lcccc}",
  "\\hline",
  "Method & N=500 & N=1000 & N=5000 & N=10000 \\\\",
  "\\hline"
)
for (method in methods) {
  z <- ms_N_summary[as.character(ms_N_summary$method) == method, ]
  z <- z[match(N_grid, z$N), ]
  ms_latex <- c(ms_latex, sprintf(
    "%s & %.2f & %.2f & %.2f & %.2f \\\\",
    method, z$mean_model_size[1], z$mean_model_size[2],
    z$mean_model_size[3], z$mean_model_size[4]))
}
ms_latex <- c(ms_latex, "\\hline", "\\end{tabular}", "\\end{table}")

aic_latex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Median AIC by variable-selection method}",
  "\\begin{tabular}{lcccc}",
  "\\hline",
  "Method & N=500 & N=1000 & N=5000 & N=10000 \\\\",
  "\\hline"
)
for (method in methods) {
  z <- aic_N_summary[as.character(aic_N_summary$method) == method, ]
  z <- z[match(N_grid, z$N), ]
  aic_latex <- c(aic_latex, sprintf(
    "%s & %.3f & %.3f & %.3f & %.3f \\\\",
    method, z$median_AIC[1], z$median_AIC[2],
    z$median_AIC[3], z$median_AIC[4]))
}
aic_latex <- c(aic_latex, "\\hline", "\\end{tabular}", "\\end{table}")

writeLines(c(sce_latex, "", metrics_latex, "", ms_latex, "", aic_latex),
           file.path(outdir, "latex_tables.tex"))

save(exclude_vars, candidate_names, methods, N_grid, q_grid_ms, q_grid_aic,
     ms_N_records, aic_N_records, ms_q_records, aic_q_records,
     ms_N_summary, aic_N_summary, ms_q_summary, aic_q_summary,
     fig_ms_N_line, fig_ms_q_line, fig_aic_N_box, fig_aic_q_box,
     stab_full, knock_full, lasso_cv, lasso_fit, selected_stab,
     sce_table, metrics_table, smd_summary, smd_unweighted, smd_ipw,
     ipw, ps_hat, high_swb_cutoff,
     file = file.path(outdir, "fig3_and_results.RData"))

print(ms_N_summary)
print(aic_N_summary)
print(sce_table[, c("model", "estimate", "se", "ci_lower", "ci_upper", "p", "selected_controls")])
print(metrics_table)
