## Mediation analysis with standardized covariates.
## This version standardizes the covariate matrix X before mediation models.
## HEA is retained as the binary treatment, and the current mediator is kept
## on its original scale to preserve the mediation interpretation.

set.seed(20260617)

library(readxl)

outdir <- "results/mediation"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

data <- as.data.frame(read_xlsx("data/data.xlsx", sheet = 1))
colnames(data) <- c("AG", "GE", "HRS", "HEA", "EHEA",
                    "ES", "JS", "IRS", "SWB", "LS",
                    "FC", "SSES", "SSS", "HS", "DD",
                    "FA", "FMO", "FHEA", "FPP", "MA",
                    "MMO", "MHEA", "MPP", "FZ", "NFI",
                    "REL_GDP", "REL_GDP_Per_Cap")

## Controls selected by the Stab-Knock main model.
selected_controls <- c("GE", "HRS", "JS", "IRS", "LS", "FC", "DD",
                       "FA", "FMO", "FHEA", "FPP", "MA", "MMO",
                       "MHEA", "MPP", "FZ", "NFI")

## Standardize the X matrix used for covariate adjustment. We keep original
## HEA and mediator variables separately, and use z_* variables for controls.
X <- data[, setdiff(names(data), "SWB"), drop = FALSE]
X_scaled <- as.data.frame(scale(X))
names(X_scaled) <- paste0("z_", names(X_scaled))
dat_scaled <- cbind(data, X_scaled)

fit_mediation_once <- function(dat, mediator, mediator_family = c("linear", "binomial")) {
  mediator_family <- match.arg(mediator_family)

  controls <- setdiff(selected_controls, mediator)
  z_controls <- paste0("z_", controls)
  names(z_controls) <- controls
  z_controls <- z_controls[sapply(dat[, z_controls, drop = FALSE],
                              function(x) length(unique(x[!is.na(x)])) > 1)]
  control_terms <- unname(z_controls)

  d_m <- dat[, c("HEA", mediator, control_terms), drop = FALSE]
  d_y <- dat[, c("SWB", "HEA", mediator, control_terms), drop = FALSE]
  names(d_m)[names(d_m) == mediator] <- "M"
  names(d_y)[names(d_y) == mediator] <- "M"

  m_formula <- as.formula(paste("M ~ HEA",
                                if (length(control_terms) > 0)
                                  paste("+", paste(control_terms, collapse = " + "))
                                else ""))
  y_formula <- as.formula(paste("SWB ~ HEA + M",
                                if (length(control_terms) > 0)
                                  paste("+", paste(control_terms, collapse = " + "))
                                else ""))

  if (mediator_family == "binomial") {
    m_fit <- glm(m_formula, data = d_m, family = binomial())
    d1 <- d_m
    d0 <- d_m
    d1$HEA <- 1
    d0$HEA <- 0
    m_diff <- mean(predict(m_fit, newdata = d1, type = "response") -
                     predict(m_fit, newdata = d0, type = "response"))
  } else {
    m_fit <- lm(m_formula, data = d_m)
    m_diff <- coef(m_fit)["HEA"]
  }

  y_fit <- lm(y_formula, data = d_y)
  b <- coef(y_fit)["M"]
  ade <- coef(y_fit)["HEA"]
  acme <- as.numeric(m_diff * b)
  total <- as.numeric(ade + acme)
  prop <- as.numeric(acme / total)
  c(ACME = acme, ADE = as.numeric(ade),
    Total_Effect = total, Prop_Mediated = prop)
}

bootstrap_mediation <- function(dat, mediator, mediator_family, B = 1000) {
  point <- fit_mediation_once(dat, mediator, mediator_family)
  boot <- matrix(NA_real_, nrow = B, ncol = length(point))
  colnames(boot) <- names(point)
  n <- nrow(dat)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    boot[b, ] <- tryCatch(
      fit_mediation_once(dat[idx, , drop = FALSE], mediator, mediator_family),
      error = function(e) rep(NA_real_, length(point))
    )
  }
  out <- data.frame(
    Mediator = mediator,
    Effect = names(point),
    Estimate = as.numeric(point),
    LowerCI = apply(boot, 2, quantile, probs = 0.025, na.rm = TRUE),
    UpperCI = apply(boot, 2, quantile, probs = 0.975, na.rm = TRUE),
    Boot_valid = colSums(!is.na(boot)),
    stringsAsFactors = FALSE
  )
  list(summary = out, boot = boot, point = point)
}

res_ehea <- bootstrap_mediation(dat_scaled, "EHEA", "binomial", B = 1000)
res_irs <- bootstrap_mediation(dat_scaled, "IRS", "linear", B = 1000)

mediation_results_scaled_X <- rbind(res_ehea$summary, res_irs$summary)

write.csv(mediation_results_scaled_X,
          file.path(outdir, "updated_mediation_results_scaled_X.csv"),
          row.names = FALSE)

save(mediation_results_scaled_X, res_ehea, res_irs, selected_controls,
     file = file.path(outdir, "updated_mediation_results_scaled_X.RData"))

print(mediation_results_scaled_X)
