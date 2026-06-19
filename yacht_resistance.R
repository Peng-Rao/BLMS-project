rm(list=ls())
# Yacht hydrodynamics — UCI #243
yh = read.csv("data/yacht_hydro.csv")
dim(yh)
head(yh)
summary(yh)

# The power-law in Fn is best seen on log-log axes;
# resistance > 0 everywhere (min = 0.01), so plain log() needs no offset
plot(log(yh$Froude), log(yh$resistance),
     pch = 16, cex = 0.5,
     xlab = "log(Froude)", ylab = "log(resistance)",
     main = "Yacht hydrodynamic resistance — power-law scaling")
abline(lm(log(resistance) ~ log(Froude), data = yh),
       col = "red", lwd = 2)

set.seed(2026)
yh$log_res <- log(yh$resistance)   # no offset: min(resistance) = 0.01 > 0
yh$log_Fr  <- log(yh$Froude)

fit <- lm(log_res ~ log_Fr + LCB + prismatic + L_disp + beam_draught + L_beam, data = yh)

## reference prior => credible intervals = OLS confidence intervals
round(cbind(estimate = coef(fit), confint(fit)), 4)
cat("R^2 =", round(summary(fit)$r.squared, 4),
    "  residual sd =", round(summary(fit)$sigma, 4), "\n")

cf <- summary(fit)$coefficients
bF <- cf["log_Fr", "Estimate"]; se <- cf["log_Fr", "Std. Error"]; df <- fit$df.residual
ci <- confint(fit)["log_Fr", ]
## beta_F | y ~ t_df centred at bF with scale se  =>  P(beta_F > 4)
cat(sprintf("beta_F: mean = %.3f, 95%% CI = [%.3f, %.3f], P(beta_F > 4) = %.3f\n",
            bF, ci[1], ci[2], pt((bF - 4) / se, df = df)))

## log-log scatter with posterior-mean line and 95% credible band
##   (= OLS confidence band under the reference prior)
hull_means <- colMeans(yh[, c("LCB","prismatic","L_disp","beam_draught","L_beam")])
xg   <- seq(min(yh$log_Fr), max(yh$log_Fr), length.out = 120)
band <- predict(fit, data.frame(log_Fr = xg, as.list(hull_means)), interval = "confidence")
band_col <- rgb(.85, .1, .1, .30)
plot(yh$log_Fr, yh$log_res, pch = 16, cex = 0.6, col = rgb(0,0,0,0.45),
     xlab = "log(Froude)", ylab = "log(resistance)", main = "Power-law fit")
polygon(c(xg, rev(xg)), c(band[, "lwr"], rev(band[, "upr"])), col = band_col, border = NA)
lines(xg, band[, "fit"], col = "firebrick", lwd = 2.4)
legend("topleft", bty = "n",
       legend = c("data", "posterior mean fit", "95% credible band"),
       pch = c(16, NA, 15), lty = c(NA, 1, NA), lwd = c(NA, 2.4, NA),
       pt.cex = c(1, NA, 1.8), col = c(rgb(0,0,0,0.45), "firebrick", band_col))

## posterior of beta_F: Student-t centred at the OLS estimate
xb <- seq(bF - 4 * se, bF + 4 * se, length.out = 200)
plot(xb, dt((xb - bF) / se, df) / se, type = "l", col = "firebrick", lwd = 2.4,
     main = expression("Posterior of "*beta[F]),
     xlab = expression(beta[F]), ylab = "posterior density")
abline(v = ci, col = "firebrick", lty = 2)
abline(v = 4, col = "navy", lwd = 2)
legend("topright", bty = "n", lwd = c(2.4, 1, 2), lty = c(1, 2, 1),
       col = c("firebrick", "firebrick", "navy"),
       legend = c("posterior", "95% CI", expression(beta[F] == 4)))

op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
plot(fitted(fit), resid(fit), pch = 16, cex = 0.5, col = rgb(0,0,0,0.45),
     xlab = "fitted log(resistance)", ylab = "residual", main = "Residuals vs fitted")
abline(h = 0, col = "firebrick"); lines(lowess(fitted(fit), resid(fit)), col = "navy", lwd = 2)
qqnorm(resid(fit), pch = 16, cex = 0.5, col = rgb(0,0,0,0.45), main = "Normal Q-Q")
qqline(resid(fit), col = "firebrick")
par(op)

## Full conditionals of the independent (semi-conjugate) prior:
##   beta | sigma^2, y    ~ N( V (S0^{-1} m0 + X'y / sigma^2), V ),
##                          with V = (S0^{-1} + X'X / sigma^2)^{-1}
##   sigma^2 | beta, y    ~ Inv-Gamma( a0 + n/2, b0 + ||y - X beta||^2 / 2 )
set.seed(2026)
X <- model.matrix(fit); y <- yh$log_res
n <- nrow(X); p <- ncol(X)
fr  <- which(colnames(X) == "log_Fr")
XtX <- crossprod(X); Xty <- crossprod(X, y)

gibbs_betaF <- function(tau, n_iter = 10000, burn = 2000) {
  m0  <- rep(0, p); m0[fr] <- 4            # physics prior centred at 4 on beta_F only
  psd <- rep(100, p); psd[fr] <- tau       # vague-but-proper N(0, 100^2) elsewhere
  S0_inv <- diag(1 / psd^2)
  a0 <- 1e-3; b0 <- 1e-3
  beta <- coef(fit); sig2 <- summary(fit)$sigma^2
  keep <- numeric(n_iter - burn)
  for (it in seq_len(n_iter)) {
    V    <- chol2inv(chol(S0_inv + XtX / sig2))
    mu   <- V %*% (S0_inv %*% m0 + Xty / sig2)
    beta <- as.vector(mu + t(chol(V)) %*% rnorm(p))
    sig2 <- 1 / rgamma(1, a0 + n / 2, b0 + sum((y - X %*% beta)^2) / 2)
    if (it > burn) keep[it - burn] <- beta[fr]
  }
  keep
}

taus  <- c(0.05, 0.1, 0.2, 0.5)
draws <- lapply(taus, gibbs_betaF)

tab <- t(sapply(draws, function(d)
  c(mean = mean(d), q2.5 = unname(quantile(d, .025)),
    q97.5 = unname(quantile(d, .975)), P_gt4 = mean(d > 4))))
rownames(tab) <- sprintf("tau = %.2f", taus)
## reference row: closed-form Student-t posterior (= OLS confint)
ref <- c(bF, ci[1], ci[2], 1 - pt((4 - bF) / se, df))
round(rbind(tab, reference = ref), 3)

## crude ACF-based effective sample size, as an MCMC sanity check
ess <- sapply(draws, function(d) {
  ac <- acf(d, lag.max = 50, plot = FALSE)$acf[-1]
  length(d) / (1 + 2 * sum(ac[ac > 0]))
})
cat("approx ESS for beta_F:", round(ess), "\n")

draw_overlay <- function() {
  cols <- c("#d7301f", "#fc8d59", "#f1a340", "#2b8cbe")
  xb <- seq(3.9, 4.95, length.out = 400)
  plot(xb, dt((xb - bF) / se, df) / se, type = "l", lwd = 2.6, col = "grey30",
       main = expression("Posterior of " * beta[F] * " vs prior strength"),
       xlab = expression(beta[F]), ylab = "posterior density",
       ylim = c(0, max(sapply(draws, function(d) max(density(d)$y)),
                       dt(0, df) / se)))
  for (i in seq_along(taus)) lines(density(draws[[i]]), col = cols[i], lwd = 2)
  abline(v = 4, col = "navy", lty = 2, lwd = 2)
  legend("topleft", bty = "n", lwd = 2, col = c("grey30", cols),
         legend = c(expression("reference " * pi(sigma^2) %prop% sigma^-2),
                    sprintf("N(4, %.2f²)", taus)))
}
draw_overlay()
## keep the figure used by the report in sync with the notebook
pdf("figures/posterior_betaF_prior.pdf", width = 7, height = 4.4)
draw_overlay(); invisible(dev.off())

if (!requireNamespace("BAS", quietly = TRUE)) {
  stop("Package 'BAS' is required. Install it with install.packages('BAS') before knitting.")
}
library(BAS)

# Optional hull-geometry covariates. These are the five variables requested in the project.
hull_covariates <- c("LCB", "prismatic", "L_disp", "beam_draught", "L_beam")

# Sanity check: fail loudly if a requested covariate is missing or misspelled.
stopifnot(all(c("log_res", "log_Fr", hull_covariates) %in% names(yh)))

# In the model-selection step, log_Fr must be present in every model.
# The five hull covariates are optional, so bas.lm enumerates 2^5 = 32 submodels.
bas_fit <- bas.lm(
  log_res ~ log_Fr + LCB + prismatic + L_disp + beam_draught + L_beam,
  data = yh,
  prior = "BIC",
  modelprior = uniform(),
  method = "deterministic",
  n.models = 2^length(hull_covariates),
  include.always = ~ log_Fr,
  force.heredity = FALSE
)

# Check that the model space has really been fully enumerated.
n_models_expected <- 2^length(hull_covariates)
n_models_observed <- length(bas_fit$postprobs)
cat("Number of models enumerated:", n_models_observed, "\n")
stopifnot(n_models_observed == n_models_expected)

# BAS may store the selected variables either as a list of indices or as a matrix,
# depending on the package version. The following code converts both cases into
# a stable 0/1 model-inclusion matrix.
term_names <- c("(Intercept)", "log_Fr", hull_covariates)

make_which_matrix <- function(which_obj, term_names) {
  if (is.list(which_obj) && !is.matrix(which_obj)) {
    out <- t(sapply(which_obj, function(idx) {
      # BAS indexes terms from 0: 0 = intercept, 1 = log_Fr, 2:6 = hull covariates.
      as.integer((seq_along(term_names) - 1) %in% idx)
    }))
    colnames(out) <- term_names
    return(out)
  }

  out <- as.matrix(which_obj)
  if (ncol(out) == length(term_names)) {
    colnames(out) <- term_names
  } else if (ncol(out) == length(term_names) - 1) {
    colnames(out) <- term_names[-1]
  } else {
    stop("Unexpected format for bas_fit$which. Inspect str(bas_fit$which).")
  }
  out
}

model_matrix <- as.data.frame(make_which_matrix(bas_fit$which, term_names))

# Normalise posterior model probabilities and attach them to the model matrix.
post_prob <- bas_fit$postprobs / sum(bas_fit$postprobs)
model_matrix$post_prob <- post_prob

# Check that log_Fr is indeed included in all models and that the five hull terms exist.
stopifnot(all(model_matrix$log_Fr == 1))
stopifnot(all(hull_covariates %in% names(model_matrix)))

# Posterior inclusion probabilities computed directly from the enumerated models.
# This avoids relying on version-dependent names in bas_fit$probne0.
pip_hull <- sapply(hull_covariates, function(v) {
  sum(model_matrix[[v]] * model_matrix$post_prob)
})

pip_sorted <- sort(pip_hull, decreasing = TRUE)
pip_table <- data.frame(
  covariate = names(pip_sorted),
  PIP = round(as.numeric(pip_sorted), 3),
  evidence = ifelse(pip_sorted >= 0.75, "strong",
                    ifelse(pip_sorted >= 0.50, "moderate", "weak")),
  row.names = NULL
)

pip_table

# Explicit answer to the project question:
# "which of LCB, prismatic, L/disp, B/draught, L/beam carry signal once Froude is in the model?"
label_map <- c(
  LCB = "LCB",
  prismatic = "prismatic",
  L_disp = "L/disp",
  beam_draught = "B/draught",
  L_beam = "L/beam"
)

signal_table <- data.frame(
  requested_name = unname(label_map[pip_table$covariate]),
  R_variable = pip_table$covariate,
  PIP = pip_table$PIP,
  evidence = pip_table$evidence,
  answer = ifelse(
    pip_table$evidence == "moderate",
    "carries moderate signal",
    "weak evidence; not retained as carrying signal"
  ),
  row.names = NULL
)

signal_table

best_covariate <- pip_table$covariate[1]
best_requested_name <- unname(label_map[best_covariate])
best_pip <- pip_table$PIP[1]

cat("Largest PIP:", best_requested_name, "(", best_covariate, ") =", best_pip, "\n")
cat(
  "Answer:",
  best_requested_name, "(", best_covariate, ") is the only hull-geometry covariate with moderate evidence of carrying signal once Froude is included; the remaining hull covariates show weak evidence.\n"
)

barplot(
  rev(pip_sorted),
  horiz = TRUE,
  las = 1,
  xlim = c(0, 1),
  xlab = "Posterior inclusion probability",
  main = "BAS/BMA posterior inclusion probabilities"
)
abline(v = 0.5, lty = 2)
text(
  x = rev(pip_sorted) + 0.04,
  y = seq_along(pip_sorted),
  labels = round(rev(pip_sorted), 3),
  cex = 0.8
)

top_models <- model_matrix[
  order(-model_matrix$post_prob),
  c(hull_covariates, "post_prob")
]

top_models$model <- apply(
  top_models[, hull_covariates, drop = FALSE],
  1,
  function(z) {
    chosen <- hull_covariates[as.logical(z)]
    if (length(chosen) == 0) {
      "log_Fr only"
    } else {
      paste("log_Fr +", paste(chosen, collapse = " + "))
    }
  }
)

top_models <- top_models[, c("model", "post_prob")]
head(transform(top_models, post_prob = round(post_prob, 3)), 10)

bma_coef <- coef(bas_fit)
print(bma_coef)

fit_bas_prior <- function(prior_name, alpha = NULL) {
  args <- list(
    formula = log_res ~ log_Fr + LCB + prismatic + L_disp + beam_draught + L_beam,
    data = yh,
    prior = prior_name,
    modelprior = uniform(),
    method = "deterministic",
    n.models = n_models_expected,
    include.always = ~ log_Fr,
    force.heredity = FALSE
  )
  if (!is.null(alpha)) {
    args$alpha <- alpha
  }
  do.call(bas.lm, args)
}

extract_pip_table <- function(fit_obj, prior_label) {
  mm <- as.data.frame(make_which_matrix(fit_obj$which, term_names))
  pp <- fit_obj$postprobs / sum(fit_obj$postprobs)
  mm$post_prob <- pp

  # Same checks as in the baseline analysis.
  stopifnot(length(pp) == nrow(mm))
  stopifnot(length(pp) == n_models_expected)
  stopifnot(all(mm$log_Fr == 1))

  pip <- sapply(hull_covariates, function(v) {
    sum(mm[[v]] * mm$post_prob)
  })

  out <- data.frame(
    prior = prior_label,
    requested_name = unname(label_map[names(pip)]),
    R_variable = names(pip),
    PIP = round(as.numeric(pip), 3),
    evidence = ifelse(pip >= 0.75, "strong",
                      ifelse(pip >= 0.50, "moderate", "weak")),
    row.names = NULL
  )
  out[order(-out$PIP), ]
}

bas_fits_sensitivity <- list(
  "BIC / reference" = bas_fit,
  "g-prior (g = n)" = fit_bas_prior("g-prior", alpha = nrow(yh)),
  "JZS" = fit_bas_prior("JZS")
)

sensitivity_table <- do.call(
  rbind,
  lapply(names(bas_fits_sensitivity), function(pr) {
    extract_pip_table(bas_fits_sensitivity[[pr]], pr)
  })
)
rownames(sensitivity_table) <- NULL

# Full sensitivity table: PIP for every hull covariate under every prior.
sensitivity_table

# Most supported hull covariate under each prior.
top_by_prior <- do.call(
  rbind,
  lapply(split(sensitivity_table, sensitivity_table$prior), function(d) {
    d[which.max(d$PIP), ]
  })
)
rownames(top_by_prior) <- NULL

top_by_prior

# Variables that reach at least moderate evidence under every prior.
moderate_sets <- lapply(split(sensitivity_table, sensitivity_table$prior), function(d) {
  d$R_variable[d$PIP >= 0.50]
})
robust_moderate <- Reduce(intersect, moderate_sets)
robust_moderate_names <- unname(label_map[robust_moderate])

if (length(robust_moderate) == 0) {
  sensitivity_answer <- paste(
    "No hull covariate reaches PIP >= 0.50 under all priors.",
    "The safest conclusion is that hull effects are weak and prior-sensitive once Froude is included."
  )
} else {
  sensitivity_answer <- paste(
    paste(robust_moderate_names, collapse = ", "),
    "is/are the hull covariate(s) with at least moderate evidence under all priors.",
    "The remaining hull covariates have weak or prior-sensitive support."
  )
}

cat("Sensitivity answer:", sensitivity_answer, "\n")

y_res <- yh$resistance
n_obs <- length(y_res)

## centre covariates at a y^2-weighted mean: removes the strong posterior
## correlation between b0 and gamma_F that ruins mixing with plain centring
w_info <- y_res^2 / sum(y_res^2)
mF <- sum(w_info * log(yh$Froude))
mB <- sum(w_info * yh$beam_draught)
xF_c <- log(yh$Froude)  - mF
xB_c <- yh$beam_draught - mB

X0 <- cbind(intercept = 1, log_Fr = xF_c)                    # M0: Froude only
X1 <- cbind(X0, beam_draught = xB_c)                         # M1: + beam_draught

save_current_pdf <- function(filename, width, height) {
  dir.create("figures", showWarnings = FALSE)
  dev.copy(pdf, file = file.path("figures", filename),
           width = width, height = height)
  invisible(dev.off())
}

if (!requireNamespace("rjags", quietly = TRUE)) {
  stop("Package 'rjags' and the JAGS program are required. Install JAGS (https://mcmc-jags.sourceforge.io), then install.packages('rjags').")
}
library(rjags)
library(coda)

## Heteroscedastic nonlinear model, written for JAGS.
##   y_i ~ Gamma(shape_i, rate_i) moment-matched to mean mu_i and variance sigma_i^2,
##   mu_i = exp(X_i %*% theta),  sigma_i = sigma0 * mu_i^delta,
##   shape_i = mu_i^2 / sigma_i^2,  rate_i = mu_i / sigma_i^2.
## The Gamma keeps the response on its positive support (resistance > 0).
## Nothing here is conjugate (the exp link and the mean-dependent variance break
## conjugacy), so JAGS automatically assigns its univariate SLICE sampler
## (base::RealSlicer) to every parameter -- verified below with list.samplers().
## theta is kept as a vector node so the monitored columns are theta[1..p],
## matching the rest of the analysis.
jags_hetero <- "
model {
  for (i in 1:n) {
    mu[i]    <- exp(inprod(X[i, ], theta))
    sigma[i] <- sigma0 * pow(mu[i], delta)
    shape[i] <- (mu[i] * mu[i]) / (sigma[i] * sigma[i])
    rate[i]  <- mu[i] / (sigma[i] * sigma[i])
    y[i] ~ dgamma(shape[i], rate[i])
  }
  for (j in 1:p) {
    theta[j] ~ dnorm(prior_mean[j], 1 / (prior_sd[j] * prior_sd[j]))
  }
  logs0  ~ dnorm(0, 1 / 100)   # log sigma0 ~ N(0, 10^2)
  sigma0 <- exp(logs0)
  delta  ~ dnorm(0, 1)         # delta ~ N(0, 1^2);  0 = homoscedastic, 1 = multiplicative
}"

## delta pinned to 0 -> constant-variance (homoscedastic) baseline for comparison.
jags_homo <- "
model {
  for (i in 1:n) {
    mu[i]    <- exp(inprod(X[i, ], theta))
    shape[i] <- (mu[i] * mu[i]) / (sigma0 * sigma0)
    rate[i]  <- mu[i] / (sigma0 * sigma0)
    y[i] ~ dgamma(shape[i], rate[i])
  }
  for (j in 1:p) {
    theta[j] ~ dnorm(prior_mean[j], 1 / (prior_sd[j] * prior_sd[j]))
  }
  logs0  ~ dnorm(0, 1 / 100)
  sigma0 <- exp(logs0)
  delta  <- 0
}"

## Fit the model with JAGS. Chains are over-dispersed in gamma_F (and delta) so
## that Gelman-Rubin is a genuine convergence test; reproducible via .RNG.seed.
run_jags <- function(X, prior_mean, prior_sd, homoscedastic = FALSE,
                     n_chains = 4, n_adapt = 1000, n_burn = 2000,
                     n_iter = 8000, thin = 2,
                     gamma_inits = c(3, 6, 5, 4), delta_inits = c(0, 1, 0.5, 1.2),
                     seed0 = 2026) {
  p <- ncol(X)
  stopifnot(length(prior_mean) == p, length(prior_sd) == p)
  dat <- list(y = y_res, X = X, n = nrow(X), p = p,
              prior_mean = prior_mean, prior_sd = prior_sd)
  theta0 <- coef(lm(log(y_res) ~ X - 1))           # crude start from the log fit

  inits <- lapply(seq_len(n_chains), function(ch) {
    th <- theta0; th[2] <- gamma_inits[ch]
    out <- list(theta = unname(th),
                logs0 = log(sd(y_res - exp(as.vector(X %*% th)))),
                .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed0 + ch)
    if (!homoscedastic) out$delta <- delta_inits[ch]
    out
  })

  jm <- jags.model(textConnection(if (homoscedastic) jags_homo else jags_hetero),
                   data = dat, inits = inits, n.chains = n_chains,
                   n.adapt = n_adapt, quiet = TRUE)
  update(jm, n_burn)                               # burn-in
  samp <- coda.samples(jm, c("theta", "sigma0", "delta"),
                       n.iter = n_iter, thin = thin)
  attr(samp, "samplers") <- unique(names(list.samplers(jm)))
  samp
}

## main prior, 4 chains each
fit_M0 <- run_jags(X0, prior_mean = c(0, 4),    prior_sd = c(10, 2))
fit_M1 <- run_jags(X1, prior_mean = c(0, 4, 0), prior_sd = c(10, 2, 2))

## gamma_F prior sensitivity (M1, 2 chains each)
fit_M1_vague  <- run_jags(X1, prior_mean = c(0, 0, 0), prior_sd = c(10, 10, 2),
                          n_chains = 2)
fit_M1_inform <- run_jags(X1, prior_mean = c(0, 4, 0), prior_sd = c(10, 0.5, 2),
                          n_chains = 2)

## Confirm JAGS used the slice sampler for every parameter of this model.
cat("Samplers JAGS assigned (M1):", paste(attr(fit_M1, "samplers"), collapse = ", "), "\n")

par_labels <- c("theta[1]" = "b0", "theta[2]" = "gamma_F",
                "theta[3]" = "beta_B", "sigma0" = "sigma0", "delta" = "delta")

diag_table <- function(fit, run) {
  keep <- intersect(varnames(fit), names(par_labels))
  rhat <- gelman.diag(fit[, keep], multivariate = FALSE)$psrf[, 1]
  ess  <- effectiveSize(fit[, keep])
  data.frame(run = run, parameter = unname(par_labels[keep]),
             Rhat = round(rhat, 4), ESS = round(ess), row.names = NULL)
}

rbind(diag_table(fit_M0, "M0 (main prior)"),
      diag_table(fit_M1, "M1 (main prior)"),
      diag_table(fit_M1_vague,  "M1 (vague gamma_F prior)"),
      diag_table(fit_M1_inform, "M1 (informative gamma_F prior)"))

## traceplots + autocorrelation for the headline run (M1, main prior)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
traceplot(fit_M1[, c("theta[2]", "theta[3]", "sigma0", "delta")])
par(op)
save_current_pdf("task3-diagnostics-trace.pdf", width = 10, height = 6)

op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
acf(as.matrix(fit_M1[, "theta[2]"]), main = "ACF gamma_F (pooled chains)")
acf(as.matrix(fit_M1[, "delta"]),    main = "ACF delta (pooled chains)")
par(op)
save_current_pdf("task3-diagnostics-acf.pdf", width = 10, height = 6)

## back-transform the centred intercept to the physical k for every draw
post_summary <- function(fit, hull = TRUE) {
  d <- as.matrix(fit)
  k <- exp(d[, "theta[1]"] - d[, "theta[2]"] * mF -
             (if (hull) d[, "theta[3]"] * mB else 0))
  out <- cbind(k = k, gamma_F = d[, "theta[2]"],
               beta_B = if (hull) d[, "theta[3]"] else NULL,
               sigma0 = d[, "sigma0"], delta = d[, "delta"])
  t(apply(out, 2, function(z)
    c(mean = mean(z), sd = sd(z),
      `2.5%` = unname(quantile(z, .025)), `50%` = unname(quantile(z, .5)),
      `97.5%` = unname(quantile(z, .975)))))
}

cat("M0 (Froude only):\n")
round(post_summary(fit_M0, hull = FALSE), 4)
cat("\nM1 (Froude + beam_draught):\n")
round(post_summary(fit_M1), 4)

gF_draws <- as.matrix(fit_M1)[, "theta[2]"]
bB_draws <- as.matrix(fit_M1)[, "theta[3]"]
dl_draws <- as.matrix(fit_M1)[, "delta"]
cat(sprintf("\nP(gamma_F > 4 | y) = %.4f   P(beta_B < 0 | y) = %.4f   P(delta > 0 | y) = %.4f\n",
            mean(gF_draws > 4), mean(bB_draws < 0), mean(dl_draws > 0)))
cat("log-log (part 1) for comparison: beta_F =", round(bF, 3),
    ", 95% CI [", round(ci[1], 3), ",", round(ci[2], 3), "]\n")
cat("log-log beam_draught 95% CI: [",
    paste(round(confint(fit)["beam_draught", ], 4), collapse = ", "), "]\n")

op <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2, 1))

## gamma_F (heteroscedastic M1)
dg <- density(gF_draws)
plot(dg, col = "firebrick", lwd = 2.4, xlim = c(3.9, max(dg$x)),
     main = expression("Posterior of "*gamma[F]*" (heteroscedastic, M1)"),
     xlab = expression(gamma[F]), ylab = "posterior density")
abline(v = 4, lty = 2)
legend("topright", bty = "n", lwd = c(2.4, 1), col = c("firebrick", "black"),
       lty = c(1, 2),
       legend = c(expression(gamma[F]*" (heteroscedastic, M1)"),
                  expression(gamma[F] == 4)))

## beam_draught effect: credible once the error model is corrected
plot(density(bB_draws), col = "firebrick", lwd = 2.4,
     main = expression("Posterior of "*beta[B]*" (beam-draught)"),
     xlab = expression(beta[B]), ylab = "posterior density")
abline(v = 0, lty = 2)
abline(v = quantile(bB_draws, c(.025, .975)), col = "firebrick", lty = 3)

## variance-growth exponent delta: 1 = multiplicative (delta = 0 marker removed)
plot(density(dl_draws), col = "firebrick", lwd = 2.4,
     main = expression("Posterior of "*delta*" (variance growth)"),
     xlab = expression(delta), ylab = "posterior density")
abline(v = 1, col = "grey50", lty = 3, lwd = 1)
legend("topright", bty = "n", lwd = 1, lty = 3,
       col = "grey50",
       legend = expression(delta == 1*" (multiplicative)"))
par(op)
save_current_pdf("task3-gamma-plot.pdf", width = 13, height = 4.6)

## gamma_F posterior under the three priors (M1)
sens <- rbind(
  "N(4, 2^2)   (main)"        = post_summary(fit_M1)["gamma_F", ],
  "N(0, 10^2)  (vague)"       = post_summary(fit_M1_vague)["gamma_F", ],
  "N(4, 0.5^2) (informative)" = post_summary(fit_M1_inform)["gamma_F", ])
round(sens, 4)

if (!requireNamespace("loo", quietly = TRUE)) {
  stop("Package 'loo' is required for PSIS-LOO. Install it with install.packages('loo') before knitting.")
}
library(loo)

log_lik_from_fit <- function(fit, X, n_draws = 4000, seed = 2026) {
  d <- as.matrix(fit)
  set.seed(seed)
  d <- d[sample(nrow(d), n_draws), ]
  th  <- d[, paste0("theta[", seq_len(ncol(X)), "]"), drop = FALSE]
  eta <- th %*% t(X)                                        # n_draws x n
  mu  <- exp(eta)
  sigma <- d[, "sigma0"] * exp(d[, "delta"] * eta)          # per-obs sd, recycles by draw
  dgamma(matrix(y_res, n_draws, n_obs, byrow = TRUE),
         shape = (mu / sigma)^2, rate = mu / (sigma * sigma), log = TRUE)
}

ic_from_log_lik <- function(ll, loo_obj = loo(ll, r_eff = NA)) {
  ## WAIC: lppd via log-sum-exp, p_waic = sum of pointwise posterior variances
  mx     <- apply(ll, 2, max)
  lppd   <- sum(mx + log(colMeans(exp(sweep(ll, 2, mx)))))
  p_waic <- sum(apply(ll, 2, var))

  ## PSIS-LOO: approximate leave-one-out CV from the same log-likelihood matrix
  k <- pareto_k_values(loo_obj)

  c(WAIC = -2 * (lppd - p_waic), p_WAIC = p_waic,
    LOOIC = loo_obj$estimates["looic", "Estimate"],
    elpd_LOO = loo_obj$estimates["elpd_loo", "Estimate"],
    p_LOO = loo_obj$estimates["p_loo", "Estimate"],
    max_pareto_k = max(k),
    n_pareto_k_gt_0.7 = sum(k > 0.7))
}

ll_M0 <- log_lik_from_fit(fit_M0, X0)
ll_M1 <- log_lik_from_fit(fit_M1, X1)
loo_M0 <- loo(ll_M0, r_eff = NA)
loo_M1 <- loo(ll_M1, r_eff = NA)

ic <- rbind("M0: Froude only"           = ic_from_log_lik(ll_M0, loo_M0),
            "M1: Froude + beam_draught" = ic_from_log_lik(ll_M1, loo_M1))
round(cbind(ic,
            dWAIC = ic[, "WAIC"] - min(ic[, "WAIC"]),
            dLOOIC = ic[, "LOOIC"] - min(ic[, "LOOIC"])), 1)

loo_compare(list(
  "M0: Froude only" = loo_M0,
  "M1: Froude + beam_draught" = loo_M1
))

## same M1 mean and priors, but delta pinned to 0 (constant variance)
fit_M1_homo <- run_jags(
  X1, prior_mean = c(0, 4, 0), prior_sd = c(10, 2, 2),
  homoscedastic = TRUE
)

ll_M1_homo <- log_lik_from_fit(fit_M1_homo, X1)
loo_M1_homo <- loo(ll_M1_homo, r_eff = NA)

ic_var <- rbind(
  "M1, constant variance (delta = 0)" = ic_from_log_lik(ll_M1_homo, loo_M1_homo),
  "M1, learned variance  (delta free)" = ic_from_log_lik(ll_M1, loo_M1))
round(cbind(ic_var,
            dWAIC  = ic_var[, "WAIC"]  - min(ic_var[, "WAIC"]),
            dLOOIC = ic_var[, "LOOIC"] - min(ic_var[, "LOOIC"])), 1)

## Froude exponent under the two variance assumptions (same mean structure)
cat(sprintf("gamma_F: constant variance = %.2f   learned variance = %.2f\n",
            mean(as.matrix(fit_M1_homo)[, "theta[2]"]),
            mean(as.matrix(fit_M1)[, "theta[2]"])))

loo_compare(list(
  "M1 constant variance" = loo_M1_homo,
  "M1 learned variance"  = loo_M1
))

## posterior-mean fit and posterior-mean variance function (sigma_i = sigma0 * mu_i^delta)
d_M1      <- as.matrix(fit_M1)
th_hat    <- colMeans(d_M1[, paste0("theta[", seq_len(ncol(X1)), "]")])
sigma0_hat <- mean(d_M1[, "sigma0"])
delta_hat  <- mean(d_M1[, "delta"])

fitted_M1 <- as.vector(exp(X1 %*% th_hat))   # fitted resistance, original scale
sigma_i   <- sigma0_hat * fitted_M1^delta_hat  # per-observation fitted sd
resid_M1  <- y_res - fitted_M1               # raw residuals
std_M1    <- resid_M1 / sigma_i              # standardized by the learned variance

## quantify heteroscedasticity: residual SD within fitted-value tertiles
ter <- cut(fitted_M1, quantile(fitted_M1, c(0, 1/3, 2/3, 1)),
           include.lowest = TRUE, labels = c("low", "mid", "high"))
sd_by_tertile <- tapply(resid_M1, ter, sd)
cat("RAW residual SD by fitted-value tertile (the heteroscedasticity to absorb):\n")
print(round(sd_by_tertile, 3))
cat(sprintf("  raw high/low spread ratio = %.1f\n",
            sd_by_tertile["high"] / sd_by_tertile["low"]))

## after standardizing by the learned sigma_i, the spread should be ~flat and ~1
std_by_tertile <- tapply(std_M1, ter, sd)
cat("STANDARDIZED residual SD by fitted-value tertile (should be ~1 and flat):\n")
print(round(std_by_tertile, 3))
cat(sprintf("  cor(|standardized residual|, fitted) = %.3f   (was strongly + under constant variance)\n",
            cor(abs(std_M1), fitted_M1)))
cat(sprintf("  fraction of standardized residuals within +/-2 = %.3f\n",
            mean(abs(std_M1) < 2)))

op <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.2, 1))

## (1) residuals vs fitted
plot(fitted_M1, resid_M1, pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.45),
     xlab = "fitted resistance", ylab = "residual",
     main = "Residuals vs fitted (original scale)")
abline(h = 0, col = "firebrick"); lines(lowess(fitted_M1, resid_M1), col = "navy", lwd = 2)

## (2) scale-location: sqrt|standardized residual| vs fitted -> homoscedasticity check
plot(fitted_M1, sqrt(abs(std_M1)), pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.45),
     xlab = "fitted resistance", ylab = "sqrt(|standardized residual|)",
     main = "Scale-location")
lines(lowess(fitted_M1, sqrt(abs(std_M1))), col = "navy", lwd = 2)

## (3) standardized residuals vs the physical driver, Froude
plot(yh$Froude, std_M1, pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.45),
     xlab = "Froude", ylab = "standardized residual",
     main = "Standardized residuals vs Froude")
abline(h = c(-2, 0, 2), col = c("grey50", "firebrick", "grey50"), lty = c(3, 1, 3))

## (4) normal Q-Q of the standardized residuals
qqnorm(std_M1, pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.45),
       main = "Normal Q-Q (standardized residuals)")
qqline(std_M1, col = "firebrick")
par(op)
save_current_pdf("task3-resid-diag.pdf", width = 10, height = 7.5)

## log-scale residuals for the part-3 fit, on the same axes as the part-1 plot
logfit_M1   <- log(fitted_M1)              # = X1 %*% th_hat  (fitted log-resistance)
logresid_M1 <- log(y_res) - logfit_M1      # log(y) - log(mu_hat)

cat(sprintf("cor(|log residual|, fitted log-resistance) = %.3f\n",
            cor(abs(logresid_M1), logfit_M1)))

op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1))
plot(logfit_M1, logresid_M1, pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.45),
     xlab = "fitted log(resistance)", ylab = "log residual",
     main = "Residuals vs fitted (log scale, part 3)")
abline(h = 0, col = "firebrick")
lines(lowess(logfit_M1, logresid_M1), col = "navy", lwd = 2)

qqnorm(logresid_M1, pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.45),
       main = "Normal Q-Q (log residuals, part 3)")
qqline(logresid_M1, col = "firebrick")
par(op)
save_current_pdf("task3-resid-log.pdf", width = 10, height = 4.4)

set.seed(7)
## Posterior predictive bands for the nonlinear mean under either variance model.
## It reads the 'delta' column, so the SAME function serves the learned-variance
## fit (delta free) and the constant-variance fit (delta = 0 => sigma_i = sigma0).
predict_nonlinear <- function(fit, X_new, n_draws = 5000) {
  d <- as.matrix(fit)
  d <- d[sample(nrow(d), n_draws), ]
  theta <- d[, paste0("theta[", seq_len(ncol(X_new)), "]"), drop = FALSE]
  eta   <- theta %*% t(X_new)
  mu    <- exp(eta)
  sigma <- d[, "sigma0"] * exp(d[, "delta"] * eta)
  yrep  <- matrix(rgamma(length(mu), shape = (mu / sigma)^2, rate = mu / (sigma * sigma)),
                  nrow = n_draws)
  list(mean  = colMeans(mu),
       lower = apply(yrep, 2, quantile, 0.025),
       upper = apply(yrep, 2, quantile, 0.975))
}

## Froude grid, hull fixed at the representative beam_draught (sample median)
fr_grid <- seq(min(yh$Froude), max(yh$Froude), length.out = 160)
bd_ref  <- median(yh$beam_draught)
Xgrid   <- cbind(intercept = 1, log_Fr = log(fr_grid) - mF, beam_draught = bd_ref - mB)

pp_homo <- predict_nonlinear(fit_M1_homo, Xgrid)   # constant variance (delta = 0)
pp_het  <- predict_nonlinear(fit_M1,      Xgrid)   # learned variance  (delta free)

## 95% predictive coverage and average band width at the observed design points
cover <- function(fit) {
  pp <- predict_nonlinear(fit, X1, n_draws = 3000)
  c(coverage = mean(y_res >= pp$lower & y_res <= pp$upper),
    mean_band_width = mean(pp$upper - pp$lower))
}
round(rbind(`M1 constant variance (delta = 0)`  = cover(fit_M1_homo),
            `M1 learned variance  (delta free)` = cover(fit_M1)), 3)

ylim_top <- max(yh$resistance, pp_homo$upper, pp_het$upper)
plot(yh$Froude, yh$resistance, pch = 16, cex = 0.65, col = rgb(0, 0, 0, 0.42),
     xlab = "Froude number", ylab = "residuary resistance",
     ylim = c(0, ylim_top * 1.04),
     main = "Posterior predictive: constant vs learned variance (Task 3)")
col_homo <- rgb(0.1, 0.2, 0.7, 0.16)
col_het  <- rgb(0.8, 0.1, 0.1, 0.16)
polygon(c(fr_grid, rev(fr_grid)), c(pmax(pp_homo$lower, 0), rev(pp_homo$upper)),
        col = col_homo, border = NA)
polygon(c(fr_grid, rev(fr_grid)), c(pmax(pp_het$lower, 0), rev(pp_het$upper)),
        col = col_het, border = NA)
points(yh$Froude, yh$resistance, pch = 16, cex = 0.65, col = rgb(0, 0, 0, 0.42))
lines(fr_grid, pp_homo$mean, col = "navy",      lwd = 2.5)
lines(fr_grid, pp_het$mean,  col = "firebrick", lwd = 2.5)
legend("topleft", bty = "n", lwd = c(2.5, 2.5, NA), pch = c(NA, NA, 16),
       pt.cex = c(NA, NA, 0.8),
       col = c("navy", "firebrick", rgb(0, 0, 0, 0.42)),
       legend = c("constant variance (delta = 0): mean + 95% band",
                  "learned variance (delta free): mean + 95% band",
                  "observed data"))
save_current_pdf("task3-ppc-var.pdf", width = 9, height = 5.8)

jags_normal_hetero <- "
model {
  for (i in 1:n) {
    mu[i]    <- exp(inprod(X[i, ], theta))
    sigma[i] <- sigma0 * pow(mu[i], delta)
    y[i] ~ dnorm(mu[i], 1 / (sigma[i] * sigma[i]))
  }
  for (j in 1:p) {
    theta[j] ~ dnorm(prior_mean[j], 1 / (prior_sd[j] * prior_sd[j]))
  }
  logs0  ~ dnorm(0, 1 / 100)
  sigma0 <- exp(logs0)
  delta  ~ dnorm(0, 1)
}"

run_jags_normal <- function(X, prior_mean, prior_sd,
                            n_chains = 4, n_adapt = 1000, n_burn = 2000,
                            n_iter = 8000, thin = 2,
                            gamma_inits = c(3, 6, 5, 4),
                            delta_inits = c(0, 1, 0.5, 1.2),
                            seed0 = 3026) {
  p <- ncol(X)
  stopifnot(length(prior_mean) == p, length(prior_sd) == p)
  dat <- list(y = y_res, X = X, n = nrow(X), p = p,
              prior_mean = prior_mean, prior_sd = prior_sd)
  theta0 <- coef(lm(log(y_res) ~ X - 1))

  inits <- lapply(seq_len(n_chains), function(ch) {
    th <- theta0; th[2] <- gamma_inits[ch]
    list(theta = unname(th),
         logs0 = log(sd(y_res - exp(as.vector(X %*% th)))),
         delta = delta_inits[ch],
         .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed0 + ch)
  })

  jm <- jags.model(textConnection(jags_normal_hetero),
                   data = dat, inits = inits, n.chains = n_chains,
                   n.adapt = n_adapt, quiet = TRUE)
  update(jm, n_burn)
  samp <- coda.samples(jm, c("theta", "sigma0", "delta"),
                       n.iter = n_iter, thin = thin)
  attr(samp, "samplers") <- unique(names(list.samplers(jm)))
  samp
}

fit_M1_normal <- run_jags_normal(
  X1, prior_mean = c(0, 4, 0), prior_sd = c(10, 2, 2)
)

cat("Samplers JAGS assigned (Normal sensitivity):",
    paste(attr(fit_M1_normal, "samplers"), collapse = ", "), "\n")

lik_sens_summary <- rbind(
  `main M1: Gamma`               = post_summary(fit_M1),
  `sensitivity: additive Normal` = post_summary(fit_M1_normal)
)
round(lik_sens_summary[, c("mean", "sd", "2.5%", "50%", "97.5%")], 4)

d_normal  <- as.matrix(fit_M1_normal)
gF_normal <- d_normal[, "theta[2]"]
bB_normal <- d_normal[, "theta[3]"]
dl_normal <- d_normal[, "delta"]

cat(sprintf(
  "additive Normal sensitivity: P(gamma_F > 4 | y) = %.4f   P(beta_B < 0 | y) = %.4f   P(delta > 0 | y) = %.4f\n",
  mean(gF_normal > 4), mean(bB_normal < 0), mean(dl_normal > 0)
))

op <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2, 1))
plot(density(gF_draws), col = "firebrick", lwd = 2.4,
     xlim = range(gF_draws, gF_normal),
     main = expression("Likelihood sensitivity: "*gamma[F]),
     xlab = expression(gamma[F]), ylab = "posterior density")
lines(density(gF_normal), col = "navy", lwd = 2.4)
abline(v = 4, lty = 2)
legend("topright", bty = "n", lwd = c(2.4, 2.4, 1),
       col = c("firebrick", "navy", "black"), lty = c(1, 1, 2),
       legend = c("Gamma M1 (main)", "additive Normal", expression(gamma[F] == 4)))

plot(density(bB_draws), col = "firebrick", lwd = 2.4,
     xlim = range(bB_draws, bB_normal),
     main = expression("Likelihood sensitivity: "*beta[B]),
     xlab = expression(beta[B]), ylab = "posterior density")
lines(density(bB_normal), col = "navy", lwd = 2.4)
abline(v = 0, lty = 2)

plot(density(dl_draws), col = "firebrick", lwd = 2.4,
     xlim = range(dl_draws, dl_normal),
     main = expression("Likelihood sensitivity: "*delta),
     xlab = expression(delta), ylab = "posterior density")
lines(density(dl_normal), col = "navy", lwd = 2.4)
abline(v = c(0, 1), lty = c(2, 3), col = c("black", "grey50"))
par(op)
save_current_pdf("task3-gamma-likelihood-sensitivity.pdf", width = 13, height = 4.6)

## additive-Normal pointwise log-likelihood, same draws structure as log_lik_from_fit().
log_lik_normal_from_fit <- function(fit, X, n_draws = 4000, seed = 2026) {
  d <- as.matrix(fit)
  set.seed(seed)
  d   <- d[sample(nrow(d), n_draws), ]
  th  <- d[, paste0("theta[", seq_len(ncol(X)), "]"), drop = FALSE]
  eta <- th %*% t(X)                                        # n_draws x n
  mu  <- exp(eta)
  sigma <- d[, "sigma0"] * exp(d[, "delta"] * eta)          # per-obs sd
  dnorm(matrix(y_res, n_draws, n_obs, byrow = TRUE),
        mean = mu, sd = sigma, log = TRUE)
}

ll_M1_normal  <- log_lik_normal_from_fit(fit_M1_normal, X1)
loo_M1_normal <- loo(ll_M1_normal, r_eff = NA)

ic_lik <- rbind(
  "M1: Gamma (main)"    = ic_from_log_lik(ll_M1,        loo_M1),
  "M1: additive Normal" = ic_from_log_lik(ll_M1_normal, loo_M1_normal))
round(cbind(ic_lik,
            dWAIC  = ic_lik[, "WAIC"]  - min(ic_lik[, "WAIC"]),
            dLOOIC = ic_lik[, "LOOIC"] - min(ic_lik[, "LOOIC"])), 1)

## elpd difference with Monte-Carlo standard error (decides if the gap is real)
loo_compare(list(
  "M1 Gamma"           = loo_M1,
  "M1 additive Normal" = loo_M1_normal
))

## additive-Normal posterior-predictive analogue of predict_nonlinear().
predict_normal <- function(fit, X_new, n_draws = 3000) {
  d <- as.matrix(fit)
  d <- d[sample(nrow(d), n_draws), ]
  theta <- d[, paste0("theta[", seq_len(ncol(X_new)), "]"), drop = FALSE]
  eta   <- theta %*% t(X_new)
  mu    <- exp(eta)
  sigma <- d[, "sigma0"] * exp(d[, "delta"] * eta)
  yrep  <- mu + matrix(rnorm(length(mu), sd = sigma), nrow = n_draws)
  list(mean  = colMeans(mu),
       lower = apply(yrep, 2, quantile, 0.025),
       upper = apply(yrep, 2, quantile, 0.975),
       p_neg = mean(yrep < 0))
}

pp_n <- predict_normal(fit_M1_normal, X1)
cmp <- rbind(
  "Gamma (main)"    = c(cover(fit_M1), P_yrep_below_0 = 0),
  "additive Normal" = c(coverage        = mean(y_res >= pp_n$lower & y_res <= pp_n$upper),
                        mean_band_width = mean(pp_n$upper - pp_n$lower),
                        P_yrep_below_0  = pp_n$p_neg))
round(cmp, 4)

set.seed(42)
pp_gamma_grid <- predict_nonlinear(fit_M1,        Xgrid)   # Gamma M1 (main)
pp_norm_grid  <- predict_normal(fit_M1_normal,    Xgrid)   # additive Normal sensitivity

ylim_top <- max(yh$resistance, pp_norm_grid$upper, pp_gamma_grid$upper)
plot(yh$Froude, yh$resistance, pch = 16, cex = 0.65, col = rgb(0, 0, 0, 0.42),
     xlab = "Froude number", ylab = "residuary resistance",
     ylim = c(min(0, pp_norm_grid$lower), ylim_top * 1.04),
     main = "Posterior predictive: Gamma (main) vs additive Normal (Task 3)")
col_gamma <- rgb(0.0, 0.5, 0.0, 0.15)
col_norm  <- rgb(0.8, 0.1, 0.1, 0.15)
polygon(c(fr_grid, rev(fr_grid)),
        c(pp_gamma_grid$lower, rev(pp_gamma_grid$upper)),
        col = col_gamma, border = NA)
polygon(c(fr_grid, rev(fr_grid)),
        c(pp_norm_grid$lower, rev(pp_norm_grid$upper)),
        col = col_norm, border = NA)
abline(h = 0, col = "grey60", lty = 3)
## un-clipped Normal lower band, to make the negative excursion explicit
lines(fr_grid, pp_norm_grid$lower, col = "firebrick", lwd = 1, lty = 2)
points(yh$Froude, yh$resistance, pch = 16, cex = 0.65, col = rgb(0, 0, 0, 0.42))
lines(fr_grid, pp_gamma_grid$mean, col = "darkgreen",  lwd = 2.5)
lines(fr_grid, pp_norm_grid$mean,  col = "firebrick",  lwd = 2.5)
legend("topleft", bty = "n",
       lwd = c(2.5, 2.5, 1, NA), lty = c(1, 1, 2, NA), pch = c(NA, NA, NA, 16),
       pt.cex = c(NA, NA, NA, 0.8),
       col = c("darkgreen", "firebrick", "firebrick", rgb(0, 0, 0, 0.42)),
       legend = c("Gamma (main): mean + 95% band",
                  "additive Normal: mean + 95% band",
                  "Normal lower band (un-clipped, dips < 0)",
                  "observed data"))
save_current_pdf("task3-gamma-predictive-overlay.pdf", width = 9, height = 5.8)

set.seed(42)

## Froude grid and representative hull values
fr_grid <- seq(min(yh$Froude), max(yh$Froude), length.out = 160)
bd_ref <- median(yh$beam_draught)
hull_ref <- as.list(colMeans(yh[, c("LCB", "prismatic", "L_disp", "beam_draught", "L_beam")]))
hull_ref$beam_draught <- bd_ref

new_loglog <- data.frame(log_Fr = log(fr_grid), hull_ref)
new_M1 <- cbind(intercept = 1,
                log_Fr = log(fr_grid) - mF,
                beam_draught = bd_ref - mB)

## Draw posterior predictive samples for the reference-prior log-log model.
## beta | sigma2, y is normal; sigma2 | y is inverse-gamma.
draw_loglog_pred <- function(fit, newdata, n_draws = 5000) {
  X <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  p <- ncol(X)
  df_res <- fit$df.residual
  beta_hat <- coef(fit)
  s2 <- sum(resid(fit)^2) / df_res
  XtX_inv <- chol2inv(chol(crossprod(X)))
  X_new <- model.matrix(delete.response(terms(fit)), newdata)

  sigma2 <- df_res * s2 / rchisq(n_draws, df = df_res)
  z <- matrix(rnorm(n_draws * p), nrow = n_draws)
  beta <- sweep(z %*% chol(XtX_inv), 1, sqrt(sigma2), `*`)
  beta <- sweep(beta, 2, beta_hat, `+`)

  eta <- beta %*% t(X_new)
  yrep <- exp(eta + matrix(rnorm(length(eta), sd = sqrt(sigma2)),
                          nrow = n_draws, ncol = nrow(X_new)))
  mu <- exp(eta + sigma2 / 2)   # posterior mean of resistance, not median

  list(mean = colMeans(mu),
       lower = apply(yrep, 2, quantile, 0.025),
       upper = apply(yrep, 2, quantile, 0.975))
}

## Draw posterior predictive samples for the heteroscedastic nonlinear M1 model.
## The main Task 3 likelihood is Gamma, parameterized by mean mu and
## sd sigma_i = sigma0 * mu_i^delta.
draw_M1_pred <- function(fit, X_new, n_draws = 5000) {
  d <- as.matrix(fit)
  d <- d[sample(seq_len(nrow(d)), n_draws), ]
  theta <- d[, paste0("theta[", seq_len(ncol(X_new)), "]"), drop = FALSE]

  eta <- theta %*% t(X_new)
  mu  <- exp(eta)
  sigma <- d[, "sigma0"] * exp(d[, "delta"] * eta)   # per-draw, per-point sd
  yrep <- matrix(rgamma(length(mu), shape = (mu / sigma)^2,
                        rate = mu / (sigma * sigma)),
                 nrow = n_draws, ncol = nrow(X_new))

  list(mean = colMeans(mu),
       lower = apply(yrep, 2, quantile, 0.025),
       upper = apply(yrep, 2, quantile, 0.975))
}

pp_loglog <- draw_loglog_pred(fit, new_loglog)
pp_M1_homo <- draw_M1_pred(fit_M1_homo, new_M1)
pp_M1 <- draw_M1_pred(fit_M1, new_M1)

## Coverage is evaluated at the observed design points, using the same predictive rule.
obs_loglog <- data.frame(log_Fr = yh$log_Fr,
                         LCB = yh$LCB, prismatic = yh$prismatic,
                         L_disp = yh$L_disp, beam_draught = yh$beam_draught,
                         L_beam = yh$L_beam)
obs_M1 <- X1
cover_loglog <- draw_loglog_pred(fit, obs_loglog, n_draws = 3000)
cover_M1_homo <- draw_M1_pred(fit_M1_homo, obs_M1, n_draws = 3000)
cover_M1 <- draw_M1_pred(fit_M1, obs_M1, n_draws = 3000)
coverage <- rbind(
  `log-log, multiplicative errors` = c(
    coverage = mean(yh$resistance >= cover_loglog$lower & yh$resistance <= cover_loglog$upper),
    mean_band_width = mean(cover_loglog$upper - cover_loglog$lower)),
  `Gamma M1, constant variance` = c(
    coverage = mean(yh$resistance >= cover_M1_homo$lower & yh$resistance <= cover_M1_homo$upper),
    mean_band_width = mean(cover_M1_homo$upper - cover_M1_homo$lower)),
  `nonlinear M1, Gamma learned variance` = c(
    coverage = mean(yh$resistance >= cover_M1$lower & yh$resistance <= cover_M1$upper),
    mean_band_width = mean(cover_M1$upper - cover_M1$lower))
)
round(coverage, 3)

pp_widths <- data.frame(
  Froude = fr_grid,
  loglog = pp_loglog$upper - pp_loglog$lower,
  gamma_constant = pp_M1_homo$upper - pp_M1_homo$lower,
  gamma_learned = pp_M1$upper - pp_M1$lower
)

plot_task4_widths <- function(file = NULL) {
  if (!is.null(file)) pdf(file, width = 8.5, height = 5.2)
  yr <- range(pp_widths[, -1], finite = TRUE)
  plot(pp_widths$Froude, pp_widths$loglog, type = "l", lwd = 2.4,
       col = "navy", ylim = yr, xlab = "Froude number",
       ylab = "95% posterior predictive band width",
       main = "Predictive uncertainty across the Froude range")
  lines(pp_widths$Froude, pp_widths$gamma_constant, lwd = 2.4,
        col = "darkorange3")
  lines(pp_widths$Froude, pp_widths$gamma_learned, lwd = 2.4,
        col = "firebrick")
  legend("topleft", bty = "n", lwd = 2.4,
         col = c("navy", "darkorange3", "firebrick"),
         legend = c("log-log", "Gamma M1 constant variance",
                    "Gamma M1 learned variance"))
  if (!is.null(file)) dev.off()
}

plot_task4_widths()
plot_task4_widths("figures/posterior_predictive_widths.pdf")

## Plot on the original scale. Both predictive families have positive support:
## log-log is lognormal after back-transformation and M1 is Gamma.
ylim_top <- max(yh$resistance, pp_loglog$upper, pp_M1$upper)
plot(yh$Froude, yh$resistance, pch = 16, cex = 0.65,
     col = rgb(0, 0, 0, 0.42), xlab = "Froude number",
     ylab = "residuary resistance", ylim = c(0, ylim_top * 1.04),
     main = "Posterior predictive curves on the original scale")

col_log <- rgb(0.1, 0.2, 0.7, 0.18)
col_add <- rgb(0.8, 0.1, 0.1, 0.18)
polygon(c(fr_grid, rev(fr_grid)),
        c(pmax(pp_loglog$lower, 0), rev(pp_loglog$upper)),
        col = col_log, border = NA)
polygon(c(fr_grid, rev(fr_grid)),
        c(pmax(pp_M1$lower, 0), rev(pp_M1$upper)),
        col = col_add, border = NA)
points(yh$Froude, yh$resistance, pch = 16, cex = 0.65,
       col = rgb(0, 0, 0, 0.42))
lines(fr_grid, pp_loglog$mean, col = "navy", lwd = 2.5)
lines(fr_grid, pp_M1$mean, col = "firebrick", lwd = 2.5)
legend("topleft", bty = "n", lwd = c(2.5, 2.5, NA), pch = c(NA, NA, 16),
       pt.cex = c(NA, NA, 0.8),
       col = c("navy", "firebrick", rgb(0, 0, 0, 0.42)),
       legend = c("log-log predictive mean + 95% band",
                  "heteroscedastic nonlinear M1 predictive mean + 95% band",
                  "observed data"))

pdf("figures/posterior_predictive_curves.pdf", width = 9, height = 5.8)
plot(yh$Froude, yh$resistance, pch = 16, cex = 0.65,
     col = rgb(0, 0, 0, 0.42), xlab = "Froude number",
     ylab = "residuary resistance", ylim = c(0, ylim_top * 1.04),
     main = "Posterior predictive curves on the original scale")
polygon(c(fr_grid, rev(fr_grid)),
        c(pmax(pp_loglog$lower, 0), rev(pp_loglog$upper)),
        col = col_log, border = NA)
polygon(c(fr_grid, rev(fr_grid)),
        c(pmax(pp_M1$lower, 0), rev(pp_M1$upper)),
        col = col_add, border = NA)
points(yh$Froude, yh$resistance, pch = 16, cex = 0.65,
       col = rgb(0, 0, 0, 0.42))
lines(fr_grid, pp_loglog$mean, col = "navy", lwd = 2.5)
lines(fr_grid, pp_M1$mean, col = "firebrick", lwd = 2.5)
legend("topleft", bty = "n", lwd = c(2.5, 2.5, NA), pch = c(NA, NA, 16),
       pt.cex = c(NA, NA, 0.8),
       col = c("navy", "firebrick", rgb(0, 0, 0, 0.42)),
       legend = c("log-log predictive mean + 95% band",
                  "heteroscedastic nonlinear M1 predictive mean + 95% band",
                  "observed data"))
dev.off()

yh_train <- read.csv("data/yacht_train.csv")
yh_test  <- read.csv("data/yacht_test.csv")

yh_train$log_res <- log(yh_train$resistance)
yh_train$log_Fr  <- log(yh_train$Froude)
yh_test$log_res  <- log(yh_test$resistance)
yh_test$log_Fr   <- log(yh_test$Froude)

log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

## Test-set posterior predictive log density for the log-log model.
## The density is evaluated on the original resistance scale, so the predictive
## density is lognormal after back-transforming the Gaussian log model.
loglog_test_score <- function(fit, newdata, y_new, n_draws = 5000, seed = 2026) {
  X <- model.matrix(fit)
  p <- ncol(X)
  df_res <- fit$df.residual
  beta_hat <- coef(fit)
  s2 <- sum(resid(fit)^2) / df_res
  XtX_inv <- chol2inv(chol(crossprod(X)))
  X_new <- model.matrix(delete.response(terms(fit)), newdata)

  set.seed(seed)
  sigma2 <- df_res * s2 / rchisq(n_draws, df = df_res)
  z <- matrix(rnorm(n_draws * p), nrow = n_draws)
  beta <- sweep(z %*% chol(XtX_inv), 1, sqrt(sigma2), `*`)
  beta <- sweep(beta, 2, beta_hat, `+`)
  eta <- beta %*% t(X_new)

  ll <- dlnorm(matrix(y_new, n_draws, length(y_new), byrow = TRUE),
               meanlog = eta, sdlog = sqrt(sigma2), log = TRUE)
  pp <- draw_loglog_pred(fit, newdata, n_draws = 3000)

  c(total_lpd = sum(apply(ll, 2, log_mean_exp)),
    mean_lpd = mean(apply(ll, 2, log_mean_exp)),
    coverage_95 = mean(y_new >= pp$lower & y_new <= pp$upper),
    mean_band_width = mean(pp$upper - pp$lower))
}

## Test-set posterior predictive log density for the nonlinear Gamma models.
## For learned variance, sigma_i = sigma0 * mu_i^delta; for the homoscedastic
## comparison delta is pinned to zero in fit_M1_homo.
nonlinear_test_score <- function(fit, X_new, y_new, n_draws = 5000, seed = 2026) {
  d <- as.matrix(fit)
  set.seed(seed)
  d <- d[sample(nrow(d), n_draws), ]
  theta <- d[, paste0("theta[", seq_len(ncol(X_new)), "]"), drop = FALSE]
  eta <- theta %*% t(X_new)
  mu <- exp(eta)
  sigma <- d[, "sigma0"] * exp(d[, "delta"] * eta)

  ll <- dgamma(matrix(y_new, n_draws, length(y_new), byrow = TRUE),
               shape = (mu / sigma)^2, rate = mu / (sigma * sigma), log = TRUE)
  pp <- predict_nonlinear(fit, X_new, n_draws = 3000)

  c(total_lpd = sum(apply(ll, 2, log_mean_exp)),
    mean_lpd = mean(apply(ll, 2, log_mean_exp)),
    coverage_95 = mean(y_new >= pp$lower & y_new <= pp$upper),
    mean_band_width = mean(pp$upper - pp$lower))
}

test_loglog <- data.frame(log_Fr = yh_test$log_Fr,
                          LCB = yh_test$LCB,
                          prismatic = yh_test$prismatic,
                          L_disp = yh_test$L_disp,
                          beam_draught = yh_test$beam_draught,
                          L_beam = yh_test$L_beam)
test_X1 <- cbind(intercept = 1,
                 log_Fr = log(yh_test$Froude) - mF,
                 beam_draught = yh_test$beam_draught - mB)

validation_table <- rbind(
  `log-log posterior predictive` = loglog_test_score(fit, test_loglog, yh_test$resistance),
  `Gamma M1 constant variance posterior predictive` = nonlinear_test_score(fit_M1_homo, test_X1, yh_test$resistance),
  `Gamma M1 learned variance posterior predictive` = nonlinear_test_score(fit_M1, test_X1, yh_test$resistance)
)

round(validation_table, 3)

## Visual validation on the held-out rows: observed test responses against the
## posterior predictive intervals from each model. Predictive means near the
## black points and intervals containing them indicate calibrated test prediction.
test_pred_loglog <- draw_loglog_pred(fit, test_loglog, n_draws = 3000)
test_pred_homo <- predict_nonlinear(fit_M1_homo, test_X1, n_draws = 3000)
test_pred_het <- predict_nonlinear(fit_M1, test_X1, n_draws = 3000)

plot_validation <- function(file = NULL) {
  if (!is.null(file)) pdf(file, width = 8.5, height = 5.5)

  pred <- list(test_pred_loglog, test_pred_homo, test_pred_het)
  means <- vapply(pred, `[[`, numeric(length(yh_test$resistance)), "mean")
  lowers <- vapply(pred, `[[`, numeric(length(yh_test$resistance)), "lower")
  uppers <- vapply(pred, `[[`, numeric(length(yh_test$resistance)), "upper")

  ord <- order(yh_test$resistance)
  x <- seq_along(ord)
  y_obs <- yh_test$resistance[ord]
  cols <- c("navy", "darkorange3", "firebrick")

  yr <- range(y_obs, lowers, uppers, finite = TRUE)
  plot(x, y_obs, pch = 16, col = rgb(0, 0, 0, 0.55), ylim = yr,
       xlab = "test observation (ordered by observed resistance)",
       ylab = "residuary resistance",
       main = "Train/test posterior predictive validation")
  grid(col = "grey88")

  offset <- c(-0.18, 0, 0.18)
  for (j in seq_along(pred)) {
    segments(x + offset[j], lowers[ord, j], x + offset[j], uppers[ord, j],
             col = adjustcolor(cols[j], alpha.f = 0.35), lwd = 1.4)
    points(x + offset[j], means[ord, j], pch = 16, cex = 0.55, col = cols[j])
  }
  points(x, y_obs, pch = 16, col = rgb(0, 0, 0, 0.55))
  legend("topleft", bty = "n", pch = c(16, 16, 16, 16),
         col = c(rgb(0, 0, 0, 0.55), cols),
         legend = c("observed test response", "log-log predictive mean + 95% interval",
                    "Gamma M1 constant variance predictive mean + 95% interval",
                    "Gamma M1 learned variance predictive mean + 95% interval"))

  if (!is.null(file)) dev.off()
}

plot_validation()
plot_validation("figures/train_test_validation.pdf")

plot_validation_widths <- function(file = NULL) {
  if (!is.null(file)) pdf(file, width = 8.5, height = 5.2)
  pred <- list(test_pred_loglog, test_pred_homo, test_pred_het)
  widths <- vapply(pred, function(p) p$upper - p$lower,
                   numeric(length(yh_test$resistance)))
  ord <- order(yh_test$resistance)
  x <- seq_along(ord)
  yr <- range(widths, finite = TRUE)
  plot(x, widths[ord, 1], type = "l", lwd = 2.2, col = "navy",
       ylim = yr, xlab = "test observation (ordered by observed resistance)",
       ylab = "95% posterior predictive interval width",
       main = "Held-out predictive interval widths")
  lines(x, widths[ord, 2], lwd = 2.2, col = "darkorange3")
  lines(x, widths[ord, 3], lwd = 2.2, col = "firebrick")
  legend("topleft", bty = "n", lwd = 2.2,
         col = c("navy", "darkorange3", "firebrick"),
         legend = c("log-log", "Gamma M1 constant variance",
                    "Gamma M1 learned variance"))
  if (!is.null(file)) dev.off()
}

plot_validation_widths()
plot_validation_widths("figures/train_test_validation_widths.pdf")
