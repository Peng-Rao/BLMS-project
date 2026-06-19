# Section 3: Gamma nonlinear model and sensitivity checks

## centre covariates at a y^2-weighted mean: removes the strong posterior
## correlation between b0 and gamma_F that ruins mixing with plain centring
w_info <- y_res^2 / sum(y_res^2)
mF <- sum(w_info * log(yh$Froude))
mB <- sum(w_info * yh$beam_draught)
xF_c <- log(yh$Froude)  - mF
xB_c <- yh$beam_draught - mB

X0 <- cbind(intercept = 1, log_Fr = xF_c)                    # M0: Froude only
X1 <- cbind(X0, beam_draught = xB_c)                         # M1: + beam_draught

dir.create("figures", showWarnings = FALSE)

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
pdf("figures/task3-diagnostics-trace.pdf", width = 10, height = 6)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
traceplot(fit_M1[, c("theta[2]", "theta[3]", "sigma0", "delta")])
par(op)
dev.off()

pdf("figures/task3-diagnostics-acf.pdf", width = 10, height = 6)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
acf(as.matrix(fit_M1[, "theta[2]"]), main = "ACF gamma_F (pooled chains)")
acf(as.matrix(fit_M1[, "delta"]),    main = "ACF delta (pooled chains)")
par(op)
dev.off()

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

pdf("figures/task3-gamma-plot.pdf", width = 13, height = 4.6)
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
dev.off()

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

pdf("figures/task3-resid-diag.pdf", width = 10, height = 7.5)
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
dev.off()

## log-scale residuals for the part-3 fit, on the same axes as the part-1 plot
logfit_M1   <- log(fitted_M1)              # = X1 %*% th_hat  (fitted log-resistance)
logresid_M1 <- log(y_res) - logfit_M1      # log(y) - log(mu_hat)

cat(sprintf("cor(|log residual|, fitted log-resistance) = %.3f\n",
            cor(abs(logresid_M1), logfit_M1)))

pdf("figures/task3-resid-log.pdf", width = 10, height = 4.4)
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
dev.off()

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
pdf("figures/task3-ppc-var.pdf", width = 9, height = 5.8)
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
dev.off()

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

pdf("figures/task3-gamma-likelihood-sensitivity.pdf", width = 13, height = 4.6)
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
dev.off()

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
pdf("figures/task3-gamma-predictive-overlay.pdf", width = 9, height = 5.8)
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
dev.off()

set.seed(42)

## Froude grid and representative hull values
fr_grid <- seq(min(yh$Froude), max(yh$Froude), length.out = 160)
bd_ref <- median(yh$beam_draught)
