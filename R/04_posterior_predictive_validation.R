# Sections 4-5: posterior predictive curves and train/test validation

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
