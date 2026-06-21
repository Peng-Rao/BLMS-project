# Section 1: log-log model and prior sensitivity

# Yacht hydrodynamics — UCI #243
yh = read.csv("data/yacht_hydro.csv")
dim(yh)
head(yh)
summary(yh)

# The power-law in Fn is best seen on log-log axes;
# resistance > 0 everywhere (min = 0.01), so plain log() needs no offset
draw_power_law_scaling <- function() {
  plot(log(yh$Froude), log(yh$resistance),
       pch = 16, cex = 0.5,
       xlab = "log(Froude)", ylab = "log(resistance)",
       main = "Yacht hydrodynamic resistance — power-law scaling")
  abline(lm(log(resistance) ~ log(Froude), data = yh),
         col = "red", lwd = 2)
}

if (interactive() || dev.cur() > 1) {
  draw_power_law_scaling()
}

pdf("figures/yacht_power_law_scaling.pdf", width = 7.6, height = 4.6)
draw_power_law_scaling()
invisible(dev.off())

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
