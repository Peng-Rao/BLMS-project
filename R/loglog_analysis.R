# =====================================================================
#  Part 1 — Bayesian log-log regression for yacht residuary resistance
#
#  Model:   log(resistance) = b0 + bF*log(Froude) + sum_k bk*hull_k + eps
#  Prior:   reference (Jeffreys) prior  p(beta, sigma^2) ~ 1/sigma^2
#  Posterior is Normal-Inverse-Gamma; we draw from it by Monte Carlo.
#
#  Goal: verify the physical expectation bF ~ 4.
# =====================================================================

set.seed(2026)

## ---- Data ----------------------------------------------------------
yh <- read.csv("data/yacht_hydro.csv")
yh$log_res <- log(yh$resistance)     # min(resistance) = 0.01 > 0, no offset needed
yh$log_Fr  <- log(yh$Froude)

predictors <- c("log_Fr", "LCB", "prismatic", "L_disp", "beam_draught", "L_beam")
form <- as.formula(paste("log_res ~", paste(predictors, collapse = " + ")))

## ---- Design matrix and reference-prior posterior -------------------
X  <- model.matrix(form, data = yh)
y  <- yh$log_res
n  <- nrow(X)
p  <- ncol(X)

XtX_inv <- solve(crossprod(X))
beta_hat <- as.vector(XtX_inv %*% crossprod(X, y))   # OLS = posterior mode/mean
names(beta_hat) <- colnames(X)
resid <- y - X %*% beta_hat
s2 <- as.numeric(crossprod(resid)) / (n - p)          # unbiased error variance

## ---- Monte Carlo draws from the posterior --------------------------
##  sigma^2 | y      ~ Inv-Gamma( (n-p)/2 , (n-p)s2/2 )
##  beta | sigma^2,y ~ N( beta_hat , sigma^2 (X'X)^{-1} )
M <- 20000
sigma2_draws <- (n - p) * s2 / rchisq(M, df = n - p)   # scaled inverse-chi^2

L <- chol(XtX_inv)                                     # XtX_inv = L'L
beta_draws <- matrix(0, M, p, dimnames = list(NULL, colnames(X)))
for (m in seq_len(M)) {
  z <- rnorm(p)
  beta_draws[m, ] <- beta_hat + sqrt(sigma2_draws[m]) * as.vector(crossprod(L, z))
}

## ---- Posterior summaries -------------------------------------------
post_mean <- colMeans(beta_draws)
post_sd   <- apply(beta_draws, 2, sd)
post_ci   <- t(apply(beta_draws, 2, quantile, probs = c(0.025, 0.975)))
prob_pos  <- colMeans(beta_draws > 0)

summ <- data.frame(
  mean   = round(post_mean, 4),
  sd     = round(post_sd, 4),
  q2.5   = round(post_ci[, 1], 4),
  q97.5  = round(post_ci[, 2], 4),
  P_gt0  = round(prob_pos, 3)
)
cat("=== Posterior summary (reference prior, ", M, " MC draws) ===\n", sep = "")
print(summ)
cat("\nR^2 =", round(1 - sum(resid^2) / sum((y - mean(y))^2), 4),
    "  residual sd s =", round(sqrt(s2), 4), "\n")

bF <- beta_draws[, "log_Fr"]
cat(sprintf("\nbeta_F : mean=%.3f  95%% CI=[%.3f, %.3f]  P(beta_F>4)=%.3f\n",
            mean(bF), quantile(bF, .025), quantile(bF, .975), mean(bF > 4)))

## ====================================================================
##  Figures
## ====================================================================
dir.create("figures", showWarnings = FALSE)

## (1) log-log scatter with posterior-mean line and 95% credible band
xg <- seq(min(yh$log_Fr), max(yh$log_Fr), length.out = 120)
hull_means <- colMeans(yh[, c("LCB","prismatic","L_disp","beam_draught","L_beam")])
Xg <- cbind(1, xg, matrix(rep(hull_means, each = length(xg)), nrow = length(xg)))
colnames(Xg) <- colnames(X)
lin_draws <- Xg %*% t(beta_draws)                 # length(xg) x M
band <- apply(lin_draws, 1, quantile, probs = c(0.025, 0.5, 0.975))

pdf("figures/loglog_scatter.pdf", width = 6.5, height = 4.8)
par(mar = c(4.2, 4.2, 1, 1))
plot(yh$log_Fr, yh$log_res, pch = 16, cex = 0.55, col = rgb(0,0,0,0.45),
     xlab = "log(Froude)", ylab = "log(resistance)")
polygon(c(xg, rev(xg)), c(band[1, ], rev(band[3, ])),
        col = rgb(0.85, 0.1, 0.1, 0.18), border = NA)
lines(xg, band[2, ], col = "firebrick", lwd = 2.2)
legend("topleft", bty = "n",
       legend = c("data", "posterior mean fit", "95% credible band"),
       pch = c(16, NA, NA), lwd = c(NA, 2.2, 8),
       col = c(rgb(0,0,0,0.45), "firebrick", rgb(0.85,0.1,0.1,0.18)))
dev.off()

## (2) Monte Carlo posterior of beta_F
pdf("figures/posterior_betaF.pdf", width = 6.0, height = 4.4)
par(mar = c(4.2, 4.2, 1, 1))
h <- hist(bF, breaks = 50, plot = FALSE)
plot(h, freq = FALSE, col = "grey85", border = "grey60",
     main = "", xlab = expression(beta[F]~"(coefficient on log Froude)"), ylab = "posterior density")
lines(density(bF), col = "firebrick", lwd = 2)
ci <- quantile(bF, c(.025, .975))
abline(v = ci, col = "firebrick", lty = 2)
abline(v = 4, col = "navy", lwd = 2)
abline(v = mean(bF), col = "firebrick", lwd = 1.5)
legend("topright", bty = "n",
       legend = c("posterior density", "95% credible interval",
                  expression("physical expectation "*beta[F]==4)),
       lwd = c(2, 1, 2), lty = c(1, 2, 1),
       col = c("firebrick", "firebrick", "navy"))
dev.off()

## (3) Residual diagnostics (motivate heteroscedasticity)
fitted_vals <- as.vector(X %*% beta_hat)
pdf("figures/residual_diag.pdf", width = 7.5, height = 3.6)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
plot(fitted_vals, resid, pch = 16, cex = 0.5, col = rgb(0,0,0,0.45),
     xlab = "fitted log(resistance)", ylab = "residual", main = "Residuals vs fitted")
abline(h = 0, col = "firebrick", lwd = 1.5)
lines(lowess(fitted_vals, resid), col = "navy", lwd = 2)
qqnorm(resid, pch = 16, cex = 0.5, col = rgb(0,0,0,0.45), main = "Normal Q-Q")
qqline(resid, col = "firebrick", lwd = 1.5)
dev.off()

## (4) Posterior 95% intervals for the hull covariates (forest plot)
hull <- c("LCB","prismatic","L_disp","beam_draught","L_beam")
fp_mean <- post_mean[hull]; fp_ci <- post_ci[hull, ]
pdf("figures/coef_forest.pdf", width = 6.2, height = 3.8)
par(mar = c(4.2, 6.5, 1, 1))
yy <- seq_along(hull)
plot(fp_mean, yy, xlim = range(fp_ci), ylim = c(0.5, length(hull) + 0.5),
     pch = 16, yaxt = "n", ylab = "", xlab = "posterior coefficient (log-log model)")
abline(v = 0, lty = 2, col = "grey50")
segments(fp_ci[, 1], yy, fp_ci[, 2], yy, lwd = 2, col = "firebrick")
points(fp_mean, yy, pch = 16)
axis(2, at = yy, labels = hull, las = 1)
dev.off()

cat("\nFigures written to figures/: loglog_scatter.pdf, posterior_betaF.pdf,",
    "residual_diag.pdf, coef_forest.pdf\n")
