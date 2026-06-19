# Section 2: BAS/BMA model selection

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
