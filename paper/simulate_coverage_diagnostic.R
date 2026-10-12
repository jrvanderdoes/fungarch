# Diagnose lower-tail coverage using data generated from known parameters.
#
# The oracle result tests the simulator and prediction scaling. The fitted
# result tests whether .fit_diag_model_new() can recover the forecast scale.

if (!requireNamespace("fungarch", quietly = TRUE)) {
  stop("Install fungarch before running this diagnostic.", call. = FALSE)
}

simulate_coverage_diagnostic <- function(
    N = 1000,
    resolution = 20,
    burnin = 100,
    training = 500,
    alpha_levels = c(0.05, 0.01),
    alpha_value = 0.20,
    w = 0.10,
    seed = 12345) {
  if (!is.null(seed)) set.seed(seed)
  if (training < 3L || training >= N) {
    stop("training must be between 3 and N - 1.", call. = FALSE)
  }

  C_eps <- fungarch:::.generate_ceps(resolution)
  eig_Ce <- eigen(C_eps)
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  hkk <- as.vector(basis[, 1] %o% basis[, 1])
  alpha_operator <- alpha_value * (hkk %o% hkk)
  true_alpha <- matrix(alpha_value, nrow = 1, ncol = 1)

  simulated <- fungarch:::generateFARCH_use(
    N = N,
    resolution = resolution,
    alpha = alpha_operator,
    delta = C_eps,
    basis = basis,
    burnin = burnin,
    p = 1,
    C_eps = C_eps,
    kl_errors = FALSE
  )

  # Use every tenth forecast origin so the diagnostic remains quick.
  forecast_indices <- seq.int(training + 1L, N, by = 10L)
  oracle_violations <- matrix(
    NA_real_, nrow = length(forecast_indices),
    ncol = length(alpha_levels)
  )
  fitted_violations <- matrix(
    NA_real_, nrow = length(forecast_indices),
    ncol = length(alpha_levels)
  )
  colnames(oracle_violations) <- paste0("alpha_", alpha_levels)
  colnames(fitted_violations) <- paste0("alpha_", alpha_levels)
  prediction_difference <- numeric(length(forecast_indices))

  for (i in seq_along(forecast_indices)) {
    t <- forecast_indices[i]
    X_train <- simulated$data[, seq_len(t - 1L), drop = FALSE]
    X_next <- simulated$data[, t]

    oracle_sigma <- fungarch:::.pred_diag_new(
      resolution = resolution,
      alphas = true_alpha,
      basis = basis,
      delta = C_eps,
      data = X_train,
      C_eps = C_eps
    )

    # Reconstruct the covariance actually used by the simulator for X_t.
    generated_sigma <- matrix(
      simulated$sigma_mat[, t], nrow = resolution, ncol = resolution
    )
    generated_svd <- svd(generated_sigma)
    generated_sqrt <- sqrt(resolution) *
      generated_svd$u %*% diag(sqrt(generated_svd$d)) %*%
      t(generated_svd$v)
    generated_covariance <- generated_sqrt %*% C_eps %*%
      generated_sqrt / resolution^2
    prediction_difference[i] <- max(abs(oracle_sigma - generated_covariance))

    fitted <- fungarch:::.fit_diag_model_new(
      X = X_train,
      C_eps = C_eps,
      basis = basis,
      eigenvalues = eigenvalues,
      w = w,
      p = 1
    )

    fitted_sigma <- fungarch:::.pred_diag_new(
      resolution = resolution,
      alphas = fitted$alpha,
      basis = basis,
      delta = fitted$delta,
      data = X_train,
      C_eps = C_eps
    )

    for (a in seq_along(alpha_levels)) {
      z <- qnorm(alpha_levels[a])
      oracle_threshold <- z * sqrt(pmax(0, diag(oracle_sigma)))
      fitted_threshold <- z * sqrt(pmax(0, diag(fitted_sigma)))

      oracle_violations[i, a] <- mean(X_next <= oracle_threshold)
      fitted_violations[i, a] <- mean(X_next <= fitted_threshold)
    }
  }

  list(
    nominal = alpha_levels,
    oracle = colMeans(oracle_violations),
    fitted = colMeans(fitted_violations),
    max_prediction_difference = max(prediction_difference),
    oracle_by_forecast = oracle_violations,
    fitted_by_forecast = fitted_violations
  )
}

diagnostic <- simulate_coverage_diagnostic()
print(diagnostic$nominal)
print(diagnostic$max_prediction_difference)
print(diagnostic$oracle)
print(diagnostic$fitted)
