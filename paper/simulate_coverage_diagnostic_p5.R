# Diagnostic size check for the p = 5 diagonal FARCH model.
#
# This is intentionally separate from simulate_coverage_diagnostic.R so the
# validated p = 1 experiment is preserved unchanged.

simulate_coverage_diagnostic_p5 <- function(
    N = 3000L,
    resolution = 20L,
    burnin = 100L,
    training = 1500L,
    alpha_levels = c(0.05, 0.01),
    lag_weights = c(0.30, 0.14, 0.08, 0.05, 0.03),
    w = 0.10,
    seed = 12345L,
    forecast_step = 10L) {

  if (!requireNamespace("fungarch", quietly = TRUE)) {
    stop("Install fungarch before running this diagnostic.")
  }

  p <- 5L

  stopifnot(length(lag_weights) == p)
  if (sum(lag_weights) >= 1) {
    stop("The lag weights must sum to less than one for this diagnostic.")
  }
  if (training <= burnin + p || training >= N) {
    stop("Use burnin + p < training < N.")
  }

  set.seed(seed)

  C_eps <- fungarch:::.generate_ceps(resolution)
  eig_Ce <- eigen(C_eps, symmetric = TRUE)
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  h <- lapply(seq_len(resolution), function(k) {
    hkk <- as.vector(basis[, k] %o% basis[, k])
    hkk %o% hkk
  })

  # Five temporal lags, all acting in the first spatial basis direction.
  alpha <- lapply(lag_weights, function(weight) weight * h[[1L]])
  alpha_diag <- matrix(lag_weights, nrow = 1L, ncol = p)

  simulated <- fungarch:::generateFARCH_use(
    N = N,
    resolution = resolution,
    alpha = alpha,
    delta = C_eps,
    basis = basis,
    burnin = burnin,
    p = p,
    eigs = eigenvalues[1L],
    C_eps = C_eps,
    kl_errors = FALSE
  )

  X <- simulated$data
  forecast_indices <- seq.int(training + 1L, N, by = forecast_step)

  # generateFARCH_use() may retain burn-in covariance entries even when the
  # returned data have already had burn-in removed. Align the final N
  # covariance entries with the N returned observations.
  sigma_n <- ncol(simulated$sigma_mat)
  data_n <- ncol(X)
  sigma_offset <- sigma_n - data_n
  if (sigma_offset < 0L) {
    stop("sigma_mat has fewer time points than simulated$data.")
  }

  # Estimate once on the training sample. The fitted parameters remain fixed
  # while predictions use the observed history through each forecast origin.
  fit <- fungarch:::.fit_diag_model_new(
    X = X[, seq_len(training), drop = FALSE],
    C_eps = C_eps,
    basis = basis,
    eigenvalues = eigenvalues,
    w = w,
    p = p
  )

  get_diag_variance <- function(sigma) {
    if (is.matrix(sigma)) {
      return(diag(sigma))
    }
    if (length(sigma) == resolution^2) {
      return(diag(matrix(sigma, nrow = resolution)))
    }
    if (length(sigma) == resolution) {
      return(as.numeric(sigma))
    }
    stop("Unexpected sigma_mat shape.")
  }

  get_true_variance <- function(index) {
    # Use the same prediction route as the validated p = 1 diagnostic, but
    # supply the known simulation parameters rather than fitted parameters.
    # This is the oracle conditional covariance for the p = 5 process.
    oracle_prediction <- fungarch:::.pred_diag_new(
      data = X[, seq_len(index - 1L), drop = FALSE],
      resolution = resolution,
      C_eps = C_eps,
      basis = basis,
      alpha = alpha_diag,
      delta = C_eps
    )
    get_diag_variance(oracle_prediction)
  }

  oracle_rates <- matrix(
    0,
    nrow = length(forecast_indices),
    ncol = length(alpha_levels),
    dimnames = list(NULL, paste0("alpha_", alpha_levels))
  )
  fitted_rates <- oracle_rates
  fitted_variance_differences <- numeric(length(forecast_indices))
  fitted_variance_ratios <- numeric(length(forecast_indices))

  for (j in seq_along(forecast_indices)) {
    index <- forecast_indices[[j]]
    variance <- get_true_variance(index)
    oracle_thresholds <- vapply(
      alpha_levels,
      function(level) qnorm(level) * sqrt(pmax(variance, 0)),
      numeric(resolution)
    )
    future_value <- X[, index]

    oracle_rates[j, ] <- vapply(
      seq_along(alpha_levels),
      function(k) mean(future_value < oracle_thresholds[, k]),
      numeric(1)
    )

    fitted_prediction <- fungarch:::.pred_diag_new(
      data = X[, seq_len(index - 1L), drop = FALSE],
      resolution = resolution,
      C_eps = C_eps,
      basis = basis,
      alpha = fit$alpha,
      delta = fit$delta
    )

    fitted_variance <- if (is.matrix(fitted_prediction)) {
      diag(fitted_prediction)
    } else if (length(fitted_prediction) == resolution^2) {
      diag(matrix(fitted_prediction, nrow = resolution))
    } else {
      as.numeric(fitted_prediction)
    }

    fitted_thresholds <- vapply(
      alpha_levels,
      function(level) qnorm(level) * sqrt(pmax(fitted_variance, 0)),
      numeric(resolution)
    )
    fitted_rates[j, ] <- vapply(
      seq_along(alpha_levels),
      function(k) mean(future_value < fitted_thresholds[, k]),
      numeric(1)
    )

    fitted_variance_differences[[j]] <- max(abs(fitted_variance - variance))
    fitted_variance_ratios[[j]] <- mean(
      fitted_variance / pmax(variance, .Machine$double.eps)
    )
  }

  list(
    p = p,
    nominal = alpha_levels,
    oracle = colMeans(oracle_rates),
    fitted = colMeans(fitted_rates),
    max_fitted_variance_difference = max(fitted_variance_differences),
    mean_fitted_variance_ratio = mean(fitted_variance_ratios),
    median_fitted_variance_ratio = median(fitted_variance_ratios),
    fitted_dimension = fit$K,
    true_alpha = alpha_diag,
    fitted_alpha = fit$alpha,
    oracle_by_forecast = oracle_rates,
    fitted_by_forecast = fitted_rates,
    settings = list(
      N = N,
      resolution = resolution,
      burnin = burnin,
      training = training,
      lag_weights = lag_weights,
      w = w,
      seed = seed,
      forecast_step = forecast_step,
      sigma_offset = sigma_offset
    )
  )
}
