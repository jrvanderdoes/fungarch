#' Evaluate forecast tightness for competing functional ARCH models
#'
#' Fits diagonal ARCH(1), diagonal ARCH(5), pointwise ARCH(1), and historical
#' benchmark forecasts over an expanding training window. The returned
#' estimates can be used to summarize how often each forecast bound exceeds
#' the realized curve.
#'
#' @param data Numeric matrix of functional observations, with curves in
#'   columns.
#' @param N_fit Number of initial curves used for the first fit.
#' @param w Smoothing/regularization weight passed to the fitted models.
#' @param alp Quantile level for the historical benchmark and model bounds.
#' @return A list with components `estimates`, `mean_estimate`, and
#'   `viol_estimate`.
#' @export
check_model_tightness <- function(data, N_fit = 450, w = 0.1, alp = 0.05) {
  if (!is.numeric(data) || is.null(dim(data)) || nrow(data) < 2) {
    stop("data must be a numeric matrix with at least two rows.")
  }
  if (!is.numeric(N_fit) || length(N_fit) != 1 ||
      N_fit != as.integer(N_fit)) {
    stop("N_fit must be one integer.")
  }
  if (N_fit <= 5 || N_fit >= ncol(data)) {
    stop("N_fit must be greater than 5 and leave forecast observations.")
  }
  if (!is.numeric(alp) || length(alp) != 1 || alp <= 0 || alp >= 1) {
    stop("alp must be a single value strictly between 0 and 1.")
  }

  resolution <- nrow(data)
  forecast <- ncol(data) - N_fit
  C_eps <- .generate_ceps(resolution)
  eig_Ce <- eigen(C_eps)
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution
  estimates <- array(dim = c(resolution, 5, forecast))

  cat(paste0("(", forecast, ")"))
  for (n in N_fit:(N_fit + forecast - 1)) {
    cat(paste0(n - N_fit + 1, ", "))
    training_data <- data[, seq_len(n), drop = FALSE]
    sds <- apply(training_data, MARGIN = 1, sd)

    diag_mod <- .fit_diag_model_new(
      X = training_data, C_eps = C_eps, basis = basis,
      eigenvalues = eigenvalues, w = w, p = 1
    )
    diag_pred <- .pred_diag_new(
      resolution, alphas = diag_mod$alpha, basis = basis,
      delta = diag_mod$delta, data = training_data, C_eps = C_eps
    )

    diag_mod5 <- .fit_diag_model_new(
      X = training_data, C_eps = C_eps, basis = basis,
      eigenvalues = eigenvalues, w = w, p = 5
    )
    diag_pred5 <- .pred_diag_new(
      resolution, alphas = diag_mod5$alpha, basis = basis,
      delta = diag_mod5$delta, data = training_data, C_eps = C_eps
    )

    pw_mod <- .fit_pw(X = training_data, w = w)
    pw_pred1 <- .pred_pw_simple(
      X = training_data, pw_mod$alpha, pw_mod$delta, C_eps
    )
    hist_pred <- apply(training_data, MARGIN = 1, quantile, probs = alp)

    estimates[, , n - N_fit + 1] <- as.matrix(data.frame(
      data[, n + 1],
      sqrt(pmax(0, diag(diag_pred))) * qnorm(alp, sd = sds),
      sqrt(pmax(0, diag(diag_pred5))) * qnorm(alp, sd = sds),
      pw_pred1 * qnorm(alp, sd = sds),
      hist_pred
    ))
  }

  mean_estimate <- apply(estimates, c(1, 2), mean)
  viol_estimate <- data.frame(
    D = rowMeans(matrix(estimates[, 2, ], nrow = resolution) >=
                   matrix(estimates[, 1, ], nrow = resolution)),
    D5 = rowMeans(matrix(estimates[, 3, ], nrow = resolution) >=
                    matrix(estimates[, 1, ], nrow = resolution)),
    P = rowMeans(matrix(estimates[, 4, ], nrow = resolution) >=
                   matrix(estimates[, 1, ], nrow = resolution)),
    H = rowMeans(matrix(estimates[, 5, ], nrow = resolution) >=
                   matrix(estimates[, 1, ], nrow = resolution))
  )

  list(
    estimates = estimates,
    mean_estimate = mean_estimate,
    viol_estimate = viol_estimate
  )
}
