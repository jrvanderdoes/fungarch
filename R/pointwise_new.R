#' Fit a multi-lag pointwise ARCH model
#'
#' Estimate pointwise autoregressive operators for a squared functional
#' process using one coefficient matrix per lag.
#'
#' @param X Finite numeric matrix of functional observations, with observations
#'   in columns.
#' @param p Positive integer number of autoregressive lags.
#' @param w Truncation threshold for the covariance inverse when
#'   `inv_type = "M"`.
#' @param inv_type Inversion method: `"M"` for a truncated eigendecomposition
#'   or `"T"` for Tikhonov regularization.
#' @param lambda Positive Tikhonov parameter used when `inv_type = "T"`.
#'
#' @return A list containing `alpha`, a list of lag coefficient matrices,
#'   `delta`, pointwise innovation variances, `mean`, the training mean used
#'   for centering, and `p`.
#'
#' @noRd
.fit_pw_new <- function(X, p = 1, w = 0.01,
                        inv_type = c("M", "T"), lambda = 1e-6){
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix", call. = FALSE)
  }
  if (nrow(X) < 1L) {
    stop("X must have at least one functional coordinate", call. = FALSE)
  }
  if (length(p) != 1L || !is.numeric(p) || !is.finite(p) ||
      p < 1 || p != floor(p)) {
    stop("p must be a single positive integer", call. = FALSE)
  }
  p <- as.integer(p)
  if (ncol(X) <= p + 1L) {
    stop("X must contain more than p + 1 observations", call. = FALSE)
  }
  if (length(w) != 1L || !is.numeric(w) || !is.finite(w) || w < 0 || w > 1) {
    stop("w must be a single finite value between 0 and 1", call. = FALSE)
  }
  inv_type <- match.arg(inv_type, c("M", "T"))
  if (length(lambda) != 1L || !is.numeric(lambda) ||
      !is.finite(lambda) || lambda <= 0) {
    stop("lambda must be a single positive finite value", call. = FALSE)
  }

  n <- nrow(X)
  N <- ncol(X)
  training_mean <- rowMeans(X)
  Y <- sweep(X, 1L, training_mean, "-")^2
  n_pairs <- N - p

  predictors <- do.call(rbind, lapply(seq_len(p), function(lag) {
    Y[, (p + 1L - lag):(N - lag), drop = FALSE]
  }))
  response <- Y[, (p + 1L):N, drop = FALSE]

  predictor_centered <- sweep(predictors, 1L, rowMeans(predictors), "-")
  response_centered <- sweep(response, 1L, rowMeans(response), "-")
  predictor_cov <- predictor_centered %*% t(predictor_centered) /
    (n_pairs - 1)
  response_predictor_cov <- response_centered %*% t(predictor_centered) /
    (n_pairs - 1)

  inverse <- .pw_new_inverse(predictor_cov, w, inv_type, lambda)
  alpha_block <- response_predictor_cov %*% inverse / n
  alpha <- lapply(seq_len(p), function(lag) {
    columns <- ((lag - 1L) * n + 1L):(lag * n)
    alpha_block[, columns, drop = FALSE]
  })

  mean_squared <- rowMeans(Y)
  fitted_mean <- Reduce("+", lapply(alpha, function(operator) {
    as.vector(operator %*% mean_squared)
  })) / n

  list(
    alpha = alpha,
    delta = pmax(0, mean_squared - fitted_mean),
    mean = training_mean,
    p = p
  )
}

#' Compute a regularized inverse for a pointwise lagged covariance matrix.
#'
#' @noRd
.pw_new_inverse <- function(matrix, w, inv_type, lambda){
  matrix <- (matrix + t(matrix)) / 2

  if (inv_type == "T") {
    return(solve(matrix + diag(lambda, nrow(matrix), ncol(matrix))))
  }

  decomposition <- eigen(matrix, symmetric = TRUE)
  values <- pmax(decomposition$values, 0)
  total <- sum(values)
  if (!is.finite(total) || total <= 0) {
    stop("lagged covariance has no positive eigenvalues", call. = FALSE)
  }

  cutoff <- which(cumsum(values) / total >= 1 - w)[1]
  tolerance <- .Machine$double.eps * max(1, max(values))
  keep <- seq_len(cutoff)
  keep <- keep[values[keep] > tolerance]
  if (!length(keep)) {
    stop("lagged covariance has no numerically usable eigenvalues", call. = FALSE)
  }

  vectors <- decomposition$vectors[, keep, drop = FALSE]
  sweep(vectors, 2L, 1 / values[keep], "*") %*% t(vectors)
}

#' Predict a multi-lag pointwise conditional covariance.
#'
#' @param X Recent functional observations, with observations in columns.
#' @param alpha List of lag coefficient matrices from `.fit_pw_new()`.
#' @param delta Pointwise innovation variances.
#' @param C_eps Innovation covariance matrix.
#' @param mean Training mean used to center the observations. Defaults to zero.
#'
#' @return A predicted conditional covariance matrix.
#'
#' @noRd
.pred_pw_new <- function(X, alpha, delta, C_eps, mean = NULL){
  estimate <- .pointwise_variance_new(X, alpha, delta, mean)
  sqrt(estimate) %*% t(sqrt(estimate)) * C_eps
}

#' Predict multi-lag pointwise conditional standard deviations.
#'
#' @inheritParams .pred_pw_new
#'
#' @return A vector of predicted pointwise standard deviations.
#'
#' @noRd
.pred_pw_new_simple <- function(X, alpha, delta, C_eps = NULL, mean = NULL){
  sqrt(.pointwise_variance_new(X, alpha, delta, mean))
}

#' Compute multi-lag pointwise conditional variances.
#'
#' @noRd
.pointwise_variance_new <- function(X, alpha, delta, mean = NULL){
  if (!is.list(alpha) || length(alpha) < 1L) {
    stop("alpha must be a nonempty list of lag matrices", call. = FALSE)
  }
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix", call. = FALSE)
  }
  if (is.null(mean)) {
    mean <- rep(0, nrow(X))
  }
  if (length(mean) != nrow(X) || any(!is.finite(mean))) {
    stop("mean must have one finite value per row of X", call. = FALSE)
  }
  if (length(delta) != nrow(X) || any(!is.finite(delta))) {
    stop("delta must have one finite value per row of X", call. = FALSE)
  }
  if (ncol(X) < length(alpha)) {
    stop("X must contain at least one observation for each lag", call. = FALSE)
  }

  X2 <- sweep(X, 1L, mean, "-")^2
  estimate <- as.numeric(delta)
  for (lag in seq_along(alpha)) {
    operator <- alpha[[lag]]
    if (!is.matrix(operator) || nrow(operator) != nrow(X) ||
        ncol(operator) != nrow(X)) {
      stop("each alpha entry must be a square matrix matching X", call. = FALSE)
    }
    previous <- X2[, ncol(X2) - lag + 1L]
    estimate <- estimate + vapply(seq_len(nrow(X2)), function(j) {
      fChange:::dot_integrate(as.matrix(operator[j, ] * previous))
    }, numeric(1))
  }

  pmax(estimate, 0)
}
