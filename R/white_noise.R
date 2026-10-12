#' Standardize functional residuals from a fitted diagonal model
#'
#' Reconstructs the model covariance at each time point, projects the curves
#' onto the fitted basis, and standardizes the resulting score vectors.
#'
#' @param resolution Number of grid points in each curve.
#' @param alphas Matrix of fitted autoregressive coefficient matrices.
#' @param basis Basis matrix used by the fitted model.
#' @param delta Innovation covariance matrix in vectorized form.
#' @param data Numeric matrix of observed curves, with curves in columns.
#' @param C_eps Innovation covariance operator used by the fitted model.
#' @param dimen Number of basis scores retained by the fitted model.
#' @return A list containing `e_hat`, the matrix of standardized functional
#'   residual curves.
#' @export
standardize_functional_residuals <- function(
    resolution, alphas, basis, delta, data, C_eps, dimen) {
  if (!is.numeric(data) || is.null(dim(data)) || nrow(data) != resolution) {
    stop("data must be a numeric matrix with resolution rows.")
  }
  if (!is.matrix(alphas) || !is.matrix(basis)) {
    stop("alphas and basis must be matrices.")
  }
  if (nrow(basis) != resolution || dimen < 1 || dimen > ncol(basis)) {
    stop("basis and dimen are incompatible with resolution.")
  }
  if (length(delta) != resolution^2) {
    stop("delta must have length resolution^2.")
  }

  N <- ncol(data)
  lag <- ncol(alphas)
  if (lag < 1 || N <= lag) {
    stop("data must contain more curves than the fitted lag order.")
  }

  M <- matrix(0, nrow = dimen, ncol = N)
  for (i in seq_len(N)) {
    for (j in seq_len(dimen)) {
      M[j, i] <- (data[, i] %*% basis[, j]) / resolution
    }
  }

  alpha_applied <- .build_white_noise_alpha_operators(
    alphas = alphas, basis = basis, resolution = resolution
  )

  sigma <- matrix(ncol = N, nrow = length(delta))
  for (i in seq.int(lag + 1, N)) {
    sigma[, i] <- as.vector(delta) + Reduce(
      "+",
      lapply(seq_len(ncol(alphas)), function(idx) {
        operator(
          alpha_applied[[idx]],
          tcrossprod(data[, i - idx], data[, i - idx])
        )
      })
    )
  }

  errors <- matrix(nrow = dimen, ncol = N - lag)
  for (i in seq.int(lag + 1, N)) {
    sigmaP <- matrix(0, ncol = dimen, nrow = dimen)
    for (j in seq_len(dimen)) {
      for (k in seq_len(dimen)) {
        sigmaP[j, k] <- sum(
          matrix(sigma[, i], resolution, resolution) *
            (basis[, j] %o% basis[, k])
        ) / resolution^2
      }
    }

    sigma_svd <- svd(sigmaP)
    sigma_sqrt <- sigma_svd$u %*%
      diag(sqrt(pmax(sigma_svd$d, 0))) %*% t(sigma_svd$v)
    errors[, i - lag] <- solve(sigma_sqrt, M[, i])
  }

  e_hat <- matrix(0, nrow = resolution, ncol = ncol(errors))
  for (i in seq_len(ncol(errors))) {
    for (j in seq_len(dimen)) {
      e_hat[, i] <- e_hat[, i] + errors[j, i] * basis[, j]
    }
  }

  list(e_hat = e_hat)
}

#' Build lag-specific diagonal covariance operators for white-noise checks.
#'
#' @noRd
.build_white_noise_alpha_operators <- function(alphas, basis, resolution) {
  lapply(seq_len(ncol(alphas)), function(p) {
    result <- matrix(0, nrow = resolution^2, ncol = resolution^2)
    for (k in seq_len(nrow(alphas))) {
      hkk <- as.vector(basis[, k] %o% basis[, k])
      result <- result + alphas[k, p] * (hkk %o% hkk)
    }
    result
  })
}

#' Compute white-noise diagnostic residuals
#'
#' Fits the diagonal functional ARCH models of orders one and five, then
#' returns reconstructed curves and standardized residual curves for both
#' fitted models.
#'
#' @param data Numeric matrix of observed curves, with curves in columns.
#' @param w Smoothing/regularization weight passed to the fitted models.
#' @return A list with components `x_hat`, `e_hat`, and `e_hat5`.
#' @export
compute_white_noise_residuals <- function(data, w = 0.05) {
  if (!is.numeric(data) || is.null(dim(data)) || nrow(data) < 2 ||
      ncol(data) < 6) {
    stop("data must be a numeric matrix with at least 2 rows and 6 columns.")
  }

  resolution <- nrow(data)
  C_eps <- fChange::var(fChange::dfts(data))
  eig_Ce <- eigen(C_eps)
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  diag_mod <- .fit_diag_model_new(
    X = data, C_eps = C_eps, basis = basis,
    eigenvalues = eigenvalues, w = w, p = 1
  )
  errors_tmp <- standardize_functional_residuals(
    resolution = resolution, alphas = diag_mod$alpha,
    basis = basis, delta = diag_mod$delta, data = data,
    C_eps = C_eps, dimen = diag_mod$K
  )

  diag_mod5 <- .fit_diag_model_new(
    X = data, C_eps = C_eps, basis = basis,
    eigenvalues = eigenvalues, w = w, p = 5
  )
  errors5_tmp <- standardize_functional_residuals(
    resolution = resolution, alphas = diag_mod5$alpha,
    basis = basis, delta = diag_mod5$delta, data = data,
    C_eps = C_eps, dimen = diag_mod5$K
  )

  M <- matrix(0, nrow = resolution, ncol = ncol(data))
  for (i in seq_len(ncol(data))) {
    for (j in seq_len(resolution)) {
      M[j, i] <- (data[, i] %*% basis[, j]) / resolution
    }
  }

  x_hat <- matrix(0, nrow = resolution, ncol = ncol(data))
  for (i in seq_len(ncol(data))) {
    for (j in seq_len(resolution)) {
      x_hat[, i] <- x_hat[, i] + M[j, i] * basis[, j]
    }
  }

  list(
    x_hat = x_hat,
    e_hat = errors_tmp$e_hat,
    e_hat5 = errors5_tmp$e_hat
  )
}
