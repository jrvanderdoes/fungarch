test_that("generate_ceps returns a symmetric OU covariance matrix", {
  C_eps <- fungarch:::.generate_ceps(5)

  expect_equal(dim(C_eps), c(5L, 5L))
  expect_equal(C_eps, t(C_eps))
  expect_equal(diag(C_eps), rep(1, 5))
  expect_true(all(C_eps >= 0 & C_eps <= 1))
})

test_that("generate_ceps rejects invalid inputs", {
  expect_error(fungarch:::.generate_ceps(0), "positive integer")
  expect_error(fungarch:::.generate_ceps(1.5), "positive integer")
  expect_error(fungarch:::.generate_ceps(5, denom = 0), "positive finite")
  expect_error(fungarch:::.generate_ceps(5, denom = -1), "positive finite")
  expect_error(fungarch:::.generate_ceps(5, denom = Inf), "positive finite")
})

test_that("estim_lag_cov computes lagged covariance", {
  X <- rbind(c(1, 2, 3), c(2, 4, 6))

  expect_equal(
    fungarch:::.estim_lag_cov(X),
    matrix(c(1, 2, 2, 4), nrow = 2)
  )
  expect_equal(
    fungarch:::.estim_lag_cov(X, lag = 1),
    matrix(c(0.5, 1, 1, 2), nrow = 2)
  )
})

test_that("estim_lag_cov rejects invalid lags and inputs", {
  X <- matrix(1:6, nrow = 2)

  expect_error(fungarch:::.estim_lag_cov(X, lag = -1), "nonnegative integer")
  expect_error(fungarch:::.estim_lag_cov(X, lag = 1.5), "nonnegative integer")
  expect_error(fungarch:::.estim_lag_cov(X, lag = 2), "at least two observations")
  expect_error(fungarch:::.estim_lag_cov(as.data.frame(X)), "numeric matrix")
})

test_that("CD_hat returns the matrices consumed by the diagonal fit", {
  skip("Legacy .CD_hat() is currently disabled")
  set.seed(1)
  X <- matrix(rnorm(4 * 30), nrow = 4)
  basis <- diag(sqrt(4), 4)
  C_eps_dag <- diag(4)

  result <- fungarch:::.CD_hat(
    X,
    dimen = 2,
    C_eps_dag = C_eps_dag,
    basis = basis
  )

  expect_named(result, c("C_hat", "D_hat", "m_X2"))
  expect_equal(dim(result$C_hat), c(2L, 2L))
  expect_equal(dim(result$D_hat), c(2L, 2L))
  expect_equal(dim(result$m_X2), c(4L, 4L))
})

test_that("fit_diag_model completes and exposes diagonal coefficients", {
  skip("Legacy .fit_diag_model() is currently disabled")
  set.seed(2)
  X <- matrix(rnorm(4 * 60), nrow = 4)
  C_eps <- fungarch:::.generate_ceps(4)
  eig_Ce <- eigen(C_eps)
  basis <- eig_Ce$vectors * sqrt(4)
  eigenvalues <- eig_Ce$values / 4

  result <- fungarch:::.fit_diag_model(
    X = X,
    C_eps = C_eps,
    basis = basis,
    eigenvalues = eigenvalues,
    dimen = 2
  )

  expect_length(result$alpha, 2)
  expect_equal(dim(result$delta), c(4L, 4L))
  expect_true(is.numeric(result$K_Ce))
})

