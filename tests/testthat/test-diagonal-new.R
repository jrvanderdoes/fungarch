test_that("diagonal p > 1 keeps lag predictors distinct in one dimension", {
  X <- matrix(1:4, nrow = 1L)
  out <- fungarch:::.CD_hat_new(
    X = X,
    p = 2L,
    dimen = 1L,
    C_eps_dag = matrix(1, nrow = 1L, ncol = 1L),
    basis = matrix(1, nrow = 1L, ncol = 1L)
  )

  # After centering X^2, the two lag columns are different. The old
  # X_Mp[, i] <- M[i] assignment recycled one value into both columns,
  # producing a symmetric C_hat with identical diagonal entries.
  expect_equal(
    out$C_hat[1L, 2L],
    15.125,
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(out$C_hat[1L, 1L], out$C_hat[2L, 2L])))
})
