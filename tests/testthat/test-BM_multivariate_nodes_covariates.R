set.seed(20260228)

test_that("multivariate SBM split keeps nodes covariates", {
  n <- 12
  Q <- 2
  npc <- n / Q

  Z <- diag(Q) %x% matrix(1, npc, 1)
  M1 <- matrix(rnorm(n * n), n, n)
  M2 <- matrix(rnorm(n * n), n, n)
  X <- cbind(1, runif(n), runif(n))

  model <- BM_gaussian_multivariate(
    "SBM",
    list(M1, M2),
    nodes_covariates = list(node = X),
    plotting = "",
    verbosity = 0
  )

  model$memberships[[Q]] <- methods::getRefClass("SBM")(from_cc = list(Z = Z), nodes_covar = X)
  model$model_parameters[[Q]] <- list(
    mu = list(
      matrix(0, nrow = Q, ncol = Q),
      matrix(0, nrow = Q, ncol = Q)
    ),
    Sigma = diag(2)
  )
  splits <- model$split_membership_model(Q)

  expect_gt(length(splits), 0)
  expect_true(all(vapply(splits, function(m) is.matrix(m$nodes_covariates), logical(1))))
  expect_true(all(vapply(splits, function(m) identical(m$nodes_covariates, X), logical(1))))
})

test_that("multivariate LBM split keeps row/col covariates", {
  n1 <- 10
  n2 <- 8
  Q1 <- 2
  Q2 <- 2
  Q <- Q1 + Q2

  Z1 <- diag(Q1) %x% matrix(1, n1 / Q1, 1)
  Z2 <- diag(Q2) %x% matrix(1, n2 / Q2, 1)
  M1 <- matrix(rnorm(n1 * n2), n1, n2)
  M2 <- matrix(rnorm(n1 * n2), n1, n2)
  X_row <- cbind(1, runif(n1), runif(n1))
  X_col <- cbind(1, runif(n2))

  model <- BM_gaussian_multivariate(
    "LBM",
    list(M1, M2),
    nodes_covariates = list(row = X_row, col = X_col),
    plotting = "",
    verbosity = 0
  )

  model$memberships[[Q]] <- methods::getRefClass("LBM")(from_cc = list(Z1 = Z1, Z2 = Z2), row_covar = X_row, col_covar = X_col)
  model$model_parameters[[Q]] <- list(
    mu = list(
      matrix(0, nrow = Q1, ncol = Q2),
      matrix(0, nrow = Q1, ncol = Q2)
    ),
    Sigma = diag(2)
  )
  splits <- model$split_membership_model(Q)

  expect_gt(length(splits), 0)
  expect_true(all(vapply(splits, function(m) is.matrix(m$row_covariates), logical(1))))
  expect_true(all(vapply(splits, function(m) is.matrix(m$col_covariates), logical(1))))
  expect_true(all(vapply(splits, function(m) identical(m$row_covariates, X_row), logical(1))))
  expect_true(all(vapply(splits, function(m) identical(m$col_covariates, X_col), logical(1))))
})
