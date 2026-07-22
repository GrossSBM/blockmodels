set.seed(12)

test_that("BM_poisson_covariates SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  L <- 70 * matrix(runif(Q * Q), Q, Q)
  M_in_expectation_without_covariates <- Z %*% L %*% t(Z)
  Y1 <- matrix(runif(n * n), n, n)
  Y2 <- matrix(runif(n * n), n, n)
  M_in_expectation <- M_in_expectation_without_covariates * exp(4.2 * Y1 - 1.2 * Y2)
  M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)

  model <- BM_poisson_covariates("SBM", M, list(Y1, Y2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_poisson_covariates SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  L <- 70 * matrix(runif(Q * Q), Q, Q)
  L[lower.tri(L)] <- t(L)[lower.tri(L)]
  M_in_expectation_without_covariates <- Z %*% L %*% t(Z)
  Y1 <- matrix(runif(n * n), n, n)
  Y2 <- matrix(runif(n * n), n, n)
  Y1[lower.tri(Y1)] <- t(Y1)[lower.tri(Y1)]
  Y2[lower.tri(Y2)] <- t(Y2)[lower.tri(Y2)]
  M_in_expectation <- M_in_expectation_without_covariates * exp(4.2 * Y1 - 1.2 * Y2)
  M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)
  M[lower.tri(M)] <- t(M)[lower.tri(M)]

  model <- BM_poisson_covariates("SBM_sym", M, list(Y1, Y2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_poisson_covariates LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  L <- 70 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  M_in_expectation_without_covariates <- Z1 %*% L %*% t(Z2)
  Y1 <- matrix(runif(n[1] * n[2]), n[1], n[2])
  Y2 <- matrix(runif(n[1] * n[2]), n[1], n[2])
  M_in_expectation <- M_in_expectation_without_covariates * exp(4.2 * Y1 - 1.2 * Y2)
  M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n[1], n[2])

  model <- BM_poisson_covariates("LBM", M, list(Y1, Y2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_poisson_covariates SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(121)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  L <- 50 * matrix(runif(Q * Q), Q, Q)
  Y1 <- matrix(runif(n * n), n, n)
  Y2 <- matrix(runif(n * n), n, n)
  E <- Z %*% L %*% t(Z) * exp(2.2 * Y1 - 0.8 * Y2)
  M <- matrix(rpois(length(as.vector(E)), as.vector(E)), n, n)

  off_diag <- which(row(M) != col(M))
  set.seed(122)
  M[sample(off_diag, max(1, floor(0.15 * length(off_diag))))] <- NA_real_

  set.seed(123)
  model_default <- BM_poisson_covariates("SBM", M, list(Y1, Y2), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(123)
  model_custom <- BM_poisson_covariates("SBM", M, list(Y1, Y2), na_replace_value = 1234, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$lambda, model_custom$model_parameters[[k]]$lambda, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$beta, model_custom$model_parameters[[k]]$beta, tolerance = 1e-10)
})
