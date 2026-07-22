set.seed(12)

test_that("BM_gaussian_covariates SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu <- 20 * matrix(runif(Q * Q), Q, Q)
  Y1 <- matrix(runif(n * n), n, n)
  Y2 <- matrix(runif(n * n), n, n)
  M <- matrix(rnorm(n * n, sd = 5), n, n) + Z %*% Mu %*% t(Z) + 4.2 * Y1 - 1.6 * Y2

  model <- BM_gaussian_covariates("SBM", M, list(Y1, Y2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian_covariates SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu <- 20 * matrix(runif(Q * Q), Q, Q)
  Mu[lower.tri(Mu)] <- t(Mu)[lower.tri(Mu)]
  Y1 <- matrix(runif(n * n), n, n)
  Y2 <- matrix(runif(n * n), n, n)
  Y1[lower.tri(Y1)] <- t(Y1)[lower.tri(Y1)]
  Y2[lower.tri(Y2)] <- t(Y2)[lower.tri(Y2)]
  M <- matrix(rnorm(n * n, sd = 5), n, n) + Z %*% Mu %*% t(Z) + 4.2 * Y1 - 1.6 * Y2
  M[lower.tri(M)] <- t(M)[lower.tri(M)]

  model <- BM_gaussian_covariates("SBM_sym", M, list(Y1, Y2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian_covariates LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  Mu <- 20 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  Y1 <- matrix(runif(n[1] * n[2]), n[1], n[2])
  Y2 <- matrix(runif(n[1] * n[2]), n[1], n[2])
  M <- matrix(rnorm(n[1] * n[2], sd = 5), n[1], n[2]) + Z1 %*% Mu %*% t(Z2) + 4.2 * Y1 - 1.6 * Y2

  model <- BM_gaussian_covariates("LBM", M, list(Y1, Y2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian_covariates SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(131)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu <- 15 * matrix(runif(Q * Q), Q, Q)
  Y1 <- matrix(runif(n * n), n, n)
  M <- matrix(rnorm(n * n, sd = 5), n, n) + Z %*% Mu %*% t(Z) + 2.2 * Y1

  off_diag <- which(row(M) != col(M))
  set.seed(132)
  M[sample(off_diag, max(1, floor(0.05 * length(off_diag))))] <- NA_real_

  set.seed(133)
  model_default <- BM_gaussian_covariates("SBM", M, list(Y1), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(133)
  model_custom <- BM_gaussian_covariates("SBM", M, list(Y1), na_replace_value = -123, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$mu, model_custom$model_parameters[[k]]$mu, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$beta, model_custom$model_parameters[[k]]$beta, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$sigma2, model_custom$model_parameters[[k]]$sigma2, tolerance = 1e-10)
})