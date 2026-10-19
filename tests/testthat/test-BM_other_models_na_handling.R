
test_that("BM_bernoulli_multiplex SBM handles partial NA", {
  skip("To be implemented in the future.")
  set.seed(141)
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)

  P00 <- matrix(runif(Q * Q), Q, Q)
  P10 <- matrix(runif(Q * Q), Q, Q)
  P01 <- matrix(runif(Q * Q), Q, Q)
  P11 <- matrix(runif(Q * Q), Q, Q)
  S <- P00 + P10 + P01 + P11
  P00 <- P00 / S
  P01 <- P01 / S
  P10 <- P10 / S
  P11 <- P11 / S

  U <- matrix(runif(n * n), n, n)
  M1 <- 1 * (U > Z %*% (P00 + P01) %*% t(Z))
  M2 <- 1 * ((U > Z %*% P00 %*% t(Z)) & (U < Z %*% (P00 + P01 + P11) %*% t(Z)))

  off_diag <- which(row(M1) != col(M1))
  set.seed(142)
  idx <- sample(off_diag, max(1, floor(0.15 * length(off_diag))))
  M1[idx] <- NA_real_
  M2[idx] <- NA_real_

  set.seed(143)
  model <- BM_bernoulli_multiplex("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model$estimate())
})

test_that("BM_gaussian_multivariate SBM handles partial NA", {
  skip("To be implemented in the future.")
  set.seed(151)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 4 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 4 * matrix(runif(Q * Q), Q, Q)
  N1 <- matrix(rnorm(n * n, sd = 1), n, n)
  N2 <- matrix(rnorm(n * n, sd = 1), n, n)
  M1 <- Z %*% Mu1 %*% t(Z) + N1
  M2 <- Z %*% Mu2 %*% t(Z) + 8 * N1 + N2

  off_diag <- which(row(M1) != col(M1))
  set.seed(152)
  idx <- sample(off_diag, max(1, floor(0.15 * length(off_diag))))
  M1[idx] <- NA_real_
  M2[idx] <- NA_real_

  set.seed(153)
  model <- BM_gaussian_multivariate("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model$estimate())
})

test_that("BM_gaussian_multivariate_independent SBM handles partial NA", {
  skip("To be implemented in the future.")
  set.seed(161)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 6 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 6 * matrix(runif(Q * Q), Q, Q)
  M1 <- matrix(rnorm(n * n, sd = 4), n, n) + Z %*% Mu1 %*% t(Z)
  M2 <- matrix(rnorm(n * n, sd = 7), n, n) + Z %*% Mu2 %*% t(Z)

  off_diag <- which(row(M1) != col(M1))
  set.seed(162)
  idx <- sample(off_diag, max(1, floor(0.15 * length(off_diag))))
  M1[idx] <- NA_real_
  M2[idx] <- NA_real_

  set.seed(163)
  model <- BM_gaussian_multivariate_independent("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model$estimate())
})

test_that("BM_gaussian_multivariate_independent_homoscedastic SBM handles partial NA", {
  skip("To be implemented in the future.")
  set.seed(171)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 5 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 5 * matrix(runif(Q * Q), Q, Q)
  M1 <- matrix(rnorm(n * n, sd = 4), n, n) + Z %*% Mu1 %*% t(Z)
  M2 <- matrix(rnorm(n * n, sd = 4), n, n) + Z %*% Mu2 %*% t(Z)

  off_diag <- which(row(M1) != col(M1))
  set.seed(172)
  idx <- sample(off_diag, max(1, floor(0.15 * length(off_diag))))
  M1[idx] <- NA_real_
  M2[idx] <- NA_real_

  set.seed(173)
  model <- BM_gaussian_multivariate_independent_homoscedastic("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model$estimate())
})

test_that("Gaussian multivariate SBM models fail gracefully when all values are NA", {
  skip("To be implemented in the future.")
  n <- 20
  M1 <- matrix(NA_real_, n, n)
  M2 <- matrix(NA_real_, n, n)

  model_mv <- BM_gaussian_multivariate("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_error(model_mv$estimate(), "No valid edges in adjacency matrix. They are all NAs\\.")

  model_mvi <- BM_gaussian_multivariate_independent("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_error(model_mvi$estimate(), "No valid edges in adjacency matrix. They are all NAs\\.")

  model_mvih <- BM_gaussian_multivariate_independent_homoscedastic("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_error(model_mvih$estimate(), "No valid edges in adjacency matrix. They are all NAs\\.")
})
test_that("Multivariate models reject NA in the adjacency matrices with an explicit error", {
  n <- 10
  M1 <- matrix(rnorm(n * n), n, n)
  M2 <- matrix(rnorm(n * n), n, n)
  M1[3, 4] <- NA
  B1 <- 1 * (M1 > 0)
  B2 <- 1 * (M2 > 0)

  for (model_class in c("BM_gaussian_multivariate", "BM_gaussian_multivariate_independent", "BM_gaussian_multivariate_independent_homoscedastic")) {
    expect_error(
      get(model_class)("SBM", list(M1, M2), plotting = "", ncores = 1, verbosity = 0),
      "not supported for multivariate models"
    )
  }
  expect_error(
    BM_bernoulli_multiplex("SBM", list(B1, B2), plotting = "", ncores = 1, verbosity = 0),
    "not supported for multivariate models"
  )
})

test_that("Models with covariates reject NA in the covariates with an explicit error", {
  n <- 10
  M <- 1 * (matrix(runif(n * n), n, n) < 0.3)
  C <- matrix(rnorm(n * n), n, n)
  C[2, 7] <- NA

  for (model_class in c("BM_bernoulli_covariates", "BM_bernoulli_covariates_fast", "BM_poisson_covariates", "BM_gaussian_covariates")) {
    expect_error(
      get(model_class)("SBM", M, list(C), plotting = "", ncores = 1, verbosity = 0),
      "contains missing values"
    )
  }
})
