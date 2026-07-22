set.seed(12)

test_that("BM_gaussian_multivariate SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 4 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 4 * matrix(runif(Q * Q), Q, Q)
  Noise1 <- matrix(rnorm(n * n, sd = 1), n, n)
  Noise2 <- matrix(rnorm(n * n, sd = 1), n, n)
  M1 <- Z %*% Mu1 %*% t(Z) + Noise1
  M2 <- Z %*% Mu2 %*% t(Z) + 10 * Noise1 + Noise2

  model <- BM_gaussian_multivariate("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian_multivariate SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 4 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 4 * matrix(runif(Q * Q), Q, Q)
  Noise1 <- matrix(rnorm(n * n, sd = 1), n, n)
  Noise2 <- matrix(rnorm(n * n, sd = 1), n, n)
  M1 <- Z %*% Mu1 %*% t(Z) + Noise1
  M2 <- Z %*% Mu2 %*% t(Z) + 10 * Noise1 + Noise2
  M1[lower.tri(M1)] <- t(M1)[lower.tri(M1)]
  M2[lower.tri(M2)] <- t(M2)[lower.tri(M2)]

  model <- BM_gaussian_multivariate("SBM_sym", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian_multivariate LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  Mu1 <- 4 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  Mu2 <- 4 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  Noise1 <- matrix(rnorm(n[1] * n[2], sd = 1), n[1], n[2])
  Noise2 <- matrix(rnorm(n[1] * n[2], sd = 1), n[1], n[2])
  M1 <- Z1 %*% Mu1 %*% t(Z2) + Noise1
  M2 <- Z1 %*% Mu2 %*% t(Z2) + 10 * Noise1 + Noise2

  model <- BM_gaussian_multivariate("LBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 3, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})
