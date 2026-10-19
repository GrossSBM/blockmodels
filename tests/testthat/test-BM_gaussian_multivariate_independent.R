set.seed(12)

test_that("BM_gaussian_multivariate_independent SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 8 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 8 * matrix(runif(Q * Q), Q, Q)
  M1 <- matrix(rnorm(n * n, sd = 5), n, n) + Z %*% Mu1 %*% t(Z)
  M2 <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu2 %*% t(Z)

  model <- BM_gaussian_multivariate_independent("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian_multivariate_independent SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 8 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 8 * matrix(runif(Q * Q), Q, Q)
  Mu1[lower.tri(Mu1)] <- t(Mu1)[lower.tri(Mu1)]
  Mu2[lower.tri(Mu2)] <- t(Mu2)[lower.tri(Mu2)]
  M1 <- matrix(rnorm(n * n, sd = 5), n, n) + Z %*% Mu1 %*% t(Z)
  M2 <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu2 %*% t(Z)
  M1[lower.tri(M1)] <- t(M1)[lower.tri(M1)]
  M2[lower.tri(M2)] <- t(M2)[lower.tri(M2)]

  model <- BM_gaussian_multivariate_independent("SBM_sym", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian_multivariate_independent LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  Mu1 <- 8 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  Mu2 <- 8 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  M1 <- matrix(rnorm(n[1] * n[2], sd = 5), n[1], n[2]) + Z1 %*% Mu1 %*% t(Z2)
  M2 <- matrix(rnorm(n[1] * n[2], sd = 10), n[1], n[2]) + Z1 %*% Mu2 %*% t(Z2)

  model <- BM_gaussian_multivariate_independent("LBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})
