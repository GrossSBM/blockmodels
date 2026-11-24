set.seed(12)

test_that("BM_bernoulli_multiplex SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  P00 <- matrix(runif(Q * Q), Q, Q)
  P10 <- matrix(runif(Q * Q), Q, Q)
  P01 <- matrix(runif(Q * Q), Q, Q)
  P11 <- matrix(runif(Q * Q), Q, Q)
  SumP <- P00 + P10 + P01 + P11
  P00 <- P00 / SumP
  P01 <- P01 / SumP
  P10 <- P10 / SumP
  P11 <- P11 / SumP
  MU <- matrix(runif(n * n), n, n)
  M1 <- 1 * (MU > Z %*% (P00 + P01) %*% t(Z))
  M2 <- 1 * ((MU > Z %*% P00 %*% t(Z)) & (MU < Z %*% (P00 + P01 + P11) %*% t(Z)))

  model <- BM_bernoulli_multiplex("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_bernoulli_multiplex SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  P00 <- matrix(runif(Q * Q), Q, Q)
  P10 <- matrix(runif(Q * Q), Q, Q)
  P01 <- matrix(runif(Q * Q), Q, Q)
  P11 <- matrix(runif(Q * Q), Q, Q)
  SumP <- P00 + P10 + P01 + P11
  P00 <- P00 / SumP
  P01 <- P01 / SumP
  P10 <- P10 / SumP
  P11 <- P11 / SumP
  P00[lower.tri(P00)] <- t(P00)[lower.tri(P00)]
  P01[lower.tri(P01)] <- t(P01)[lower.tri(P01)]
  P10[lower.tri(P10)] <- t(P10)[lower.tri(P10)]
  P11[lower.tri(P11)] <- t(P11)[lower.tri(P11)]
  MU <- matrix(runif(n * n), n, n)
  MU[lower.tri(MU)] <- t(MU)[lower.tri(MU)]
  M1 <- 1 * (MU > Z %*% (P00 + P01) %*% t(Z))
  M2 <- 1 * ((MU > Z %*% P00 %*% t(Z)) & (MU < Z %*% (P00 + P01 + P11) %*% t(Z)))

  model <- BM_bernoulli_multiplex("SBM_sym", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_bernoulli_multiplex LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  P00 <- matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  P10 <- matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  P01 <- matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  P11 <- matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  SumP <- P00 + P10 + P01 + P11
  P00 <- P00 / SumP
  P01 <- P01 / SumP
  P10 <- P10 / SumP
  P11 <- P11 / SumP
  MU <- matrix(runif(n[1] * n[2]), n[1], n[2])
  M1 <- 1 * (MU > Z1 %*% (P00 + P01) %*% t(Z2))
  M2 <- 1 * ((MU > Z1 %*% P00 %*% t(Z2)) & (MU < Z1 %*% (P00 + P01 + P11) %*% t(Z2)))

  model <- BM_bernoulli_multiplex("LBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_bernoulli_multiplex handles combinations never observed", {
  # the second layer is included in the first one: the combination 01 never occurs
  npc <- 15
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  P <- matrix(c(.8, .2, .2, .6), Q, Q)
  M1 <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z))
  M2 <- M1 * (matrix(runif(n * n), n, n) < 0.5)
  sym <- function(M) {M[lower.tri(M)] <- t(M)[lower.tri(M)]; M}

  for (membership in c("SBM", "SBM_sym")) {
    layers <- if (membership == "SBM_sym") list(sym(M1), sym(M2)) else list(M1, M2)
    model <- BM_bernoulli_multiplex(membership, layers, plotting = "", explore_min = 2, explore_max = 2,
                                    ncores = 1, verbosity = 0)
    expect_model_estimation(model)
    expect_true(all(is.finite(model$ICL)))
    expect_true(all(model$model_parameters[[2]]$pi$`01` == 0))
  }
  # identical layers: neither 01 nor 10 occurs
  model <- BM_bernoulli_multiplex("SBM", list(M1, M1), plotting = "", explore_min = 2, explore_max = 2,
                                  ncores = 1, verbosity = 0)
  expect_model_estimation(model)
  expect_true(all(is.finite(model$ICL)))

  # LBM
  M1 <- 1 * (matrix(runif(20 * 30), 20, 30) < 0.4)
  M2 <- M1 * (matrix(runif(20 * 30), 20, 30) < 0.5)
  model <- BM_bernoulli_multiplex("LBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2,
                                  ncores = 1, verbosity = 0)
  expect_model_estimation(model)
  expect_true(all(is.finite(model$ICL[-1])))
})
