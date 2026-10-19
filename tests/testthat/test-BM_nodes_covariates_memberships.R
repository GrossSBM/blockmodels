set.seed(20260226)

test_that("dispatcher SBM exports covariate membership matrices", {
  n <- 12
  Q <- 3
  npc <- n / Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  M <- matrix(rbinom(n * n, 1, 0.4), n, n)
  X <- cbind(1, runif(n), runif(n))

  res_sbm <- blockmodels:::dispatcher(
    membership_name = "SBM",
    membership_init = list(
      Z = Z,
      nodes_covariates = X
    ),
    model_name = "bernoulli",
    network = list(adjacency = M),
    real_EM = TRUE
  )

  expect_true("B" %in% names(res_sbm$membership))
  expect_true(is.matrix(res_sbm$membership$B))
  expect_equal(dim(res_sbm$membership$B), c(ncol(X), ncol(Z)))
  expect_true(is.matrix(res_sbm$membership$alpha))
  expect_equal(dim(res_sbm$membership$alpha), c(nrow(Z), ncol(Z)))
})

test_that("dispatcher LBM exports covariate membership matrices", {
  n1 <- 10
  n2 <- 8
  Q1 <- 2
  Q2 <- 2
  npc1 <- n1 / Q1
  npc2 <- n2 / Q2

  Z1 <- diag(Q1) %x% matrix(1, npc1, 1)
  Z2 <- diag(Q2) %x% matrix(1, npc2, 1)
  M_lbm <- matrix(rbinom(n1 * n2, 1, 0.35), n1, n2)
  X_row <- cbind(1, runif(n1), runif(n1))
  X_col <- cbind(1, runif(n2))

  res_lbm <- blockmodels:::dispatcher(
    membership_name = "LBM",
    membership_init = list(
      Z1 = Z1,
      Z2 = Z2,
      row_covariates = X_row,
      col_covariates = X_col
    ),
    model_name = "bernoulli",
    network = list(adjacency = M_lbm),
    real_EM = TRUE
  )

  expect_true("B" %in% names(res_lbm$membership))
  expect_true("G" %in% names(res_lbm$membership))
  expect_true(is.matrix(res_lbm$membership$B))
  expect_true(is.matrix(res_lbm$membership$G))
  expect_equal(dim(res_lbm$membership$B), c(ncol(X_row), ncol(Z1)))
  expect_equal(dim(res_lbm$membership$G), c(ncol(X_col), ncol(Z2)))
  expect_true(is.matrix(res_lbm$membership$alpha1))
  expect_true(is.matrix(res_lbm$membership$alpha2))
  expect_equal(dim(res_lbm$membership$alpha1), c(nrow(Z1), ncol(Z1)))
  expect_equal(dim(res_lbm$membership$alpha2), c(nrow(Z2), ncol(Z2)))
})
