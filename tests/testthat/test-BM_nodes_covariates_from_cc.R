set.seed(20260303)

test_that("SBM preserves nodes covariates through from_cc and to_cc", {
  n <- 12
  Q <- 3
  npc <- n / Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  X <- cbind(1, runif(n), runif(n))

  membership <- methods::getRefClass("SBM")(
    from_cc = list(Z = Z, nodes_covariates = X)
  )

  expect_true(is.matrix(membership$nodes_covariates))
  expect_equal(membership$nodes_covariates, X)
  expect_true("nodes_covariates" %in% names(membership$to_cc()))
  expect_equal(membership$to_cc()$nodes_covariates, X)
})

test_that("SBM_sym preserves nodes covariates through from_cc and to_cc", {
  n <- 12
  Q <- 3
  npc <- n / Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  X <- cbind(1, runif(n), runif(n))

  membership <- methods::getRefClass("SBM_sym")(
    from_cc = list(Z = Z, nodes_covariates = X)
  )

  expect_true(is.matrix(membership$nodes_covariates))
  expect_equal(membership$nodes_covariates, X)
  expect_true("nodes_covariates" %in% names(membership$to_cc()))
  expect_equal(membership$to_cc()$nodes_covariates, X)
})

test_that("LBM preserves row and col covariates through from_cc and to_cc", {
  n1 <- 10
  n2 <- 8
  Q1 <- 2
  Q2 <- 2
  npc1 <- n1 / Q1
  npc2 <- n2 / Q2
  Z1 <- diag(Q1) %x% matrix(1, npc1, 1)
  Z2 <- diag(Q2) %x% matrix(1, npc2, 1)
  X_row <- cbind(1, runif(n1), runif(n1))
  X_col <- cbind(1, runif(n2))

  membership <- methods::getRefClass("LBM")(
    from_cc = list(Z1 = Z1, Z2 = Z2, row_covariates = X_row, col_covariates = X_col)
  )

  expect_true(is.matrix(membership$row_covariates))
  expect_true(is.matrix(membership$col_covariates))
  expect_equal(membership$row_covariates, X_row)
  expect_equal(membership$col_covariates, X_col)
  expect_true("row_covariates" %in% names(membership$to_cc()))
  expect_true("col_covariates" %in% names(membership$to_cc()))
  expect_equal(membership$to_cc()$row_covariates, X_row)
  expect_equal(membership$to_cc()$col_covariates, X_col)
})