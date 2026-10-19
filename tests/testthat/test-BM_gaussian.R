set.seed(12)

test_that("BM_gaussian SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu <- 20 * matrix(runif(Q * Q), Q, Q)
  M <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu %*% t(Z)

  model <- BM_gaussian("SBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu <- 20 * matrix(runif(Q * Q), Q, Q)
  Mu[lower.tri(Mu)] <- t(Mu)[lower.tri(Mu)]
  M <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu %*% t(Z)
  M[lower.tri(M)] <- t(M)[lower.tri(M)]

  model <- BM_gaussian("SBM_sym", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  Mu <- 20 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  M <- matrix(rnorm(n[1] * n[2], sd = 10), n[1], n[2]) + Z1 %*% Mu %*% t(Z2)

  model <- BM_gaussian("LBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian SBM with nodes covariates runs and passes nodes covariates", {
    npc <- 10
    Q <- 2
    n <- npc * Q
    Z <- diag(Q) %x% matrix(1, npc, 1)
    X <- cbind(1, runif(n), runif(n))
    Mu <- 20 * matrix(runif(Q * Q), Q, Q)
    M <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu %*% t(Z)

    model <- BM_gaussian(
        "SBM",
        M,
        nodes_covariates = list(node = X),
        plotting = "",
        explore_min = 2,
        explore_max = 4,
        ncores = 1,
        verbosity = 0
    )

    expect_model_estimation(model)
    expect_model_nodes_covariates(model,X)
})

test_that("BM_gaussian SBM_sym with nodes covariates runs and passes nodes covariates", {
    npc <- 10
    Q <- 2
    n <- npc * Q
    Z <- diag(Q) %x% matrix(1, npc, 1)
    X <- cbind(1, runif(n), runif(n))
    Mu <- 20 * matrix(runif(Q * Q), Q, Q)
    M <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu %*% t(Z)
    M[lower.tri(M)] <- t(M)[lower.tri(M)]

    model <- BM_gaussian("SBM_sym", M, nodes_covariates = list(node = X), plotting = "", explore_min = 2, explore_max = 4, ncores = 2, verbosity = 0)
    expect_model_estimation(model)
    expect_model_nodes_covariates(model, X)
})

test_that("BM_bernoulli LBM with nodes covariates runs and passes nodes covariates", {
    npc <- c(20, 10)
    Q <- c(1, 2)
    n <- npc * Q
    Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
    Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
    X_row <- cbind(1, runif(n[1]), runif(n[1]))
    X_col <- cbind(1, runif(n[2]), runif(n[2]))
    Mu <- 20 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
    M <- matrix(rnorm(n[1] * n[2], sd = 10), n[1], n[2]) + Z1 %*% Mu %*% t(Z2)


    model <- BM_gaussian("LBM", M, nodes_covariates = list(row = X_row, col = X_col), plotting = "", explore_min = 2, explore_max = 4, ncores = 2, verbosity = 0)
    expect_model_estimation(model)
    expect_model_row_col_covariates(model = model, expected_row_covariates = X_row, expected_col_covariates = X_col)
})

test_that("BM_gaussian SBM handles partial NA values", {
    set.seed(12)

    npc <- 10
    Q <- 2
    n <- npc * Q
    Z <- diag(Q) %x% matrix(1, npc, 1)
    Mu <- 20 * matrix(runif(Q * Q), Q, Q)
    M <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu %*% t(Z)

    off_diag <- which(row(M) != col(M))
    set.seed(123)
    M[sample(off_diag, max(1, floor(0.1 * length(off_diag))))] <- NA_real_

    model <- BM_gaussian(
        "SBM",
        M,
        plotting = "",
        explore_min = 2,
        explore_max = 2,
        ncores = 1,
        verbosity = 0
    )

    expect_no_error(model$estimate())
    expect_true(is.numeric(model$ICL))
    expect_true(all(is.finite(model$ICL)))
})

test_that("BM_gaussian SBM fails gracefully when all values are NA", {
    n <- 20
    M_all_na <- matrix(NA_real_, n, n)

    model <- BM_gaussian(
        "SBM",
        M_all_na,
        plotting = "",
        explore_min = 2,
        explore_max = 2,
        ncores = 1,
        verbosity = 0
    )

    expect_error(
        model$estimate(),
        "No valid edges in adjacency matrix. They are all NAs\\."
    )
})

test_that("BM_gaussian SBM_sym handles partial NA values", {
    set.seed(12)

    npc <- 10
    Q <- 2
    n <- npc * Q
    Z <- diag(Q) %x% matrix(1, npc, 1)
    Mu <- 20 * matrix(runif(Q * Q), Q, Q)
    Mu[lower.tri(Mu)] <- t(Mu)[lower.tri(Mu)]
    M <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu %*% t(Z)
    M[lower.tri(M)] <- t(M)[lower.tri(M)]

    upper_off_diag <- which(upper.tri(M))
    set.seed(123)
    idx <- sample(upper_off_diag, max(1, floor(0.1 * length(upper_off_diag))))
    idx_rc <- arrayInd(idx, dim(M))
    idx_sym <- cbind(idx_rc[, 2], idx_rc[, 1])
    M[idx] <- NA_real_
    M[idx_sym] <- NA_real_

    model <- BM_gaussian(
        "SBM_sym",
        M,
        plotting = "",
        explore_min = 2,
        explore_max = 2,
        ncores = 1,
        verbosity = 0
    )

    expect_no_error(model$estimate())
    expect_true(is.numeric(model$ICL))
    expect_true(all(is.finite(model$ICL)))
})

test_that("BM_gaussian SBM_sym fails gracefully when all values are NA", {
    n <- 20
    M_all_na <- matrix(NA_real_, n, n)

    model <- BM_gaussian(
        "SBM_sym",
        M_all_na,
        plotting = "",
        explore_min = 2,
        explore_max = 2,
        ncores = 1,
        verbosity = 0
    )

    expect_error(
        model$estimate(),
        "No valid edges in adjacency matrix. They are all NAs\\."
    )
})

test_that("BM_gaussian LBM handles partial NA values", {
    set.seed(12)

    npc <- c(20, 10)
    Q <- c(1, 2)
    n <- npc * Q
    Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
    Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
    Mu <- 20 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
    M <- matrix(rnorm(n[1] * n[2], sd = 10), n[1], n[2]) + Z1 %*% Mu %*% t(Z2)

    set.seed(123)
    M[sample(seq_along(M), max(1, floor(0.1 * length(M))))] <- NA_real_

    model <- BM_gaussian(
        "LBM",
        M,
        plotting = "",
        explore_min = 2,
        explore_max = 2,
        ncores = 1,
        verbosity = 0
    )

    expect_no_error(model$estimate())
    expect_true(is.numeric(model$ICL[-1]))
    expect_true(all(is.finite(model$ICL[-1])))
})

test_that("BM_gaussian LBM fails gracefully when all values are NA", {
    n1 <- 20
    n2 <- 20
    M_all_na <- matrix(NA_real_, n1, n2)

    model <- BM_gaussian(
        "LBM",
        M_all_na,
        plotting = "",
        explore_min = 2,
        explore_max = 2,
        ncores = 1,
        verbosity = 0
    )

    expect_error(
        model$estimate(),
        "No valid edges in adjacency matrix. They are all NAs\\."
    )
})

