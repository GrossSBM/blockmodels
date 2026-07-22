set.seed(12)

test_that("BM_poisson SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  L <- 70 * matrix(runif(Q * Q), Q, Q)
  M_in_expectation <- Z %*% L %*% t(Z)
  M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)

  model <- BM_poisson("SBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_poisson SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  L <- 70 * matrix(runif(Q * Q), Q, Q)
  L[lower.tri(L)] <- t(L)[lower.tri(L)]
  M_in_expectation <- Z %*% L %*% t(Z)
  M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)
  M[lower.tri(M)] <- t(M)[lower.tri(M)]

  model <- BM_poisson("SBM_sym", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_poisson LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  L <- 70 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  M_in_expectation <- Z1 %*% L %*% t(Z2)
  M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n[1], n[2])

  model <- BM_poisson("LBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_poisson SBM with nodes covariates runs and passes nodes covariates", {
    npc <- 10
    Q <- 2
    n <- npc * Q
    Z <- diag(Q) %x% matrix(1, npc, 1)
    X <- cbind(1, runif(n), runif(n))
    L <- 70 * matrix(runif(Q * Q), Q, Q)
    M_in_expectation <- Z %*% L %*% t(Z)
    M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)

    model <- BM_poisson(
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
    expect_model_nodes_covariates(model, X)
})

test_that("BM_poisson SBM_sym with nodes covariates runs and passes nodes covariates", {
    npc <- 10
    Q <- 2
    n <- npc * Q
    Z <- diag(Q) %x% matrix(1, npc, 1)
    X <- cbind(1, runif(n), runif(n))
    L <- 70 * matrix(runif(Q * Q), Q, Q)
    M_in_expectation <- Z %*% L %*% t(Z)
    M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)
    M[lower.tri(M)] <- t(M)[lower.tri(M)]

    model <- BM_poisson(
        "SBM_sym",
        M,
        nodes_covariates = list(node = X),
        plotting = "",
        explore_min = 2,
        explore_max = 4,
        ncores = 1,
        verbosity = 0
    )

    expect_model_estimation(model)
    expect_model_nodes_covariates(model, X)
})

test_that("BM_poisson LBM with nodes covariates runs and passes nodes covariates", {
    npc <- c(20, 10)
    Q <- c(1, 2)
    n <- npc * Q
    Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
    Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
    X_row <- cbind(1, runif(n[1]), runif(n[1]))
    X_col <- cbind(1, runif(n[2]), runif(n[2]))
    L <- 70 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
    M_in_expectation <- Z1 %*% L %*% t(Z2)
    M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n[1], n[2])

    model <- BM_poisson("LBM", M,
        nodes_covariates = list(row = X_row, col = X_col), 
        plotting = "", 
        explore_min = 2, 
        explore_max = 4, 
        ncores = 2, verbosity = 0
    )
    expect_model_estimation(model)
})

test_that("BM_poisson SBM handles partial NA values", {
    set.seed(12)

    npc <- 10
    Q <- 2
    n <- npc * Q
    Z <- diag(Q) %x% matrix(1, npc, 1)
    L <- 70 * matrix(runif(Q * Q), Q, Q)
    M_in_expectation <- Z %*% L %*% t(Z)
    M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)

    off_diag <- which(row(M) != col(M))
    set.seed(123)
    M[sample(off_diag, max(1, floor(0.1 * length(off_diag))))] <- NA_real_

    model <- BM_poisson(
        "SBM",
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

test_that("BM_poisson SBM fails gracefully when all values are NA", {
    n <- 20
    M_all_na <- matrix(NA_real_, n, n)

    model <- BM_poisson(
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

test_that("BM_poisson SBM_sym handles partial NA values", {
    set.seed(12)

    npc <- 10
    Q <- 2
    n <- npc * Q
    Z <- diag(Q) %x% matrix(1, npc, 1)
    L <- 70 * matrix(runif(Q * Q), Q, Q)
    L[lower.tri(L)] <- t(L)[lower.tri(L)]
    M_in_expectation <- Z %*% L %*% t(Z)
    M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)
    M[lower.tri(M)] <- t(M)[lower.tri(M)]

    upper_off_diag <- which(upper.tri(M))
    set.seed(123)
    idx <- sample(upper_off_diag, max(1, floor(0.1 * length(upper_off_diag))))
    idx_rc <- arrayInd(idx, dim(M))
    idx_sym <- cbind(idx_rc[, 2], idx_rc[, 1])
    M[idx] <- NA_real_
    M[idx_sym] <- NA_real_

    model <- BM_poisson(
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

test_that("BM_poisson SBM_sym fails gracefully when all values are NA", {
    n <- 20
    M_all_na <- matrix(NA_real_, n, n)

    model <- BM_poisson(
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

test_that("BM_poisson LBM handles partial NA values", {
    set.seed(12)

    npc <- c(20, 10)
    Q <- c(1, 2)
    n <- npc * Q
    Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
    Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
    L <- 70 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
    M_in_expectation <- Z1 %*% L %*% t(Z2)
    M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n[1], n[2])

    set.seed(123)
    M[sample(seq_along(M), max(1, floor(0.1 * length(M))))] <- NA_real_

    model <- BM_poisson(
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

test_that("BM_poisson LBM fails gracefully when all values are NA", {
    n1 <- 20
    n2 <- 20
    M_all_na <- matrix(NA_real_, n1, n2)

    model <- BM_poisson(
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

test_that("BM_poisson computation is invariant to na_replace_value", {
    set.seed(321)

    npc <- 12
    Q <- 2
    n <- npc * Q
    Z <- diag(Q) %x% matrix(1, npc, 1)
    L <- 70 * matrix(runif(Q * Q), Q, Q)
    M_in_expectation <- Z %*% L %*% t(Z)
    M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)

    off_diag <- which(row(M) != col(M))
    set.seed(654)
    M[sample(off_diag, max(1, floor(0.15 * length(off_diag))))] <- NA_real_

    set.seed(999)
    model_default <- BM_poisson(
        "SBM",
        M,
        na_replace_value = 0,
        plotting = "",
        explore_min = 2,
        explore_max = 2,
        ncores = 1,
        verbosity = 0
    )
    model_default$estimate()

    set.seed(999)
    model_custom <- BM_poisson(
        "SBM",
        M,
        na_replace_value = 42,
        plotting = "",
        explore_min = 2,
        explore_max = 2,
        ncores = 1,
        verbosity = 0
    )
    model_custom$estimate()

    expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
    max_ICL <- which.max(model_default$ICL)
    expect_equal(model_default$model_parameters[[max_ICL]]$lambda, model_custom$model_parameters[[max_ICL]]$lambda, tolerance = 1e-10)
})

test_that("BM_poisson validates na_replace_value input", {
    n <- 20
    M <- matrix(0, n, n)

    expect_error(
        BM_poisson(
            "SBM",
            M,
            na_replace_value = NA_real_,
            plotting = "",
            explore_min = 2,
            explore_max = 2,
            ncores = 1,
            verbosity = 0
        ),
        "na_replace_value must be a single finite numeric value\\."
    )

    expect_error(
        BM_poisson(
            "SBM",
            M,
            na_replace_value = c(0, 1),
            plotting = "",
            explore_min = 2,
            explore_max = 2,
            ncores = 1,
            verbosity = 0
        ),
        "na_replace_value must be a single finite numeric value\\."
    )
})