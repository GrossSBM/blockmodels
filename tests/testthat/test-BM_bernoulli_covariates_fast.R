set.seed(12)

test_that("BM_bernoulli_covariates_fast SBM estimation runs", {
    npc <- 10
    Q <- 2
    n <- npc * Q
    sigmo <- function(x) {
        1 / (1 + exp(-x))
    }
    Z <- diag(Q) %x% matrix(1, npc, 1)
    Mg <- 8 * matrix(runif(Q * Q), Q, Q) - 4
    Y1 <- matrix(runif(n * n), n, n) - 0.5
    Y2 <- matrix(runif(n * n), n, n) - 0.5
    M_in_expectation <- sigmo(Z %*% Mg %*% t(Z) + 5 * Y1 - 3 * Y2)
    M <- 1 * (matrix(runif(n * n), n, n) < M_in_expectation)

    model <- BM_bernoulli_covariates_fast("SBM", M, list(Y1, Y2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
    expect_model_estimation(model)
})

test_that("BM_bernoulli_covariates_fast SBM_sym estimation runs", {
    npc <- 10
    Q <- 2
    n <- npc * Q
    sigmo <- function(x) {
        1 / (1 + exp(-x))
    }
    Z <- diag(Q) %x% matrix(1, npc, 1)
    Mg <- 8 * matrix(runif(Q * Q), Q, Q) - 4
    Mg[lower.tri(Mg)] <- t(Mg)[lower.tri(Mg)]
    Y1 <- matrix(runif(n * n), n, n) - 0.5
    Y2 <- matrix(runif(n * n), n, n) - 0.5
    Y1[lower.tri(Y1)] <- t(Y1)[lower.tri(Y1)]
    Y2[lower.tri(Y2)] <- t(Y2)[lower.tri(Y2)]
    M_in_expectation <- sigmo(Z %*% Mg %*% t(Z) + 5 * Y1 - 3 * Y2)
    M <- 1 * (matrix(runif(n * n), n, n) < M_in_expectation)
    M[lower.tri(M)] <- t(M)[lower.tri(M)]

    model <- BM_bernoulli_covariates_fast("SBM_sym", M, list(Y1, Y2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
    expect_model_estimation(model)
})

test_that("BM_bernoulli_covariates_fast LBM estimation runs", {
    npc <- c(20, 10)
    Q <- c(1, 2)
    n <- npc * Q
    sigmo <- function(x) {
        1 / (1 + exp(-x))
    }
    Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
    Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
    Mg <- 8 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2]) - 4
    Y1 <- matrix(runif(n[1] * n[2]), n[1], n[2]) - 0.5
    Y2 <- matrix(runif(n[1] * n[2]), n[1], n[2]) - 0.5
    M_in_expectation <- sigmo(Z1 %*% Mg %*% t(Z2) + 5 * Y1 - 3 * Y2)
    M <- 1 * (matrix(runif(n[1] * n[2]), n[1], n[2]) < M_in_expectation)

    model <- BM_bernoulli_covariates_fast("LBM", M, list(Y1, Y2), plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
    expect_model_estimation(model)
})


test_that("BM_bernoulli_covariates_fast SBM handles partial NA and is invariant to na_replace_value", {
    set.seed(111)

    npc <- 10
    Q <- 2
    n <- npc * Q
    sigmo <- function(x) 1 / (1 + exp(-x))
    Z <- diag(Q) %x% matrix(1, npc, 1)
    Mg <- 8 * matrix(runif(Q * Q), Q, Q) - 4
    Y1 <- matrix(runif(n * n), n, n) - 0.5
    Y2 <- matrix(runif(n * n), n, n) - 0.5
    P <- sigmo(Z %*% Mg %*% t(Z) + 4 * Y1 - 1.5 * Y2)
    M <- 1 * (matrix(runif(n * n), n, n) < P)

    off_diag <- which(row(M) != col(M))
    set.seed(112)
    M[sample(off_diag, max(1, floor(0.15 * length(off_diag))))] <- NA_real_

    set.seed(113)
    model_default <- BM_bernoulli_covariates_fast("SBM", M, list(Y1, Y2), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
    expect_no_error(model_default$estimate())

    set.seed(113)
    model_custom <- BM_bernoulli_covariates_fast("SBM", M, list(Y1, Y2), na_replace_value = -321, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
    expect_no_error(model_custom$estimate())

    expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
    k <- which.max(model_default$ICL)
    expect_equal(model_default$model_parameters[[k]]$m, model_custom$model_parameters[[k]]$m, tolerance = 1e-10)
    expect_equal(model_default$model_parameters[[k]]$beta, model_custom$model_parameters[[k]]$beta, tolerance = 1e-10)
})
