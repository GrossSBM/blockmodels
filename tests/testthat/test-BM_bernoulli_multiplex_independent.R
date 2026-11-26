set.seed(12)

simulate_layers <- function(Z1, Z2, P) {
  lapply(P, function(Pk) 1 * (matrix(runif(nrow(Z1) * nrow(Z2)), nrow(Z1), nrow(Z2)) < Z1 %*% Pk %*% t(Z2)))
}
sym <- function(M) {M[lower.tri(M)] <- t(M)[lower.tri(M)]; M}

test_that("BM_bernoulli_multiplex_independent SBM, SBM_sym and LBM estimations run", {
  Z <- diag(2) %x% matrix(1, 10, 1)
  layers <- simulate_layers(Z, Z, list(matrix(runif(4), 2), matrix(runif(4), 2)))
  model <- BM_bernoulli_multiplex_independent("SBM", layers, plotting = "", explore_min = 2, explore_max = 2,
                                              ncores = 1, verbosity = 0)
  expect_model_estimation(model)
  expect_length(model$model_parameters[[2]]$pi, 2)
  expect_equal(dim(model$model_parameters[[2]]$pi[[1]]), c(2, 2))
  expect_equal(model$model_parameters[[2]]$n_parameters, 2 * 2 * 2)

  model <- BM_bernoulli_multiplex_independent("SBM_sym", lapply(layers, sym), plotting = "", explore_min = 2,
                                              explore_max = 2, ncores = 1, verbosity = 0)
  expect_model_estimation(model)
  expect_equal(model$model_parameters[[2]]$n_parameters, 2 * 3)

  Z1 <- matrix(1, 20, 1)
  Z2 <- diag(2) %x% matrix(1, 10, 1)
  layers <- simulate_layers(Z1, Z2, list(matrix(runif(2), 1), matrix(runif(2), 1)))
  model <- BM_bernoulli_multiplex_independent("LBM", layers, plotting = "", explore_min = 2, explore_max = 2,
                                              ncores = 1, verbosity = 0)
  expect_model_estimation(model)
  expect_equal(length(model$prediction(2)), 2)
})

test_that("BM_bernoulli_multiplex_independent with one layer is BM_bernoulli", {
  Z <- diag(3) %x% matrix(1, 15, 1)
  M <- simulate_layers(Z, Z, list(matrix(c(.8, .1, .1, .1, .6, .2, .1, .2, .5), 3)))[[1]]
  for (membership in c("SBM", "SBM_sym")) {
    X <- if (membership == "SBM_sym") sym(M) else M
    a <- BM_bernoulli(membership, X, plotting = "", explore_max = 3, ncores = 1, verbosity = 0)
    b <- BM_bernoulli_multiplex_independent(membership, list(X), plotting = "", explore_max = 3, ncores = 1,
                                            verbosity = 0)
    a$estimate()
    b$estimate()
    # same model, up to the convergence of the EM
    expect_equal(b$ICL, a$ICL, tolerance = 1e-5)
    # same estimates, up to the labels of the blocks
    expect_equal(sort(b$model_parameters[[3]]$pi[[1]]), sort(a$model_parameters[[3]]$pi), tolerance = 1e-5)
  }
})

test_that("BM_bernoulli_multiplex_independent recovers the blocks and the parameters", {
  Z <- diag(3) %x% matrix(1, 30, 1)
  P <- list(matrix(c(.8, .1, .1, .1, .6, .2, .1, .2, .5), 3), matrix(c(.1, .7, .7, .7, .1, .7, .7, .7, .1), 3))
  layers <- simulate_layers(Z, Z, P)
  model <- BM_bernoulli_multiplex_independent("SBM", layers, plotting = "", explore_max = 5, ncores = 1,
                                              verbosity = 0)
  model$estimate()
  expect_equal(which.max(model$ICL), 3)
  cl <- apply(model$memberships[[3]]$Z, 1, which.max)
  truth <- rep(1:3, each = 30)
  expect_equal(max(table(cl, truth)), 30)
  perm <- apply(table(cl, truth), 2, which.max)
  for (k in 1:2) {
    expect_lt(max(abs(model$model_parameters[[3]]$pi[[k]][perm, perm] - P[[k]])), 0.1)
  }
})

test_that("BM_bernoulli_multiplex_independent handles empty and full layers", {
  Z <- diag(2) %x% matrix(1, 15, 1)
  M <- simulate_layers(Z, Z, list(matrix(c(.8, .2, .2, .6), 2)))[[1]]
  for (membership in c("SBM", "SBM_sym")) {
    X <- if (membership == "SBM_sym") sym(M) else M
    model <- BM_bernoulli_multiplex_independent(membership, list(X, 0 * X, 1 + 0 * X), plotting = "",
                                                explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
    expect_model_estimation(model)
    expect_true(all(is.finite(model$ICL)))
    expect_true(all(model$model_parameters[[2]]$pi[[2]] == 0))
    expect_true(all(model$model_parameters[[2]]$pi[[3]] == 1))
  }
  M <- 1 * (matrix(runif(20 * 30), 20, 30) < 0.4)
  model <- BM_bernoulli_multiplex_independent("LBM", list(M, 0 * M), plotting = "", explore_min = 2,
                                              explore_max = 2, ncores = 1, verbosity = 0)
  expect_model_estimation(model)
  expect_true(all(is.finite(model$ICL[-1])))
})

add_na <- function(M, prop, symmetric = FALSE) {
  idx <- which(row(M) != col(M) & (!symmetric | row(M) < col(M)))
  idx <- sample(idx, floor(prop * length(idx)))
  M[idx] <- NA
  if (symmetric) M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

test_that("BM_bernoulli_multiplex_independent with one layer and NA is BM_bernoulli with NA", {
  Z <- diag(3) %x% matrix(1, 15, 1)
  M <- simulate_layers(Z, Z, list(matrix(c(.8, .1, .1, .1, .6, .2, .1, .2, .5), 3)))[[1]]
  for (membership in c("SBM", "SBM_sym")) {
    X <- add_na(if (membership == "SBM_sym") sym(M) else M, 0.2, symmetric = membership == "SBM_sym")
    a <- BM_bernoulli(membership, X, plotting = "", explore_max = 3, ncores = 1, verbosity = 0)
    b <- BM_bernoulli_multiplex_independent(membership, list(X), plotting = "", explore_max = 3, ncores = 1,
                                            verbosity = 0)
    expect_equal(b$data_number(), a$data_number())
    a$estimate()
    b$estimate()
    expect_equal(b$ICL, a$ICL, tolerance = 1e-5)
    # same estimates, up to the labels of the blocks and the convergence of the EM
    expect_equal(sort(b$model_parameters[[3]]$pi[[1]]), sort(a$model_parameters[[3]]$pi), tolerance = 1e-3)
  }

  # LBM: the explorations of the two models may differ (different initializations),
  # so compare one estimation from the same initial membership
  M <- 1 * (matrix(runif(20 * 30), 20, 30) < 0.4)
  X <- add_na(M, 0.2)
  a <- BM_bernoulli("LBM", X, plotting = "", ncores = 1, verbosity = 0)
  b <- BM_bernoulli_multiplex_independent("LBM", list(X), plotting = "", ncores = 1, verbosity = 0)
  expect_equal(b$data_number(), a$data_number())
  init <- list(Z1 = t(rmultinom(20, 1, rep(1, 2))), Z2 = t(rmultinom(30, 1, rep(1, 3))))
  ra <- blockmodels:::dispatcher("LBM", init, "bernoulli", list(adjacency = X), FALSE)
  rb <- blockmodels:::dispatcher("LBM", init, "bernoulli_multiplex_independent", list(adjacency = list(X)), FALSE)
  expect_equal(rb$PL, ra$PL, tolerance = 1e-8)
  expect_equal(rb$model$pi[[1]], ra$model$pi, tolerance = 1e-8)
  expect_equal(rb$membership$Z1, ra$membership$Z1, tolerance = 1e-8)
})

test_that("BM_bernoulli_multiplex_independent handles NA on some layers only", {
  Z <- diag(3) %x% matrix(1, 30, 1)
  P <- list(matrix(c(.8, .1, .1, .1, .6, .2, .1, .2, .5), 3), matrix(c(.1, .7, .7, .7, .1, .7, .7, .7, .1), 3))
  layers <- simulate_layers(Z, Z, P)
  # 30% of missing values on the first layer, the second one is complete
  layers[[1]] <- add_na(layers[[1]], 0.3)
  model <- BM_bernoulli_multiplex_independent("SBM", layers, plotting = "", explore_max = 5, ncores = 1,
                                              verbosity = 0)
  model$estimate()
  expect_true(all(is.finite(model$ICL)))
  expect_equal(which.max(model$ICL), 3)
  cl <- apply(model$memberships[[3]]$Z, 1, which.max)
  truth <- rep(1:3, each = 30)
  expect_equal(max(table(cl, truth)), 30)
  perm <- apply(table(cl, truth), 2, which.max)
  for (k in 1:2) {
    expect_lt(max(abs(model$model_parameters[[3]]$pi[[k]][perm, perm] - P[[k]])), 0.1)
  }

  # NA on all the layers of some dyads, in SBM_sym and LBM
  layers <- lapply(simulate_layers(Z, Z, P), sym)
  idx <- is.na(add_na(layers[[1]], 0.2, symmetric = TRUE))
  layers <- lapply(layers, function(M) {M[idx] <- NA; M})
  model <- BM_bernoulli_multiplex_independent("SBM_sym", layers, plotting = "", explore_min = 3,
                                              explore_max = 3, ncores = 1, verbosity = 0)
  expect_equal(model$data_number(), 90 * 89 / 2 - sum(idx) / 2)
  expect_model_estimation(model)
  expect_true(all(is.finite(model$ICL)))

  layers <- simulate_layers(Z, diag(2) %x% matrix(1, 20, 1), list(matrix(runif(6), 3), matrix(runif(6), 3)))
  layers[[2]] <- add_na(layers[[2]], 0.3)
  model <- BM_bernoulli_multiplex_independent("LBM", layers, plotting = "", explore_min = 2, explore_max = 3,
                                              ncores = 1, verbosity = 0)
  expect_model_estimation(model)
  expect_true(all(is.finite(model$ICL[-1])))
})

test_that("BM_bernoulli_multiplex_independent fails gracefully when a layer is all NA", {
  M <- 1 * (matrix(runif(400), 20, 20) < 0.3)
  model <- BM_bernoulli_multiplex_independent("SBM", list(M, NA * M), plotting = "", explore_min = 2,
                                              explore_max = 2, ncores = 1, verbosity = 0)
  expect_error(model$estimate(), "No valid edges in adjacency matrix. They are all NAs\\.")
})

test_that("BM_bernoulli_multiplex_independent ICL penalty sums log(observed dyads) over the layers", {
  Z <- diag(2) %x% matrix(1, 15, 1)
  P <- list(matrix(c(.8, .2, .2, .6), 2), matrix(c(.3, .6, .6, .3), 2), matrix(c(.5, .1, .1, .5), 2))
  layers <- lapply(simulate_layers(Z, Z, P), sym)
  # different proportions of missing values on the layers
  layers <- mapply(function(M, prop) add_na(M, prop, symmetric = TRUE), layers, c(0, 0.2, 0.5), SIMPLIFY = FALSE)
  n <- 30
  n_obs <- list(
    SBM = vapply(layers, function(M) sum(!is.na(M)) - n, numeric(1)),
    SBM_sym = vapply(layers, function(M) (sum(!is.na(M)) - n) / 2, numeric(1))
  )
  expect_equal(n_obs$SBM[1], n * (n - 1))
  for (membership in c("SBM", "SBM_sym")) {
    m <- BM_bernoulli_multiplex_independent(membership, layers, plotting = "", explore_min = 3, explore_max = 3,
                                            ncores = 1, verbosity = 0)
    m$estimate()
    for (Q in seq_along(m$ICL)) {
      n_param_layer <- if (membership == "SBM") Q^2 else Q * (Q + 1) / 2
      expect_equal(m$model_parameters[[Q]]$n_parameters, 3 * n_param_layer)
      expected <- m$PL[Q] - .5 * (sum(n_param_layer * log(n_obs[[membership]])) + m$memberships[[Q]]$ICL_penalty())
      expect_equal(m$ICL[Q], expected, tolerance = 1e-8)
    }
  }

  layers <- simulate_layers(Z, diag(2) %x% matrix(1, 10, 1), list(matrix(runif(4), 2), matrix(runif(4), 2)))
  layers[[2]] <- add_na(layers[[2]], 0.4)
  m <- BM_bernoulli_multiplex_independent("LBM", layers, plotting = "", explore_min = 4, explore_max = 4,
                                          ncores = 1, verbosity = 0)
  m$estimate()
  for (Q in which(is.finite(m$ICL))) {
    Q1 <- ncol(m$memberships[[Q]]$Z1); Q2 <- ncol(m$memberships[[Q]]$Z2)
    expected <- m$PL[Q] - .5 * (Q1 * Q2 * (log(30 * 20) + log(sum(!is.na(layers[[2]])))) +
                                  m$memberships[[Q]]$ICL_penalty())
    expect_equal(m$ICL[Q], expected, tolerance = 1e-8)
  }
})
