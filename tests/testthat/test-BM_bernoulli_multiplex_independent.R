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
