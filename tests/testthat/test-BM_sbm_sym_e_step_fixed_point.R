## Regression tests for issue #8: for undirected memberships (SBM_sym), the
## posteriors returned by the E step must satisfy the mean-field fixed point
## that counts each dyad *once*, not the one that counts it twice (which is
## what SBM::e_step() computes for a directed membership).
##
## The expected fixed point is computed independently in R from the model
## parameters returned by the (internal) dispatcher, following the same
## methodology as the reproduction in
## https://github.com/GrossSBM/blockmodels/issues/8

simulate_sym_bernoulli <- function(n, Q, within = .45, between = .3) {
  cl <- sample(Q, n, TRUE)
  P <- matrix(between, Q, Q)
  diag(P) <- within
  A <- 1 * (matrix(runif(n * n), n) < P[cl, cl])
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  list(A = A, cl = cl)
}

initial_Z <- function(cl, Q) {
  Zi <- diag(Q)[cl, ] * 0.6 + 0.2
  Zi / rowSums(Zi)
}

fixed_point_from_data_term <- function(data_term, alpha) {
  l <- sweep(data_term, 2, log(alpha), "+")
  l <- exp(l - apply(l, 1, max))
  l / rowSums(l)
}

test_that("BM_bernoulli SBM_sym E step matches the single-count fixed point", {
  set.seed(3)
  n <- 40
  Q <- 2
  sim <- simulate_sym_bernoulli(n, Q)
  A <- sim$A
  Zi <- initial_Z(sim$cl, Q)

  r <- blockmodels:::dispatcher("SBM_sym", list(Z = Zi), "bernoulli", list(adjacency = A), TRUE)
  Z <- r$membership$Z
  alpha <- r$membership$alpha
  pi <- r$model$pi

  # Each dyad counted once (the correct, single-count data term).
  data_term_once <- A %*% Z %*% log(pi) + ((1 - A) * (1 - diag(n))) %*% Z %*% log(1 - pi)
  # Each dyad counted twice (the pre-fix, buggy behaviour).
  data_term_twice <- 2 * data_term_once

  residual_once <- max(abs(Z - fixed_point_from_data_term(data_term_once, alpha)))
  residual_twice <- max(abs(Z - fixed_point_from_data_term(data_term_twice, alpha)))

  expect_lt(residual_once, 0.05)
  # The buggy ("twice") fixed point must now be clearly off, otherwise the
  # test would not be discriminating.
  expect_gt(residual_twice, 0.2)
})

test_that("BM_poisson SBM_sym E step matches the single-count fixed point", {
  set.seed(3)
  n <- 30
  Q <- 2
  cl <- sample(Q, n, TRUE)
  Lambda <- matrix(c(2.5, 3, 3, 2.5), Q, Q)
  A <- matrix(rpois(n * n, Lambda[cl, cl]), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  Zi <- initial_Z(cl, Q)

  r <- blockmodels:::dispatcher("SBM_sym", list(Z = Zi), "poisson", list(adjacency = A), TRUE)
  Z <- r$membership$Z
  alpha <- r$membership$alpha
  lambda <- r$model$lambda

  off_diag <- 1 - diag(n)
  data_term_once <- A %*% Z %*% log(lambda) - (off_diag %*% Z) %*% lambda
  data_term_twice <- 2 * data_term_once

  residual_once <- max(abs(Z - fixed_point_from_data_term(data_term_once, alpha)))
  residual_twice <- max(abs(Z - fixed_point_from_data_term(data_term_twice, alpha)))

  expect_lt(residual_once, 0.05)
  # The buggy ("twice") fixed point must now be clearly off, otherwise the
  # test would not be discriminating.
  expect_gt(residual_twice, 0.15)
})

test_that("BM_gaussian SBM_sym E step matches the single-count fixed point", {
  set.seed(3)
  n <- 28
  Q <- 2
  cl <- sample(Q, n, TRUE)
  Mu <- matrix(c(0.3, -0.3, -0.3, 0.3), Q, Q)
  A <- matrix(rnorm(n * n, Mu[cl, cl], 1), n, n)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  Zi <- initial_Z(cl, Q)

  r <- blockmodels:::dispatcher("SBM_sym", list(Z = Zi), "gaussian", list(adjacency = A), TRUE)
  Z <- r$membership$Z
  alpha <- r$membership$alpha
  mu <- r$model$mu
  sigma2 <- r$model$sigma2

  off_diag <- 1 - diag(n)
  data_term_once <- (A %*% Z %*% mu) / sigma2 - (off_diag %*% Z %*% (mu * mu)) / (2 * sigma2)
  data_term_twice <- 2 * data_term_once

  residual_once <- max(abs(Z - fixed_point_from_data_term(data_term_once, alpha)))
  residual_twice <- max(abs(Z - fixed_point_from_data_term(data_term_twice, alpha)))

  expect_lt(residual_once, 0.05)
  # The buggy ("twice") fixed point must now be clearly off, otherwise the
  # test would not be discriminating.
  expect_gt(residual_twice, 0.15)
})
