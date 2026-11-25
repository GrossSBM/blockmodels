# The ICL penalty of multivariate models must use log(number of dyads),
# not log(number of layers * number of dyads).

test_that("multivariate data_number counts dyads, not layer entries", {
  set.seed(1)
  n1 <- 12; n2 <- 9; L <- 3
  adj_sq <- lapply(1:L, function(l) matrix(rbinom(n1 * n1, 1, .3), n1, n1))
  adj_sym <- lapply(adj_sq, function(m) { m[lower.tri(m)] <- t(m)[lower.tri(m)]; m })
  adj_rect <- lapply(1:L, function(l) matrix(rbinom(n1 * n2, 1, .3), n1, n2))

  for (ctor in c("BM_bernoulli_multiplex", "BM_bernoulli_multiplex_independent")) {
    f <- get(ctor)
    expect_equal(f("SBM", adj_sq, verbosity = 0)$data_number(), n1 * (n1 - 1))
    expect_equal(f("SBM_sym", adj_sym, verbosity = 0)$data_number(), n1 * (n1 - 1) / 2)
    expect_equal(f("LBM", adj_rect, verbosity = 0)$data_number(), n1 * n2)
  }
})

test_that("multiplex ICL uses log(n(n-1)) in its penalty", {
  set.seed(2)
  n <- 30; L <- 2
  adj <- lapply(1:L, function(l) matrix(rbinom(n * n, 1, .3), n, n))
  for (ctor in c("BM_bernoulli_multiplex", "BM_bernoulli_multiplex_independent")) {
    m <- get(ctor)("SBM", adj, verbosity = 0, plotting = "", explore_min = 3, explore_max = 3)
    m$estimate()
    for (Q in seq_along(m$ICL)) {
      expected <- m$PL[Q] - .5 * (m$model_parameters[[Q]]$n_parameters * log(n * (n - 1)) +
                                    m$memberships[[Q]]$ICL_penalty())
      expect_equal(m$ICL[Q], expected, tolerance = 1e-8)
    }
  }
})

test_that("gaussian_multivariate counts L(L+1)/2 parameters for Sigma", {
  set.seed(3)
  n1 <- 20; n2 <- 15; L <- 3
  adj_sq <- lapply(1:L, function(l) matrix(rnorm(n1 * n1), n1, n1))
  adj_sym <- lapply(adj_sq, function(m) (m + t(m)) / 2)
  adj_rect <- lapply(1:L, function(l) matrix(rnorm(n1 * n2), n1, n2))
  n_mu <- list(
    SBM = function(mb) ncol(mb$Z)^2 * L,
    SBM_sym = function(mb) ncol(mb$Z) * (ncol(mb$Z) + 1) / 2 * L,
    LBM = function(mb) ncol(mb$Z1) * ncol(mb$Z2) * L
  )
  adjs <- list(SBM = adj_sq, SBM_sym = adj_sym, LBM = adj_rect)
  for (membership in names(adjs)) {
    m <- BM_gaussian_multivariate(membership, adjs[[membership]], verbosity = 0,
                                  plotting = "", explore_min = 3, explore_max = 3)
    m$estimate()
    fitted <- which(!vapply(m$memberships, is.null, logical(1)))
    expect_gt(length(fitted), 0)
    for (Q in fitted) {
      expect_equal(m$model_parameters[[Q]]$n_parameters,
                   n_mu[[membership]](m$memberships[[Q]]) + L * (L + 1) / 2)
    }
  }
})
