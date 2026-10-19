set.seed(20261019)

sample_sbm_with_nodes_covariates <- function(n, Q, p, symmetric) {
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  B <- cbind(matrix(rnorm(p * (Q - 1), sd = 1.5), p, Q - 1), 0)
  Z <- t(apply(X %*% B, 1, function(eta) rmultinom(1, 1, exp(eta) / sum(exp(eta)))))
  P <- matrix(0.05, Q, Q)
  diag(P) <- 0.5
  M <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z))
  if (symmetric) {
    M[lower.tri(M)] <- t(M)[lower.tri(M)]
  }
  diag(M) <- 0
  list(M = M, X = X)
}

for (membership_type in c("SBM", "SBM_sym")) {
  test_that(paste("ICL penalty of", membership_type, "accounts for nodes covariates"), {
    n <- 60
    p <- 3
    data <- sample_sbm_with_nodes_covariates(n, Q = 3, p = p, symmetric = membership_type == "SBM_sym")

    model <- BM_bernoulli(
      membership_type, data$M,
      nodes_covariates = list(node = data$X),
      plotting = "", explore_min = 3, explore_max = 3, ncores = 1, verbosity = 0
    )
    model$estimate()

    for (Q in 2:3) {
      membership <- model$memberships[[Q]]
      expect_equal(nrow(membership$B), p)
      expect_equal(membership$ICL_penalty(), p * (Q - 1) * log(n))
    }
  })
}

test_that("ICL penalty of SBM and SBM_sym is unchanged without nodes covariates", {
  n <- 30
  Z <- diag(2) %x% matrix(1, n / 2, 1)
  for (membership_type in c("SBM", "SBM_sym")) {
    membership <- getRefClass(membership_type)(from_cc = list(Z = Z))
    expect_equal(membership$ICL_penalty(), (2 - 1) * log(n))
  }
})
