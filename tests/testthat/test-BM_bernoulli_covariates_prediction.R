set.seed(20261020)

sigmo <- function(x) 1 / (1 + exp(-x))

# E[Y_ij] = sum_{q,l} Z1[i,q] Z2[j,l] sigmo(m[q,l] + B[i,j]), entry by entry
brute_force_prediction <- function(Z1, Z2, m, B) {
  P <- matrix(0, nrow(B), ncol(B))
  for (i in seq_len(nrow(B))) {
    for (j in seq_len(ncol(B))) {
      P[i, j] <- sum(outer(Z1[i, ], Z2[j, ]) * sigmo(m + B[i, j]))
    }
  }
  P
}

soft_memberships <- function(n, Q) {
  Z <- matrix(rgamma(n * Q, shape = 0.5), n, Q)
  Z / rowSums(Z)
}

test_that("bernoulli covariates prediction is exact for soft memberships", {
  Z1 <- soft_memberships(15, 3)
  Z2 <- soft_memberships(12, 4)
  m <- matrix(rnorm(12, sd = 2), 3, 4)
  B <- matrix(rnorm(15 * 12, sd = 2), 15, 12)

  expected <- brute_force_prediction(Z1, Z2, m, B)
  expect_equal(blockmodels:::bernoulli_covariates_prediction(Z1, Z2, m, B), expected, tolerance = 1e-12)

  # the previous formula, sigmo(Z1 m Z2' + B), is not the expectation
  expect_gt(max(abs(sigmo(Z1 %*% m %*% t(Z2) + B) - expected)), 1e-2)
})

test_that("bernoulli covariates prediction matches the plug-in formula for hard memberships", {
  Z1 <- diag(3)[sample(3, 15, replace = TRUE), ]
  Z2 <- diag(2)[sample(2, 12, replace = TRUE), ]
  m <- matrix(rnorm(6, sd = 2), 3, 2)
  B <- matrix(rnorm(15 * 12), 15, 12)

  expect_equal(
    blockmodels:::bernoulli_covariates_prediction(Z1, Z2, m, B),
    sigmo(Z1 %*% m %*% t(Z2) + B),
    tolerance = 1e-12
  )
})

for (model_class in c("BM_bernoulli_covariates", "BM_bernoulli_covariates_fast")) {
  for (membership_type in c("SBM", "LBM")) {
    test_that(paste(model_class, membership_type, "prediction() is the expectation under the posterior"), {
      n1 <- 20
      n2 <- if (membership_type == "SBM") n1 else 16
      M <- 1 * (matrix(runif(n1 * n2), n1, n2) < 0.3)
      Y <- matrix(rnorm(n1 * n2), n1, n2)

      model <- get(model_class)(membership_type, M, list(Y), plotting = "", explore_min = 4, explore_max = 4, ncores = 1, verbosity = 0)
      model$estimate()

      # the largest explored model (at least 2 blocks on one side), whose posterior
      # probabilities (almost indicators after the estimation) are replaced by soft
      # ones, for which the plug-in formula is not exact
      Q <- max(which(vapply(model$model_parameters, function(p) length(p$m) > 1, logical(1))))
      parameters <- model$model_parameters[[Q]]
      Z1 <- soft_memberships(n1, nrow(parameters$m))
      Z2 <- soft_memberships(n2, ncol(parameters$m))
      if (membership_type == "SBM") {
        Z2 <- Z1
        model$memberships[[Q]] <- getRefClass("SBM")(from_cc = list(Z = Z1))
      } else {
        model$memberships[[Q]] <- getRefClass("LBM")(from_cc = list(Z1 = Z1, Z2 = Z2))
      }

      expect_equal(
        model$prediction(Q),
        brute_force_prediction(Z1, Z2, parameters$m, parameters$beta[1] * Y),
        tolerance = 1e-12
      )
    })
  }
}
