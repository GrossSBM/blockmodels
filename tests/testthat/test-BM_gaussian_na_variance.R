set.seed(20261021)

# Expected sum of squares under the posterior, on the observed entries only:
#   sum_{(i,j) observed} sum_{q,l} Z1[i,q] Z2[j,l] (Y[i,j] - mu[q,l])^2
# divided by the number of observed entries (the diagonal is excluded for SBM)
expected_sigma2 <- function(Y, Z1, Z2, mu, zero_diag) {
  O <- 1 * !is.na(Y)
  if (zero_diag) diag(O) <- 0
  Y[is.na(Y)] <- 0
  SS <- sum(Y^2 * O) + sum(mu^2 * (t(Z1) %*% O %*% Z2)) - 2 * sum(mu * (t(Z1) %*% (Y * O) %*% Z2))
  SS / sum(O)
}

sample_gaussian <- function(n1, n2, Q1, Q2) {
  Z1 <- diag(Q1)[sample(Q1, n1, replace = TRUE), ]
  Z2 <- diag(Q2)[sample(Q2, n2, replace = TRUE), ]
  Mu <- matrix(runif(Q1 * Q2, 0, 5), Q1, Q2)
  Z1 %*% Mu %*% t(Z2) + matrix(rnorm(n1 * n2), n1, n2)
}

for (membership_type in c("SBM", "LBM")) {
  for (with_na in c(FALSE, TRUE)) {
    test_that(paste("BM_gaussian", membership_type, "sigma2 is the mean squared error over the observed entries", if (with_na) "(with NA)" else "(without NA)"), {
      n1 <- 30
      n2 <- if (membership_type == "SBM") n1 else 24
      Y <- sample_gaussian(n1, n2, 2, 2)
      if (with_na) {
        Y[sample(n1 * n2, floor(0.2 * n1 * n2))] <- NA
      }

      model <- BM_gaussian(membership_type, Y, plotting = "", explore_min = 3, explore_max = 3, ncores = 1, verbosity = 0)
      model$estimate()

      for (Q in which(!is.na(model$ICL))) {
        membership <- model$memberships[[Q]]
        parameters <- model$model_parameters[[Q]]
        Z1 <- if (membership_type == "SBM") membership$Z else membership$Z1
        Z2 <- if (membership_type == "SBM") membership$Z else membership$Z2
        expect_equal(
          parameters$sigma2,
          expected_sigma2(Y, Z1, Z2, parameters$mu, zero_diag = membership_type == "SBM"),
          tolerance = 1e-10
        )
      }
    })
  }
}

for (membership_type in c("SBM", "SBM_sym", "LBM")) {
  test_that(paste("BM_gaussian_covariates", membership_type, "ignores the covariates on the missing entries"), {
    n1 <- 30
    n2 <- if (membership_type == "LBM") 24 else n1
    C <- matrix(rnorm(n1 * n2), n1, n2)
    Y <- sample_gaussian(n1, n2, 2, 2) + 2 * C
    if (membership_type == "SBM_sym") {
      Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)]
      C[lower.tri(C)] <- t(C)[lower.tri(C)]
    }
    na <- sample(n1 * n2, floor(0.2 * n1 * n2))
    Y[na] <- NA
    if (membership_type == "SBM_sym") {
      Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)]
    }

    # same data, except the values of the covariate where the network is missing
    C_other <- C
    C_other[is.na(Y)] <- rnorm(sum(is.na(Y)), sd = 10)
    if (membership_type == "SBM_sym") {
      C_other[lower.tri(C_other)] <- t(C_other)[lower.tri(C_other)]
    }

    fit <- function(covariate) {
      set.seed(1)
      model <- BM_gaussian_covariates(membership_type, Y, list(covariate), plotting = "", explore_min = 3, explore_max = 3, ncores = 1, verbosity = 0)
      model$estimate()
      model
    }
    model <- fit(C)
    model_other <- fit(C_other)

    expect_equal(model$ICL, model_other$ICL, tolerance = 1e-10)
    Q <- which.max(model$ICL)
    expect_equal(model$model_parameters[[Q]]$beta, model_other$model_parameters[[Q]]$beta, tolerance = 1e-8)
    expect_equal(model$model_parameters[[Q]]$sigma2, model_other$model_parameters[[Q]]$sigma2, tolerance = 1e-8)
  })
}
