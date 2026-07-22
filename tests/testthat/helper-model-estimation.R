library(blockmodels)

expect_model_estimation <- function(model) {
  testthat::expect_no_error(model$estimate())
  testthat::expect_true(is.numeric(model$ICL))
  testthat::expect_gte(length(model$ICL), 1)
  testthat::expect_type(which.max(model$ICL), "integer")
}

expect_model_nodes_covariates <- function(model, expected_nodes_covariates) {
  testthat::expect_true(length(model$memberships) > 0)
  testthat::expect_true(all(vapply(
    model$memberships,
    function(membership) {
      !is.null(membership$nodes_covariates) && identical(membership$nodes_covariates, expected_nodes_covariates)
    },
    logical(1)
  )))
}

expect_model_row_col_covariates <- function(model, expected_row_covariates, expected_col_covariates) {
  testthat::expect_true(length(model$memberships) > 0)
  testthat::expect_true(all(vapply(
    model$memberships[-1], # Because for LBM the first element is not defined
    function(membership) {
      !is.null(membership$row_covariates) &&
        !is.null(membership$col_covariates) &&
        identical(membership$row_covariates, expected_row_covariates) &&
        identical(membership$col_covariates, expected_col_covariates)
    },
    logical(1)
  )))
}
