# Changelog

## blockmodels 1.2.1

- Added covariates on the nodes, driving the block memberships through a
  multinomial logistic regression (argument `nodes_covariates` of all
  the models: a list with one matrix named `node` for SBM and SBM_sym,
  or with matrices named `row` and/or `col` for LBM). The estimated
  coefficients are stored in the field `B` (and `G` for the columns of
  an LBM) of the memberships, and the ICL penalty accounts for them.
- Beginning adding handling of `NA` values in the data. Only scalar
  (with or without covariates) bernoulli, poisson and gaussian emission
  distribution models are supported for now. The `NA` values are ignored
  in the likelihood calculations, and the corresponding entries in the
  data are not used to update the parameters of the model. Multivariate
  models, and missing values in the edge covariates, now stop with an
  explicit error.
- Added a `NEWS.md` file to track changes to the package.
