# Changelog

## blockmodels 1.2.1

### New features

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
  models, and missing values in the edge covariates, stop with an
  explicit error.

### Bug fixes that change the results

- Undirected models (`SBM_sym` memberships, i.e. all
  `sbm::estimateSimpleSBM()` calls on a symmetric network, the default
  for symmetric networks): the E step actually executed was the one of
  `SBM`, which counts each dyad twice instead of once. As a consequence,
  the posterior probabilities (`Z`) were overconfident (too sharp) and
  the VEM did not exactly optimize the criterion it reports, although
  the final clustering was mostly unaffected on well-separated data.
  This also exposed (and fixes) two other, previously untested, latent
  bugs that were introduced by the same root cause (an `e_fixed_step()`
  overload resolved on the static, base `SBM` type instead of on the
  actual `SBM_sym` membership type, so the SBM_sym-specific E step code
  was never executed in practice):
  - `BM_gaussian_covariates` with `SBM_sym`: the E step used the wrong
    matrices (the variance term used the adjacency/covariates residual
    matrix instead of the “all ones” matrix, and vice-versa for the
    cross term).
  - `BM_bernoulli_covariates_fast` with `SBM_sym`: the E step was
    missing a `Z` factor, causing a matrix dimension mismatch (this
    model would error out with `SBM_sym` memberships). The `SBM`
    (directed) E step of the same model was also silently dropping half
    of its terms because of a stray semicolon, independently of the
    `SBM_sym` issue; both are fixed.

  Because of this fix, **all fits of a symmetric network (`SBM_sym`,
  including `sbm::estimateSimpleSBM()`) will return different, more
  accurate posterior probabilities and, for `BM_gaussian_covariates` and
  `BM_bernoulli_covariates_fast`, different parameter estimates as
  well.** Model selection and downstream uses of the posteriors may be
  affected, especially when the signal is weak
  ([\#8](https://github.com/GrossSBM/blockmodels/issues/8)).
- Fixed `prediction()` of `BM_bernoulli_covariates` and
  `BM_bernoulli_covariates_fast`: the expected network is now the
  expectation under the posterior probabilities, sum over (q,l) of
  Z\[i,q\] Z\[j,l\] sigmo(m\[q,l\] + covariates effect), instead of
  sigmo(Z m Z’ + covariates effect), which was exact only for indicator
  memberships ([\#2](https://github.com/GrossSBM/blockmodels/issues/2)).
  As the residuals guide the initializations by splitting groups, the
  estimates of these models may slightly differ from the previous
  versions.
- `BM_gaussian` with SBM and SBM_sym memberships: the variance is now
  the mean squared error over the n(n-1) entries out of the diagonal (it
  was divided by n^2), consistently with the likelihood. Its estimate is
  multiplied by n/(n-1) compared to the previous versions.
- Gaussian models with `NA`: the variance and the likelihood use the
  number of observed entries, and the covariates of
  `BM_gaussian_covariates` on the missing entries are ignored. The
  residuals used for the initializations are null on the missing
  entries.

### Other changes

- Faster E and M steps in presence of `NA` (the mask is applied once to
  the precomputed matrices of the network), and no n x n matrix in the M
  step of the nodes covariates.
- Julien Chiquet is the new maintainer of the package; Louis Lacoste
  joins the authors.
- Test suite with testthat, continuous integration and pkgdown website
  on GitHub (<https://github.com/GrossSBM/blockmodels>).
- Added a `NEWS.md` file to track changes to the package.
