# blockmodels 1.2.1

* Beginning adding handling of `NA` values in the data. Only scalar (with or without covariates) bernoulli, poisson and gaussian emission distribution models are supported for now. The `NA` values are ignored in the likelihood calculations, and the corresponding entries in the data are not used to update the parameters of the model.
* Added a `NEWS.md` file to track changes to the package.
