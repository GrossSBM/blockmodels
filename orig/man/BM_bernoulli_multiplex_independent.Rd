\name{BM_bernoulli_multiplex_independent}
\alias{BM_bernoulli_multiplex_independent}
\alias{\S4method{BM_bernoulli_multiplex_independent}{new}}
\title{Perform estimation on blockmodels for multiplex binary networks with
    independent layers}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership. The network has several binary layers on the
same nodes. Given the blocks of the nodes, the layers are independent, each
one following a Bernoulli blockmodel with its own connectivity matrix: there is
no interaction between the layers (see \code{\link{BM_bernoulli_multiplex}} for
a model of the dependence between the layers). The estimated connectivity
matrices are stored in \code{model_parameters[[Q]]$pi}, a list with one matrix
per layer.}
\usage{\S4method{BM_bernoulli_multiplex_independent}{new}(
    MEMBERSHIP_ARG, 
    MULTIVARIATE_MODEL_ARG, 
    MODEL_ARG)}
\arguments{
    MEMBERSHIP_ARG_ITEM
    MULTIVARIATE_MODEL_NA_ARG_ITEM
    MODEL_ARG_ITEM
}
