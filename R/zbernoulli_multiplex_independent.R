
BM_bernoulli_multiplex_independent <-
  setRefClass("BM_bernoulli_multiplex_independent",
              contains = "multivariate_model",
              methods = list(
                initialize = function(membership_type,adj,...)
                {
                  .self$initFields(membership_name = membership_type,
                                   model_name = "bernoulli_multiplex_independent",
                                   adj = adj,
                                   ...)
                  .self$postinit()
                },
                accepts_NA = function() { TRUE },
                model_penalty = function(model)
                {
                  # The layers are independent given the blocks: the BIC-like
                  # penalty is the sum over the layers of
                  # (number of parameters of the layer) * log(number of observed dyads of the layer)
                  n_parameters_per_layer <- model$n_parameters / length(adj)
                  sum(vapply(adj, function(a) n_parameters_per_layer * log(.self$layer_data_number(a)),
                             numeric(1)))
                },
                layer_data_number = function(a)
                {
                  # number of observed dyads in the layer a
                  mask <- !is.na(a)
                  if(membership_name=="SBM")
                  {
                    return(sum(mask) - sum(diag(mask)))
                  }
                  if(membership_name=="SBM_sym")
                  {
                    return((sum(mask) - sum(diag(mask)))/2)
                  }
                  return(sum(mask))
                },
                plot_parameters = function(Q)
                {
                },
                prediction = function(Q)
                {
                  pred <- list()
                  for(k in 1:length(adj))
                  {
                    if(membership_name=='LBM')
                    {
                      pred[[k]] <- memberships[[Q]]$Z1 %*% model_parameters[[Q]]$pi[[k]] %*% t(memberships[[Q]]$Z2)
                    }
                    else
                    {
                      pred[[k]] <- memberships[[Q]]$Z %*% model_parameters[[Q]]$pi[[k]] %*% t(memberships[[Q]]$Z)
                    }
                  }
                  return(pred)
                },
                residual = function(Q)
                {
                  pred <- .self$prediction(Q)
                  ret <- list()
                  for(k in 1:length(adj))
                  {
                    ret[[k]] <- adj[[k]]-pred[[k]]
                    # null residuals on the missing entries: they must not
                    # drive the initializations (splits of the groups)
                    ret[[k]][is.na(ret[[k]])] <- 0
                  }
                  return(ret)
                }
              )
  )
