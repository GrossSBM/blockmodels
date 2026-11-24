
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
                  }
                  return(ret)
                }
              )
  )
