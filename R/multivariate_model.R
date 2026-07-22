
setRefClass("multivariate_model",
    contains = "model",
    fields = list(
        adj = "list",
        nodes_covariates = "list",
        na_replace_value = "numeric"
    ),
    methods = list(
        postinit = function()
        {
            callSuper()

            if (length(na_replace_value) == 0)
            {
                na_replace_value <<- 0
            }
            if (length(na_replace_value) != 1 || !is.finite(na_replace_value))
            {
                stop("na_replace_value must be a single finite numeric value.")
            }

            if(length(adj)<1)
            {
                stop(paste("The adjacency list must have at least one matrix",
                           "Inanis vacuum est."))
            }

            for(i in 1:length(adj))
            {
                if(!all(dim(adj[[1]])==dim(adj[[i]])) || length(dim(adj[[i]]))!=2)
                {
                    stop(paste("All adjacencies matrix must have the same size"))
                }
                if(membership_name=="SBM" || membership_name=="SBM_sym")
                {
                    if(nrow(adj[[i]])!=ncol(adj[[i]]))
                    {
                        stop(paste("The adjacency matrix",i,"does not have the same number of rows and columns.",
                                   "Furibunda matrix.."))
                    }
                }

                if(membership_name=="SBM_sym")
                {
                    if(isSymmetric(adj[[i]]))
                    {
                        adj[[i]] <<- (adj[[i]]+t(adj[[i]]))/2
                    }
                    else
                    {
                        stop("Adjacency matrix",i,"is not symmetric. You need more coffee.")
                    }
                }
            }

            if (length(nodes_covariates) > 0) {
                if (any(!sapply(nodes_covariates, is.matrix))) {
                    stop(paste("Nodes covariates must be matrices."))
                }

                # TODO add checks for the values of the covariates

                if ((membership_name == "SBM_sym" || membership_name == "SBM")) {
                    if (length(nodes_covariates) > 1) {
                        stop(paste("Multiple nodes covariates given for SBM.", "Should only be a list with one matrix."))
                    }
                    if (nrow(nodes_covariates[["node"]]) != nrow(adj[[1]])) {
                        stop(paste("The number of rows of the node covariates matrix must match the number of nodes."))
                    }
                    if (is.null(names(nodes_covariates)) || names(nodes_covariates) != "node") {
                        stop(paste("For SBM node covariates matrix must be named 'node'."))
                    }
                }
                if (membership_name == "LBM") {
                    if (is.null(names(nodes_covariates))) {
                        stop(paste("For LBM node covariates matrices must be named (row, col) to indicate which nodes the covariates are on."))
                    }
                    if (any(!(names(nodes_covariates) %in% c("row", "col")))) {
                        stop(paste("For LBM node covariates matrices, the names must be either row or col."))
                    }

                    sapply(c("row", "col"), function(dim) {
                        if (dim %in% names(nodes_covariates)) {
                            number_of_nodes <- ifelse(dim == "row", nrow(adj[[1]]), ncol(adj[[1]]))
                            if (nrow(nodes_covariates[[dim]]) != number_of_nodes) {
                                stop(paste0("The number of rows for the ", dim, " nodes covariates matrix must match the number of ", dim, " nodes."))
                            }
                        }
                    })
                }
            }
        },
        number_of_nodes = function() { dim(adj[[1]]) },
        show_network = function()
        {
            paste(nrow(adj[[1]]),"x",ncol(adj[[1]]),"multivariate network in dimention",length(adj))
        },  
        network_to_cc = function() { list(adjacency = adj, na_replace_value = na_replace_value) },
        data_number = function()
        {
            if(membership_name=="SBM")
            {
                return(dim(adj[[1]])[1]*(dim(adj[[1]])[1]-1)*length(adj))
            }
            if(membership_name=="SBM_sym")
            {
                return(dim(adj[[1]])[1]*(dim(adj[[1]])[1]-1)/2*length(adj))
            }
            else
            {
                return(dim(adj[[1]])[1]*dim(adj[[1]])[2]*length(adj))
            }
        },
        split_membership_model = function(Q)
        {
            membership <- memberships[[Q]]
            error <- .self$residual(Q)

            if(membership_name == "SBM" || membership_name == "SBM_sym")
            {
                result <- list()
                for(q in 1:Q)
                {
                    allcordsbind <- cbind(error[[1]], t(error[[1]]))
                    if(length(adj)>1)
                    {
                        for(k in 2:length(adj))
                        {
                            allcordsbind <- cbind(allcordsbind, error[[k]], t(error[[k]]))
                        }
                    }
                    sub_classif <- coordinates_split(
                        allcordsbind,
                        membership$Z[,q]
                        )
                    Z <- cbind(membership$Z,membership$Z[,q])
                    Z[,q] <- Z[,q]*sub_classif
                    Z[,Q+1] <- Z[,Q+1]*(1-sub_classif)
                    result <- c(result, list(
                            getRefClass(membership_name)(from_cc=list(Z=Z), nodes_covar = .self$nodes_covariates[["node"]])
                        ))
                }
                return(result)
            }
            if(membership_name == "LBM")
            {

                Q1<-dim(membership$Z1)[2]
                Q2<-dim(membership$Z2)[2]

                split1 <- TRUE
                split2 <- TRUE
                if(length(exploration_direction)!=0)
                {
                    if(Q1<exploration_direction[1] || Q2<exploration_direction[2])
                    {
                        if(Q1/exploration_direction[1] < Q2/exploration_direction[2])
                        {
                            split1 <- TRUE
                            split2 <- FALSE
                        }
                        else
                        {
                            split1 <- FALSE
                            split2 <- TRUE
                        }
                    }
                }

                result <- list()
                if(split1)
                {
                    for(q in 1:Q1)
                    {
                        allcordsbind <- cbind(error[[1]])
                        if(length(adj)>1)
                        {
                            for(k in 2:length(adj))
                            {
                                allcordsbind <- cbind(allcordsbind, error[[k]])
                            }
                        }
                        sub_classif <- coordinates_split(
                            allcordsbind,
                            membership$Z1[,q]
                            )
                        Z1 <- cbind(membership$Z1,membership$Z1[,q])
                        Z1[,q] <- Z1[,q]*sub_classif
                        Z1[,Q1+1] <- Z1[,Q1+1]*(1-sub_classif)
                        result <- c(result, list(
                            getRefClass(membership_name)(from_cc=list(Z1=Z1,Z2=membership$Z2, row_covariates =.self$nodes_covariates[["row"]],
                            col_covariates=.self$nodes_covariates[["col"]]))
                            ))
                    }
                }

                if(split2)
                {
                    for(q in 1:Q2)
                    {
                        allcordsbind <- cbind(t(error[[1]]))
                        if(length(adj)>1)
                        {
                            for(k in 2:length(adj))
                            {
                                allcordsbind <- cbind(allcordsbind, t(error[[k]]))
                            }
                        }
                        sub_classif <- coordinates_split(
                            allcordsbind,
                            membership$Z2[,q]
                            )
                        Z2 <- cbind(membership$Z2,membership$Z2[,q])
                        Z2[,q] <- Z2[,q]*sub_classif
                        Z2[,Q2+1] <- Z2[,Q2+1]*(1-sub_classif)
                        result <- c(result, list(
                            getRefClass(membership_name)(from_cc=list(Z1=membership$Z1,Z2=Z2, row_covariates=.self$nodes_covariates[["row"]],
                            col_covariates=.self$nodes_covariates[["col"]]))
                            ))
                    }
                }
                return(result)
            }
        },
        provide_init = function(Q)
        {
            return(list())
        },
        plot_obs_pred = function(Q)
        {
        },
        plot_transform = function(x){x}
    )
)


