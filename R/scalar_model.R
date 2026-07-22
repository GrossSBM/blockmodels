setRefClass("scalar_model",
    contains = "model",
    fields = list(
        adj = "matrix",
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

            if(membership_name=="SBM" || membership_name=="SBM_sym")
            {
                if(nrow(adj)!=ncol(adj))
                {
                    stop(paste("The adjacency matrix does not have the same number of rows and columns.",
                               "Vocatus periculosum ad sanitatem est."))
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
                    if (nrow(nodes_covariates[["node"]]) != nrow(adj)) {
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
                            number_of_nodes <- ifelse(dim == "row", nrow(adj), ncol(adj))
                            if (nrow(nodes_covariates[[dim]]) != number_of_nodes) {
                                stop(paste0("The number of rows for the ", dim, " nodes covariates matrix must match the number of ", dim, " nodes."))
                            }
                        }
                    })
                }
            }

            if(membership_name=="SBM_sym")
            {
                if(isSymmetric(adj))
                {
                    adj <<- (adj+t(adj))/2
                }
                else
                {
                    stop("Adjacency matrix is not symmetric. You need more sleep.")
                }
            }
        },
        number_of_nodes = function() { dim(adj) },
        show_network = function()
        {
            paste(nrow(adj),"x",ncol(adj),"scalar network")
        },  
        network_to_cc = function() { list(adjacency = adj, na_replace_value = na_replace_value) },
        split_membership_model = function(Q)
        {
            membership <- memberships[[Q]]
            error <- .self$residual(Q)

            if(membership_name == "SBM" || membership_name == "SBM_sym")
            {
                result <- list()
                for(q in 1:Q)
                {
                    sub_classif <- coordinates_split(
                        cbind(error, t(error)),
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
                        sub_classif <- coordinates_split(
                            cbind(error),
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
                        sub_classif <- coordinates_split(
                            cbind(t(error)),
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
        data_number = function()
        {
            mask <- !is.na(adj)
            if(membership_name=="SBM")
            {
                return(sum(mask) - sum(diag(mask)))
            }
            if(membership_name=="SBM_sym")
            {
                return((sum(mask) - sum(diag(mask)))/2)
            }
            else
            {
                return(sum(mask))
            }
        },
        precompute = function()
        {
            if(membership_name == "SBM" || membership_name == "SBM_sym")
            {
                if(length(precomputed)>0)
                {
                    return()
                }
                else
                {
                    if(length(ICL)!=0)
                    {
                        say(1, "Computation of eigen decomposition used for initalizations")
                        
                        error <- .self$residual(1)
                        W<- error %*% t(error)
                        W<-1/(1+exp(-W/sd(W)))
                        D<- diag(1/sqrt(rowSums(W)))
                        L<- D %*% W %*% D

                        precomputed$eigen <<- eigen(L, symmetric=TRUE)
                    }
                }
            }
            if(membership_name == "LBM")
            {
                if(length(precomputed)>0)
                {
                    return()
                }
                else
                {
                    if(length(ICL)!=0)
                    {
                        say(1, "Computation of eigen decomposition used for initalizations")
                        error <- .self$residual(2)
                        
                        say(2,"for rows")
                        
                        W1<- error %*% t(error)
                        W1<-1/(1+exp(-W1/sd(W1)))
                        D1<- diag(1/sqrt(rowSums(W1)))
                        L1<- D1 %*% W1 %*% D1

                        precomputed$eigen1 <<- eigen(L1, symmetric=TRUE)
                        
                        say(2,"for cols")
                        
                        W2<- t(error) %*% error
                        W2<-1/(1+exp(-W2/sd(W2)))
                        D2<- diag(1/sqrt(rowSums(W2)))
                        L2<- D2 %*% W2 %*% D2

                        precomputed$eigen2 <<- eigen(L2, symmetric=TRUE)
                    }
                }
            }
        },
        provide_init = function(Q)
        {
            .self$precompute()
            if(membership_name == "SBM" || membership_name == "SBM_sym")
            {
                return(
                    list(
                        getRefClass(membership_name)(
                            classif=blockmodelskmeans(
                                as.matrix(precomputed$eigen$vectors[,1:Q]),
                                Q
                            ),
                            nodes_covar = .self$nodes_covariates[["node"]]
                        )
                    )
                )
            }
            if(membership_name == "LBM")
            {
                result <- list()
                if(length(exploration_direction)!=0)
                {
                    theo <- Q*exploration_direction/sum(exploration_direction)
                    itheo <- floor(theo)
                    if(sum(itheo)<Q)
                    {
                        added <- which.max(theo-itheo)
                        itheo[added] <- itheo[added]+1
                    }
                }
                for(Q1 in 1:(Q-1))
                {
                    do_this <- TRUE
                    if(length(exploration_direction)!=0)
                    {
                        if(Q1 == itheo[1])
                        {
                            do_this <- TRUE
                        }
                        else
                        {
                            do_this <- FALSE
                        }
                    }
                    if(do_this)
                    {
                        Q2<-Q-Q1
                        if(Q1<=nrow(adj) && Q2<=ncol(adj))
                        {
                            result[[Q1]] <- getRefClass(membership_name)(
                                classif=list(
                                    blockmodelskmeans(as.matrix(precomputed$eigen1$vectors[,1:Q1]),Q1),
                                    blockmodelskmeans(as.matrix(precomputed$eigen2$vectors[,1:Q2]),Q2)
                                ),
                                row_covar = .self$nodes_covariates[["row"]],
                                col_covar = .self$nodes_covariates[["col"]]
                            )
                        }
                    }
                }
                result <- result[!sapply(result,is.null)]

                return(result)
            }
        },
        plot_obs_pred = function(Q)
        {
            pred <- .self$prediction(Q)

            par(mfrow=c(2,1))
            if(membership_name == "LBM")
            {
                order1 <- order(memberships[[Q]]$map()$C1)
                order2 <- order(memberships[[Q]]$map()$C2)
            }
            else
            {
                order1 <- order(memberships[[Q]]$map()$C)
                order2 <- order1
            }

            rn<-rownames(adj)
            cn<-colnames(adj)

            if(is.null(rn))
            {
                rn<-1:nrow(adj)
            }
            if(is.null(cn))
            {
                cn<-1:ncol(adj)
            }

            matrixplot(.self$plot_transform(adj[order1,order2]),rowlabels=rn[order1],collabels=cn[order2])
            
            matrixplot(.self$plot_transform(pred[order1,order2]),rowlabels=rn[order1],collabels=cn[order2])
        },
        plot_transform = function(x){x},
        residual = function(Q)
        {
            non_na_adj <- adj
            non_na_adj[is.na(adj)] <- 0
            non_na_adj-.self$prediction(Q)
        }

            
    )
)


