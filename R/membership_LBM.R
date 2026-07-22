# membership definition

setRefClass("LBM",
    fields = list(
        Z1 = "matrix",
        Z2 = "matrix",
        row_covariates = "matrix",
        col_covariates = "matrix",
        B = "matrix",
        G = "matrix"
    ),
    methods = list(
        initialize = function(network_size = NULL, classif = NULL, from_cc = NULL, row_covar = NULL, col_covar = NULL) {
            if (length(classif) == 0) {
                if (length(network_size) > 0) {
                    Z1 <<- matrix(1, nrow = network_size[1], ncol = 1)
                    Z2 <<- matrix(1, nrow = network_size[2], ncol = 1)
                } else {
                    Z1 <<- from_cc$Z1
                    Z2 <<- from_cc$Z2
                    if ("row_covariates" %in% names(from_cc) && !is.null(from_cc[["row_covariates"]])) {
                        row_covariates <<- from_cc[["row_covariates"]]
                    }
                    if ("B" %in% names(from_cc) && length(from_cc[["B"]]) > 0) {
                        B <<- from_cc[["B"]]
                    }
                    if ("col_covariates" %in% names(from_cc) && !is.null(from_cc[["col_covariates"]])) {
                        col_covariates <<- from_cc[["col_covariates"]]
                    }
                    if ("G" %in% names(from_cc) && length(from_cc[["G"]]) > 0) {
                        G <<- from_cc[["G"]]
                    }
                }
            } else {
                fclassif1 <- factor(classif[[1]])
                classif1 <- as.numeric(fclassif1)

                fclassif2 <- factor(classif[[2]])
                classif2 <- as.numeric(fclassif2)

                Q1 <- length(levels(fclassif1))
                Q2 <- length(levels(fclassif2))

                Z1 <<- matrix(0, nrow = length(classif1), ncol = Q1)
                Z2 <<- matrix(0, nrow = length(classif2), ncol = Q2)

                for (i in seq_along(classif1))
                {
                    Z1[i, classif1[i]] <<- 1
                }

                for (i in seq_along(classif2))
                {
                    Z2[i, classif2[i]] <<- 1
                }
            }
            # TODO Add a check in user function to assert nodes_covariates size, where defined is correct
            if (!is.null(row_covar)) {
                row_covariates <<- row_covar
            }
            if (!is.null(col_covar)) {
                col_covariates <<- col_covar
            }

            install_membership_alpha_binding(.self, "alpha1", "Z1", "row_covariates", "B")
            install_membership_alpha_binding(.self, "alpha2", "Z2", "col_covariates", "G")
        },
        digest = function() {
            digest::digest(
                list(
                    order_round_matrix(Z1),
                    order_round_matrix(Z2)
                ),
                algo = "sha256"
            )
        },
        show_short = function() {
            paste(ncol(Z1), "row groups,", ncol(Z2), "col groups")
        },
        show = function() {
            cat("LBM membership\n")
            cat(paste("    Groups:", ncol(Z1), "row groups,", ncol(Z2), "col groups\n"))
            cat(paste("    Nodes:", nrow(Z1), "row nodes,", nrow(Z2), "col nodes\n"))
            cat("    Useful fields and methods:\n")
            cat("        $Z1 : matrix of row nodes memberships\n")
            cat("        $Z2 : matrix of col nodes memberships\n")
            if (length(row_covariates) > 0) {
                cat("        $row_covariates : matrix of row nodes covariates\n")
            }
            if (length(col_covariates) > 0) {
                cat("        $col_covariates : matrix of col nodes covariates\n")
            }
            cat("        $plot() : plot the memberships\n")
        },
        to_cc = function() {
            output_list <- list(Z1 = Z1, Z2 = Z2)

            if (length(row_covariates) > 0) {
                output_list[["row_covariates"]] <- row_covariates
            }
            if (length(col_covariates) > 0) {
                output_list[["col_covariates"]] <- col_covariates
            }
            if (length(B) > 0) {
                output_list[["B"]] <- B
            }
            if (length(G) > 0) {
                output_list[["G"]] <- G
            }
            return(output_list)
        },
        map = function() {
            list(
                C1 = apply(Z1, 1, which.max),
                C2 = apply(Z2, 1, which.max)
            )
        },
        ICL_penalty = function() {
            if (length(B) > 0) {
                p <- nrow(B)
            } else {
                p <- 1
            }
            if (length(G) > 0) {
                q <- nrow(G)
            } else {
                q <- 1
            }
            p * (dim(Z1)[2] - 1) * log(dim(Z1)[1]) + q * (dim(Z2)[2] - 1) * log(dim(Z2)[1])
        },
        merges = function(merge1, merge2) {
            result <- list()
            Q1 <- dim(Z1)[2]
            Q2 <- dim(Z2)[2]

            if (Q1 > 1 && merge1) {
                for (k1 in 1:(Q1 - 1))
                {
                    for (k2 in (k1 + 1):Q1)
                    {
                        Zn <- as.matrix(Z1[, -k2])
                        Zn[, k1] <- Z1[, k1] + Z1[, k2]
                        result <- c(result, list(getRefClass("LBM")(from_cc = list(Z1 = Zn, Z2 = Z2, row_covariates = row_covariates, col_covariates = col_covariates), row_covar = row_covariates, col_covar = col_covariates)))
                    }
                }
            }

            if (Q2 > 1 && merge2) {
                for (k1 in 1:(Q2 - 1))
                {
                    for (k2 in (k1 + 1):Q2)
                    {
                        Zn <- as.matrix(Z2[, -k2])
                        Zn[, k1] <- Z2[, k1] + Z2[, k2]
                        result <- c(result, list(getRefClass("LBM")(from_cc = list(Z1 = Z1, Z2 = Zn, row_covariates = row_covariates, col_covariates = col_covariates), row_covar = row_covariates, col_covar = col_covariates)))
                    }
                }
            }

            return(result)
        },
        plot = function() {
            par(mfrow = c(1, 2))

            rn1 <- rownames(Z1)
            if (is.null(rn1)) {
                rn1 <- seq_len(nrow(Z1))
            }
            ordering1 <- order(.self$map()$C1)
            matrixplot(as.matrix(Z1[ordering1, ]), rowlabels = rn1[ordering1])

            rn2 <- rownames(Z2)
            if (is.null(rn2)) {
                rn2 <- seq_len(nrow(Z2))
            }
            ordering2 <- order(.self$map()$C2)
            matrixplot(as.matrix(Z2[ordering2, ]), rowlabels = rn2[ordering2])
        }
    )
)
