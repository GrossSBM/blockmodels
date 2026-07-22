
matrixplot <- function(x,colorbar=FALSE,xlab="",ylab="",title=NULL,rowlabels=NULL,collabels=NULL)
{
    if(is.null(rowlabels))
    {
        rowlabels <- rownames(x)
    }

    if(is.null(collabels))
    {
        collabels <- colnames(x)
    }

    if(is.null(rowlabels))
    {
        rowlabels<-1:nrow(x)
    }

    if(is.null(collabels))
    {
        collabels<-1:ncol(x)
    }

    inds <- nrow(x):1

    image(1:ncol(x),1:nrow(x),as.matrix(t(x)[,inds]),xlab=xlab,ylab=ylab,axes=FALSE,col=hsv(0,0,1-(0:99)/99,1))

    if(!is.null(title))
    {
        title(main=title)
    }

    axis(BELOW<-1, at=1:ncol(x),labels=collabels)
    axis(LEFT<-2, at=1:nrow(x),labels=rowlabels[inds])
}   

order_round_matrix <- function(x)
{
    floor(10*as.matrix(x[,order(colSums(x),decreasing=T)]))
}

cumtime <- function()
{
    sum(
        sapply(
            strsplit(
                system(
                    'cat /proc/$PPID/stat | cut -d " " -f 14-17',
                    intern=TRUE
                ),
                ' '
            ),
            strtoi
        )
    )
}

softmax_rows <- function(x)
{
    x_shift <- x - apply(x, 1, max)
    ex <- exp(x_shift)
    ex / rowSums(ex)
}

membership_compute_alpha <- function(Z, covariates, coefficients)
{
    if (length(covariates) > 0 && length(coefficients) > 0)
    {
        return(softmax_rows(covariates %*% coefficients))
    }

    colMeans(Z)
}

install_membership_alpha_binding <- function(object_env, binding_name, z_field, covariates_field, coefficients_field)
{
    if (!exists(binding_name, envir = object_env, inherits = FALSE))
    {
        makeActiveBinding(
            binding_name,
            function(value)
            {
                if (!missing(value))
                {
                    stop(paste(binding_name, "is read-only"))
                }

                Z <- get(z_field, envir = object_env)
                covariates <- get(covariates_field, envir = object_env)
                coefficients <- get(coefficients_field, envir = object_env)

                membership_compute_alpha(Z, covariates, coefficients)
            },
            object_env
        )
    }
}
