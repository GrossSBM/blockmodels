require('blockmodels')
set.seed(12)
##
## SBM
##

## generation of one SBM network
npc <- 10 # nodes per class
Q <- 2 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
M<-matrix(rnorm(n*n,sd=10),n,n)+Z%*%Mu%*%t(Z) ## adjacency matrix

## estimation
my_model <- BM_gaussian("SBM",M , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

##
## SBM symmetric
##

## generation of one SBM_sym network
npc <- 10 # nodes per class
Q <- 2 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
Mu[lower.tri(Mu)]<-t(Mu)[lower.tri(Mu)]
M<-matrix(rnorm(n*n,sd=10),n,n)+Z%*%Mu%*%t(Z) ## adjacency matrix
M[lower.tri(M)]<-t(M)[lower.tri(M)]

## estimation
my_model <- BM_gaussian("SBM_sym",M , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

##
## LBM
##

## generation of one LBM network
npc <- c(20,10) # nodes per class
Q <- c(1,2) # classes
n <- npc * Q # nodes
Z1<-diag(Q[1])%x%matrix(1,npc[1],1)
Z2<-diag(Q[2])%x%matrix(1,npc[2],1)
Mu<-20*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
M<-matrix(rnorm(n[1]*n[2],sd=10),n[1],n[2])+Z1%*%Mu%*%t(Z2) ## adjacency matrix

## estimation
my_model <- BM_gaussian("LBM",M , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

##
## NA handling (SBM)
##

npc_na <- 10
Q_na <- 2
n_na <- npc_na * Q_na
Z_na <- diag(Q_na) %x% matrix(1, npc_na, 1)
Mu_na <- 20 * matrix(runif(Q_na * Q_na), Q_na, Q_na)
M_na <- matrix(rnorm(n_na * n_na, sd=10), n_na, n_na) + Z_na %*% Mu_na %*% t(Z_na)

off_diag <- which(row(M_na) != col(M_na))
set.seed(123)
M_na[sample(off_diag, max(1, floor(0.1 * length(off_diag))))] <- NA

my_model_na <- BM_gaussian("SBM", M_na, plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model_na$estimate()
stopifnot(length(my_model_na$ICL) > 0)
stopifnot(all(is.finite(my_model_na$ICL)))

M_all_na <- matrix(NA_real_, nrow(M_na), ncol(M_na))
has_error <- FALSE
tryCatch(
	{
		my_model_all_na <- BM_gaussian("SBM", M_all_na, plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
		my_model_all_na$estimate()
	},
	error = function(e)
	{
		has_error <<- TRUE
	}
)
stopifnot(has_error)
