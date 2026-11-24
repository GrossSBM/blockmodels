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
P1<-matrix(runif(Q*Q),Q,Q) # connectivity matrix of the first layer
P2<-matrix(runif(Q*Q),Q,Q) # connectivity matrix of the second layer
M1<-1*(matrix(runif(n*n),n,n)<Z%*%P1%*%t(Z)) ## adjacency matrices
M2<-1*(matrix(runif(n*n),n,n)<Z%*%P2%*%t(Z))

## estimation
my_model <- BM_bernoulli_multiplex_independent("SBM",list(M1,M2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
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
P1<-matrix(runif(Q*Q),Q,Q)
P2<-matrix(runif(Q*Q),Q,Q)
P1[lower.tri(P1)]<-t(P1)[lower.tri(P1)]
P2[lower.tri(P2)]<-t(P2)[lower.tri(P2)]
M1<-1*(matrix(runif(n*n),n,n)<Z%*%P1%*%t(Z)) ## adjacency matrices
M2<-1*(matrix(runif(n*n),n,n)<Z%*%P2%*%t(Z))
M1[lower.tri(M1)]<-t(M1)[lower.tri(M1)]
M2[lower.tri(M2)]<-t(M2)[lower.tri(M2)]

## estimation
my_model <- BM_bernoulli_multiplex_independent("SBM_sym",list(M1,M2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
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
P1<-matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
P2<-matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
M1<-1*(matrix(runif(n[1]*n[2]),n[1],n[2])<Z1%*%P1%*%t(Z2)) ## adjacency matrices
M2<-1*(matrix(runif(n[1]*n[2]),n[1],n[2])<Z1%*%P2%*%t(Z2))

## estimation
my_model <- BM_bernoulli_multiplex_independent("LBM",list(M1,M2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)
