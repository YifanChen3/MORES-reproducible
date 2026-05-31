library(Rcpp)
library(MASS)
library(Matrix)
library(RSpectra)
library(glmmsel)
library(MOMENT)

multi.glmmsel.sensitivity <- function(Y,X,id,lambda.path,p,q,sigmaB=NULL,beta=NULL){
  nis <- tabulate(factor(id,levels = unique(id)))
  m <- length(nis)
  N <- sum(nis)
  d <- ncol(Y)
  sigmaB.list <- NULL
  sigmaE.list <- NULL
  beta.final <- matrix(0,p,d)
  lambda.optimal <- NULL
  #group <- rep(1:m,each=ni)
  for (di in 1:d) {
    y <- Y[,di]
    
    fit <- cv.glmmsel(X, y, id, family = "gaussian",nfold = 5,lambda = lambda.path)
    
    lambda.index <- which(fit$lambda==fit$lambda.min)
    sigmaB.list[[di]] <- diag(c(fit$fit$gamma0[lambda.index],fit$fit$gamma[c(p:(p+q-2)),lambda.index]))
    beta.final[,di] <- as.matrix(c(fit$fit$beta0[lambda.index],fit$fit$beta[c(1:(p-1)),lambda.index]))
    lambda.optimal[di] <- fit$lambda.min
    sigmaE.list[[di]] <- fit$fit$sigma2[lambda.index]
  }
  sigmaB.final <- bdiag(sigmaB.list)
  sigmaE.final <- bdiag(sigmaE.list)
  
  if (is.null(sigmaB)==TRUE | (is.null(beta) == TRUE)){
    return(list(lambda.optimal=lambda.optimal,
                sigmaB.final=sigmaB.final,beta.final=beta.final,sigmaE.final=sigmaE.final))
  }
  else{
    fnorm.B <- norm(sigmaB-sigmaB.final,"F")
    diagB.true <- diag(sigmaB) != 0
    diagB.pred <- diag(sigmaB.final) != 0
    F1.B <- 2*sum(diagB.pred == 1 & diagB.true == 1)/(2*sum(diagB.pred == 1 & diagB.true == 1)+sum(diagB.pred == 1 & diagB.true == 0)+sum(diagB.pred == 0 & diagB.true == 1))
    
    fnorm.beta <- norm(beta-beta.final,"F")
    beta.true <- as.vector(beta) != 0
    beta.pred <- as.vector(beta.final) != 0
    F1.beta <- 2*sum(beta.pred == 1 & beta.true == 1)/(2*sum(beta.pred == 1 & beta.true == 1)+sum(beta.pred == 1 & beta.true == 0)+sum(beta.pred == 0 & beta.true == 1))
    
    return(list(fnorm.B=fnorm.B,F1.B=F1.B,fnorm.beta=fnorm.beta,F1.beta=F1.beta,
                lambda.optimal=lambda.optimal,
                sigmaB.final=sigmaB.final,beta.final=beta.final,
                true.sigmaB=sigmaB,true.beta=beta))
  }
}


########## MOMENT,q=15,s=6 ###########
sim.fun <- function(seed){
  d <- 5
  q <- 15
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 6
    rho1 <- 0.5
    rho2 <- 0.6
    rho3 <- 0.7
    rho4 <- 0.8
    rho5 <- 0.9
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    sigmaB <- bdiag(M1,M2,M3,M4,M5)
    rho <- 0.75
    sigmaE <- diag(c(1,0.9,0.8,0.7,0.6))
    lambda.max <- 90000
    name <- "independent"
  }
  
  if (seed>=51 & seed<=100){
    s <- 6
    rho1 <- 0.5
    rho2 <- 0.5
    rho3 <- 0.5
    rho4 <- 0.5
    rho5 <- 0.5
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    R1 <- cbind(M1,M2,M3,M4,M5)
    R2 <- cbind(M2,M1,M2,M3,M4)
    R3 <- cbind(M3,M2,M1,M2,M3)
    R4 <- cbind(M4,M3,M2,M1,M2)
    R5 <- cbind(M5,M4,M3,M2,M1)
    sigmaB <- rbind(R1,R2,R3,R4,R5)
    rho <- 0.75
    sigmaE <- toeplitz(rho^seq(0, d - 1))
    lambda.max <- 90000
    name <- "correlated"
  }
  
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  # Gaussian B
  vecB <- mvrnorm(m,rep(0,dq),sigmaB)
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  # Below is for Xi=Zi with intercept (high d)
  Z <- cbind(1,matrix(rnorm(N*(q-1)),nrow = N,ncol = q-1))
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- cbind(1,matrix(rnorm(N*(p-1)),nrow = N,ncol = p-1))
  beta <- matrix(c(-2:2,rep(0,15),2,1,1,-2,-1,rep(0,15),-1,2,1,2,-2,rep(0,15),-2,2,1,1,2,rep(0,15),1,2,2,-2,-1,rep(0,15)),nrow = p,ncol = d)
  
  Y <- X %*% beta + as.matrix(Z_long) %*% B + epsilon
  
  lambda.min <- lambda.max/100
  length.out <- 100
  lambdas <- exp(seq(log(lambda.max),log(lambda.min),length.out=length.out))
  
  tau.max <- 10
  tau.min <- tau.max/100
  length.out <- 100
  taus <- seq(tau.max,tau.min,length.out=length.out)
  
  t1 <- Sys.time()
  result.sim <- MOMENT.cv(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/sensitivity/our_full_high_1_4/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}


########## MOMENT,q=8,s=4 ###########
sim.fun <- function(seed){
  d <- 5
  q <- 8
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 4
    rho1 <- 0.5
    rho2 <- 0.6
    rho3 <- 0.7
    rho4 <- 0.8
    rho5 <- 0.9
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    sigmaB <- bdiag(M1,M2,M3,M4,M5)
    rho <- 0.75
    sigmaE <- diag(c(1,0.9,0.8,0.7,0.6))
    lambda.max <- 60000
    name <- "independent"
  }
  
  if (seed>=51 & seed<=100){
    s <- 4
    rho1 <- 0.5
    rho2 <- 0.5
    rho3 <- 0.5
    rho4 <- 0.5
    rho5 <- 0.5
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    R1 <- cbind(M1,M2,M3,M4,M5)
    R2 <- cbind(M2,M1,M2,M3,M4)
    R3 <- cbind(M3,M2,M1,M2,M3)
    R4 <- cbind(M4,M3,M2,M1,M2)
    R5 <- cbind(M5,M4,M3,M2,M1)
    sigmaB <- rbind(R1,R2,R3,R4,R5)
    rho <- 0.75
    sigmaE <- toeplitz(rho^seq(0, d - 1))
    lambda.max <- 60000
    name <- "correlated"
  }
  
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  # Gaussian B
  vecB <- mvrnorm(m,rep(0,dq),sigmaB)
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  # Below is for Xi=Zi with intercept (high d)
  Z <- cbind(1,matrix(rnorm(N*(q-1)),nrow = N,ncol = q-1))
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- cbind(1,matrix(rnorm(N*(p-1)),nrow = N,ncol = p-1))
  beta <- matrix(c(-2:2,rep(0,15),2,1,1,-2,-1,rep(0,15),-1,2,1,2,-2,rep(0,15),-2,2,1,1,2,rep(0,15),1,2,2,-2,-1,rep(0,15)),nrow = p,ncol = d)
  
  Y <- X %*% beta + as.matrix(Z_long) %*% B + epsilon
  
  lambda.min <- lambda.max/100
  length.out <- 100
  lambdas <- exp(seq(log(lambda.max),log(lambda.min),length.out=length.out))
  
  tau.max <- 10
  tau.min <- tau.max/100
  length.out <- 100
  taus <- seq(tau.max,tau.min,length.out=length.out)
  
  t1 <- Sys.time()
  result.sim <- MOMENT.cv(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/sensitivity/our_full_low_1_4/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

########## MOMENT(Marginal),q=15,s=6 ###########
sim.fun <- function(seed){
  d <- 5
  q <- 15
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 6
    rho1 <- 0.5
    rho2 <- 0.6
    rho3 <- 0.7
    rho4 <- 0.8
    rho5 <- 0.9
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    sigmaB <- bdiag(M1,M2,M3,M4,M5)
    rho <- 0.75
    sigmaE <- diag(c(1,0.9,0.8,0.7,0.6))
    lambda.max <- 90000
    name <- "independent"
  }
  
  if (seed>=51 & seed<=100){
    s <- 6
    rho1 <- 0.5
    rho2 <- 0.5
    rho3 <- 0.5
    rho4 <- 0.5
    rho5 <- 0.5
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    R1 <- cbind(M1,M2,M3,M4,M5)
    R2 <- cbind(M2,M1,M2,M3,M4)
    R3 <- cbind(M3,M2,M1,M2,M3)
    R4 <- cbind(M4,M3,M2,M1,M2)
    R5 <- cbind(M5,M4,M3,M2,M1)
    sigmaB <- rbind(R1,R2,R3,R4,R5)
    rho <- 0.75
    sigmaE <- toeplitz(rho^seq(0, d - 1))
    lambda.max <- 90000
    name <- "correlated"
  }
  
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  # Gaussian B
  vecB <- mvrnorm(m,rep(0,dq),sigmaB)
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  # Below is for Xi=Zi with intercept (high d)
  Z <- cbind(1,matrix(rnorm(N*(q-1)),nrow = N,ncol = q-1))
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- cbind(1,matrix(rnorm(N*(p-1)),nrow = N,ncol = p-1))
  beta <- matrix(c(-2:2,rep(0,15),2,1,1,-2,-1,rep(0,15),-1,2,1,2,-2,rep(0,15),-2,2,1,1,2,rep(0,15),1,2,2,-2,-1,rep(0,15)),nrow = p,ncol = d)
  
  Y <- X %*% beta + as.matrix(Z_long) %*% B + epsilon
  
  lambda.min <- lambda.max/100
  length.out <- 100
  lambdas <- exp(seq(log(lambda.max),log(lambda.min),length.out=length.out))
  
  tau.max <- 10
  tau.min <- tau.max/100
  length.out <- 100
  taus <- seq(tau.max,tau.min,length.out=length.out)
  
  t1 <- Sys.time()
  result.sim <- multi.MOMENT(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/sensitivity/our_marginal_high_1_4/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}


########## MOMENT(Marginal),q=8,s=4 ###########
sim.fun <- function(seed){
  d <- 5
  q <- 8
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 4
    rho1 <- 0.5
    rho2 <- 0.6
    rho3 <- 0.7
    rho4 <- 0.8
    rho5 <- 0.9
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    sigmaB <- bdiag(M1,M2,M3,M4,M5)
    rho <- 0.75
    sigmaE <- diag(c(1,0.9,0.8,0.7,0.6))
    lambda.max <- 60000
    name <- "independent"
  }
  
  if (seed>=51 & seed<=100){
    s <- 4
    rho1 <- 0.5
    rho2 <- 0.5
    rho3 <- 0.5
    rho4 <- 0.5
    rho5 <- 0.5
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    R1 <- cbind(M1,M2,M3,M4,M5)
    R2 <- cbind(M2,M1,M2,M3,M4)
    R3 <- cbind(M3,M2,M1,M2,M3)
    R4 <- cbind(M4,M3,M2,M1,M2)
    R5 <- cbind(M5,M4,M3,M2,M1)
    sigmaB <- rbind(R1,R2,R3,R4,R5)
    rho <- 0.75
    sigmaE <- toeplitz(rho^seq(0, d - 1))
    lambda.max <- 60000
    name <- "correlated"
  }
  
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  # Gaussian B
  vecB <- mvrnorm(m,rep(0,dq),sigmaB)
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  # Below is for Xi=Zi with intercept (high d)
  Z <- cbind(1,matrix(rnorm(N*(q-1)),nrow = N,ncol = q-1))
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- cbind(1,matrix(rnorm(N*(p-1)),nrow = N,ncol = p-1))
  beta <- matrix(c(-2:2,rep(0,15),2,1,1,-2,-1,rep(0,15),-1,2,1,2,-2,rep(0,15),-2,2,1,1,2,rep(0,15),1,2,2,-2,-1,rep(0,15)),nrow = p,ncol = d)
  
  Y <- X %*% beta + as.matrix(Z_long) %*% B + epsilon
  
  lambda.min <- lambda.max/100
  length.out <- 100
  lambdas <- exp(seq(log(lambda.max),log(lambda.min),length.out=length.out))
  
  tau.max <- 10
  tau.min <- tau.max/100
  length.out <- 100
  taus <- seq(tau.max,tau.min,length.out=length.out)
  
  t1 <- Sys.time()
  result.sim <- multi.MOMENT(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/sensitivity/our_marginal_low_1_4/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}


########## Thompson,q=15,s=6 ###########
sim.fun <- function(seed){
  d <- 5
  q <- 15
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 6
    rho1 <- 0.5
    rho2 <- 0.6
    rho3 <- 0.7
    rho4 <- 0.8
    rho5 <- 0.9
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    sigmaB <- bdiag(M1,M2,M3,M4,M5)
    rho <- 0.75
    sigmaE <- diag(c(1,0.9,0.8,0.7,0.6))
    name <- "independent"
  }
  
  if (seed>=51 & seed<=100){
    s <- 6
    rho1 <- 0.5
    rho2 <- 0.5
    rho3 <- 0.5
    rho4 <- 0.5
    rho5 <- 0.5
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    R1 <- cbind(M1,M2,M3,M4,M5)
    R2 <- cbind(M2,M1,M2,M3,M4)
    R3 <- cbind(M3,M2,M1,M2,M3)
    R4 <- cbind(M4,M3,M2,M1,M2)
    R5 <- cbind(M5,M4,M3,M2,M1)
    sigmaB <- rbind(R1,R2,R3,R4,R5)
    rho <- 0.75
    sigmaE <- toeplitz(rho^seq(0, d - 1))
    name <- "correlated"
  }
  
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  # Gaussian B
  vecB <- mvrnorm(m,rep(0,dq),sigmaB)
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  # Below is for Xi=Zi with intercept (high d)
  Z <- cbind(1,matrix(rnorm(N*(q-1)),nrow = N,ncol = q-1))
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  Z_no_intercept <- Z[,-1]
  X <- cbind(1,matrix(rnorm(N*(p-1)),nrow = N,ncol = p-1))
  U <- cbind(X,Z_no_intercept)[,-1]
  beta <- matrix(c(-2:2,rep(0,15),2,1,1,-2,-1,rep(0,15),-1,2,1,2,-2,rep(0,15),-2,2,1,1,2,rep(0,15),1,2,2,-2,-1,rep(0,15)),nrow = p,ncol = d)
  
  Y <- X %*% beta + as.matrix(Z_long) %*% B + epsilon
  
  lambda_seq <- seq(1,100,length=100)
  
  t1 <- Sys.time()
  result.sim <- multi.glmmsel.sensitivity(Y,U,id,lambda_seq,p,q,sigmaB,beta)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/sensitivity/Thompson_high/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}


########## Thompson,q=8,s=4 ###########
sim.fun <- function(seed){
  d <- 5
  q <- 8
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 4
    rho1 <- 0.5
    rho2 <- 0.6
    rho3 <- 0.7
    rho4 <- 0.8
    rho5 <- 0.9
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    sigmaB <- bdiag(M1,M2,M3,M4,M5)
    rho <- 0.75
    sigmaE <- diag(c(1,0.9,0.8,0.7,0.6))
    name <- "independent"
  }
  
  if (seed>=51 & seed<=100){
    s <- 4
    rho1 <- 0.5
    rho2 <- 0.5
    rho3 <- 0.5
    rho4 <- 0.5
    rho5 <- 0.5
    M1 <- bdiag(toeplitz(rho1^seq(0, s-1)),matrix(0,q-s,q-s))
    M2 <- bdiag(toeplitz(rho2^seq(0, s-1)),matrix(0,q-s,q-s))
    M3 <- bdiag(toeplitz(rho3^seq(0, s-1)),matrix(0,q-s,q-s))
    M4 <- bdiag(toeplitz(rho4^seq(0, s-1)),matrix(0,q-s,q-s))
    M5 <- bdiag(toeplitz(rho5^seq(0, s-1)),matrix(0,q-s,q-s))
    R1 <- cbind(M1,M2,M3,M4,M5)
    R2 <- cbind(M2,M1,M2,M3,M4)
    R3 <- cbind(M3,M2,M1,M2,M3)
    R4 <- cbind(M4,M3,M2,M1,M2)
    R5 <- cbind(M5,M4,M3,M2,M1)
    sigmaB <- rbind(R1,R2,R3,R4,R5)
    rho <- 0.75
    sigmaE <- toeplitz(rho^seq(0, d - 1))
    name <- "correlated"
  }
  
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  # Gaussian B
  vecB <- mvrnorm(m,rep(0,dq),sigmaB)
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  # Below is for Xi=Zi with intercept (high d)
  Z <- cbind(1,matrix(rnorm(N*(q-1)),nrow = N,ncol = q-1))
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  Z_no_intercept <- Z[,-1]
  X <- cbind(1,matrix(rnorm(N*(p-1)),nrow = N,ncol = p-1))
  U <- cbind(X,Z_no_intercept)[,-1]
  beta <- matrix(c(-2:2,rep(0,15),2,1,1,-2,-1,rep(0,15),-1,2,1,2,-2,rep(0,15),-2,2,1,1,2,rep(0,15),1,2,2,-2,-1,rep(0,15)),nrow = p,ncol = d)
  
  Y <- X %*% beta + as.matrix(Z_long) %*% B + epsilon
  
  lambda_seq <- seq(1,100,length=100)
  
  t1 <- Sys.time()
  result.sim <- multi.glmmsel.sensitivity(Y,U,id,lambda_seq,p,q,sigmaB,beta)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/sensitivity/Thompson_low/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

library(parallel)
seed = 1:100
ncores = detectCores()
#ncores
cl = makeCluster(ncores)
result_list = mclapply(seed, sim.fun, mc.cores = 12)
stopCluster(cl)



