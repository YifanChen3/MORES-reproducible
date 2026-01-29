######################## multiLMMsel #############################
########## Low Gaussian ###########
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
    lambda.max <- 100000 - 50000
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
    lambda.max <- 100000 - 50000
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
  
  # Below is for independent X and Z (low d)
  Z <- matrix(rnorm(N*q),nrow = N,ncol = q)
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- matrix(rnorm(N*p),nrow = N,ncol = p)
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
  result.sim <- multiLMMsel.cv(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/our_cv_final_low/gaussian/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

######## Low Laplace #########
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
    lambda.max <- 100000 - 50000
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
    lambda.max <- 100000 - 50000
    name <- "correlated"
  }
  
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  #Multivariate Laplace
  eig <- eigen(sigmaB)
  Lambda <- eig$values * (eig$values > 1e-8)
  U <- eig$vectors
  L <- U %*% diag(sqrt(Lambda))
  z <- matrix(rnorm(m*dq),m,dq)
  W <- rexp(m, rate = 1)
  vecB <- (z * sqrt(W)) %*% t(L)
  
  
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  # Below is for independent X and Z (low d)
  Z <- matrix(rnorm(N*q),nrow = N,ncol = q)
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- matrix(rnorm(N*p),nrow = N,ncol = p)
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
  result.sim <- multiLMMsel.cv(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/our_cv_final_low/laplace/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

########## High Gaussian ###########
sim.fun <- function(seed){
  d <- 5
  q <- 20
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 8
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
    lambda.max <- 100000
    name <- "independent"
  }
  
  if (seed>=51 & seed<=100){
    s <- 8
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
    lambda.max <- 100000
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
  X <- Z
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
  result.sim <- multiLMMsel.cv(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/our_cv_final_high/gaussian/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

########## High Laplace ###########
sim.fun <- function(seed){
  d <- 5
  q <- 20
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 8
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
    lambda.max <- 100000
    name <- "independent"
  }
  
  if (seed>=51 & seed<=100){
    s <- 8
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
    lambda.max <- 100000
    name <- "correlated"
  }
  
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  #Multivariate Laplace
  eig <- eigen(sigmaB)
  Lambda <- eig$values * (eig$values > 1e-8)
  U <- eig$vectors
  L <- U %*% diag(sqrt(Lambda))
  z <- matrix(rnorm(m*dq),m,dq)
  W <- rexp(m, rate = 1)
  vecB <- (z * sqrt(W)) %*% t(L)
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  # Below is for Xi=Zi with intercept (high d)
  Z <- cbind(1,matrix(rnorm(N*(q-1)),nrow = N,ncol = q-1))
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- Z
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
  result.sim <- multiLMMsel.cv(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/our_cv_final_high/laplace/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

library(parallel)
seed = 1:100
ncores = detectCores()
#ncores
cl = makeCluster(ncores)
t1 = Sys.time()
#t1
result_list = mclapply(seed, sim.fun, mc.cores = 12)
#t2 = Sys.time()
#t2
#t2-t1
stopCluster(cl)

















######################## multiLMMsel (Marginal) #############################
multi.our <- function(Y,Z,X,lambda.path,tau.path,id,sigmaB=NULL,sigmaE=NULL,beta=NULL,threshold = 0.01,alpha = 1,gamma=2){
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  q <- ncol(Z)/m
  p <- ncol(X)
  N <- sum(nis)
  d <- ncol(Y)
  sigmaB.list <- NULL
  sigmaE.final <- matrix(0,d,d)
  beta.final <- matrix(0,p,d)
  lambda.optimal <- NULL
  tau.optimal <- NULL
  #l.lambda.path <- matrix(0,length(lambda.path),d)
  #l.tau.path <- matrix(0,length(tau.path),d)
  for (di in 1:d) {
    y <- as.matrix(Y[,di])
    cv.result <- multiLMMsel.cv(y,Z,X,lambda.path,tau.path,id = id,sigmaB = NULL,sigmaE = NULL,beta=NULL,threshold = threshold,alpha = alpha,gamma=gamma)
    sigmaB.list[[di]] <- cv.result$SigmaB
    sigmaE.final[di,di] <- cv.result$SigmaE
    beta.final[,di] <- as.vector(cv.result$beta)
    lambda.optimal[di] <- cv.result$lambda.optimal
    tau.optimal[di] <- cv.result$tau.optimal
    #l.lambda.path[,di] <- as.vector(cv.result$l.lambda)
    #l.tau.path[,di] <- as.vector(cv.result$l.tau)
  }
  sigmaB.final <- bdiag(sigmaB.list)
  
  if (is.null(sigmaB)==TRUE | (is.null(sigmaE) == TRUE) | (is.null(beta) == TRUE)){
    return(list(lambda.optimal=lambda.optimal,tau.optimal=tau.optimal,
                sigmaB.final=sigmaB.final,sigmaE.final=sigmaE.final,beta.final=beta.final))
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
                lambda.optimal=lambda.optimal,tau.optimal=tau.optimal,
                sigmaB.final=sigmaB.final,sigmaE.final=sigmaE.final,beta.final=beta.final,
                true.sigmaB=sigmaB,true.sigmaE=sigmaE,true.beta=beta))
  }
}

############# High Gaussian #############
sim.our.independent <- function(seed){
  d <- 5
  q <- 20
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 8
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
    lambda.max <- 100000 
    name <- "independent"
  }

  if (seed>=51 & seed<=100){
    s <- 8
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
    lambda.max <- 100000
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
  
  #Below is for Xi=Zi with intercept (high d)
  Z <- cbind(1,matrix(rnorm(N*(q-1)),nrow = N,ncol = q-1))
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- Z
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
  result.sim <- multi.our(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/our_cv_final_high/gaussian_marginal/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

############# High Laplace #############
sim.our.independent <- function(seed){
  d <- 5
  q <- 20
  dq <- d*q
  p <- 20
  m <- 100
  ni <- 10
  
  if (seed>=1 & seed<=50){
    s <- 8
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
    lambda.max <- 100000 
    name <- "independent"
  }
  
  if (seed>=51 & seed<=100){
    s <- 8
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
    lambda.max <- 100000
    name <- "correlated"
  }
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  #Multivariate Laplace
  eig <- eigen(sigmaB)
  Lambda <- eig$values * (eig$values > 1e-8)
  U <- eig$vectors
  L <- U %*% diag(sqrt(Lambda))
  z <- matrix(rnorm(m*dq),m,dq)
  W <- rexp(m, rate = 1)
  vecB <- (z * sqrt(W)) %*% t(L)
  
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  #Below is for Xi=Zi with intercept (high d)
  Z <- cbind(1,matrix(rnorm(N*(q-1)),nrow = N,ncol = q-1))
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- Z
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
  result.sim <- multi.our(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 0.5,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/our_cv_final_high/Laplace_marginal/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

############# Low Gaussian #############
sim.our.independent <- function(seed){
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
    lambda.max <- 50000 
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
    lambda.max <- 50000
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
  
  # Below is for independent X and Z (low d)
  Z <- matrix(rnorm(N*q),nrow = N,ncol = q)
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- matrix(rnorm(N*p),nrow = N,ncol = p)
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
  result.sim <- multi.our(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/our_cv_final_low/gaussian_marginal/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

############# Low Laplace #############
sim.our.independent <- function(seed){
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
    lambda.max <- 50000 
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
    lambda.max <- 50000
    name <- "correlated"
  }
  
  set.seed(seed*1234+739)
  
  id <- rep(1:m,each=ni)
  nis <- tabulate(match(id, unique(id)))
  m <- length(nis)
  N <- sum(nis)
  epsilon <- mvrnorm(N,rep(0,d),sigmaE)
  
  #Multivariate Laplace
  eig <- eigen(sigmaB)
  Lambda <- eig$values * (eig$values > 1e-8)
  U <- eig$vectors
  L <- U %*% diag(sqrt(Lambda))
  z <- matrix(rnorm(m*dq),m,dq)
  W <- rexp(m, rate = 1)
  vecB <- (z * sqrt(W)) %*% t(L)
  B <- do.call(rbind,lapply(1:nrow(vecB), function(i) {
    matrix(vecB[i, ], nrow = q, ncol = d, byrow = FALSE)
  }))
  
  # Below is for independent X and Z (low d)
  Z <- matrix(rnorm(N*q),nrow = N,ncol = q)
  Z_list <- lapply(split(seq_len(nrow(Z)), id),
                   function(idx) Z[idx, , drop = FALSE])
  Z_long <- bdiag(lapply(Z_list, as.matrix))
  X <- matrix(rnorm(N*p),nrow = N,ncol = p)
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
  result.sim <- multi.our(Y,Z,X,lambdas,taus,id = id,sigmaB = sigmaB,sigmaE = sigmaE,beta=beta,threshold = 0.01,alpha = 1,gamma=2)
  t2 <- Sys.time()
  result.sim$time <- as.numeric(difftime(t2, t1, units = "secs"))
  
  file.name <- paste("~/Documents/project/Val_Sel_Multi_Mix_Effect_Model/sim_data/our_cv_final_low/laplace_marginal/",name,"/",seed,".RData",sep = "")
  save(result.sim,file = file.name)
}

library(parallel)
seed = 1:100
ncores = detectCores()
cl = makeCluster(ncores)
result_list = mclapply(seed, sim.our.independent, mc.cores = 12)
stopCluster(cl)










