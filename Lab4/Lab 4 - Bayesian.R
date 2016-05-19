
# Lab 4
library(mvtnorm)
library(coda)
## Assignment 1
eBay <- read.delim("C:/Users/Gustav/Documents/Bayesian-Learning/Lab4/eBay.dat", sep="", header=TRUE)

# a)
PoiGLM <- glm(nBids ~ ., family = poisson(), data = eBay[,-2])
summary(PoiGLM)
# VerifyID, Sealed, MajBlem, LogBook and MinBidShare are significant. 
# PowerSeller, Minblem and LargNeg insignificant. 

# b)
y <- as.matrix(eBay[,1])
X <- as.matrix(eBay[,-1])
Sigma <- 100 * solve(t(X) %*% X)
mu <- as.vector(rep(0,ncol(X))) # Prior mean vector

LogPostPoisson <- function(betaVect,y,X,mu,Sigma){
  nPara <- length(betaVect)
  
  logLik <- sum( y *  (X %*% betaVect) - exp(X %*% betaVect)  )
  abs(logLik) == Inf
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE)
  return(logLik + logPrior)
}

initVal <- as.vector(rep(0,dim(X)[2])); 
OptimResults<-optim(initVal,LogPostPoisson,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

approxPostStd <- sqrt(diag(-solve(OptimResults$hessian))) # Computing approximate standard deviations.
print('The posterior mode is:')
print(OptimResults$par)
print('The approximate posterior standard deviation is:')
print(approxPostStd)
# c)

logPostFunc <- function(theta, dist, ...){
  dVal <- sum(dist(theta, ...))
  return(dVal)
}
logPostFunc(theta=as.numeric(y), dist=dpois, lambda=exp((X) %*% as.numeric(posMode)), log=TRUE)

sum(dpois(x=as.numeric(y), lambda=exp((X) %*% as.numeric(posMode)), log=TRUE))
# Proposal distribution
posMode <- OptimResults$par
negHessian <- diag(sqrt(diag(-solve(OptimResults$hessian))))


betas <- data.frame(matrix(vector(), 10, 9))
betas[1,] <- posMode
MA <- function(func, c, iter, dist, ...){
  accProb <- data.frame(U=0, alpha=0)
  for(i in 1:iter){
    # First, draw a proposal set of betas
    b_point <- rmvnorm(1, mean=as.numeric(betas[i,]), sigma= c * negHessian)
    # Sample a uniform[0,1] valie
    U <- runif(1, 0, 1)
    accProb[i,1] <- U
    # Evaluates dmvnorm for the draw and the former set of betas
    # And do the same for the target distribution
    fx_1 <- logPostFunc(as.numeric(y), dist, lambda=exp((X) %*% as.numeric(b_point)), log = FALSE)
    fx_2 <- logPostFunc(as.numeric(y), dist, lambda=exp((X) %*% as.numeric(betas[i,])), log = FALSE)
    # Then, the ratio is calculated
    ratio <-exp(log(fx_1)-log(fx_2))
    alpha <- min(c(1, ratio))
    accProb[i,2] <- alpha
    
    if (U <= alpha) {
      betas[i+1,] <- b_point
    }else{
      betas[i+1,] <- betas[i,]
    } 
  }
  res <- list(betas=betas, accProb=accProb)
  return(res)
}

Sims <- MA(iter = 2000, c = 2.4/3, dist=dpois)
mean((Sims$accProb[2])[,1])

test <- data.frame(Sims$betas)

par(mfrow=c(3,3))
for(i in 1:9){
  plot(test[,i], type="l")  
}

for(i in 1:9){
  hist(test[,i])  
}

effectiveSize(test)

