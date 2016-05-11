
#### Assignment 1 ####
rainfall <- read.delim("C:/Users/Gustav/Documents/Bayesian Learning/Lab3/rainfall.dat",
                       sep="", header = TRUE)
library(ggplot2)
library(gridExtra)

colMeans(rainfall)
## a)
# priors and others
mu0 <- 1
kappa0 <- 1
v0 <- 1
sigma0 <- 1
n <- nrow(rainfall)
ybar <- colMeans(rainfall)
s2 <- var(rainfall[,1])

# Posteriors
muN <- (kappa0 / (kappa0 + n)) * mu0 + (n / (kappa0 + n)) * ybar 
kappaN <- kappa0 + n
vN <- v0 +n 
vNsigmaN <- v0*sigma0 + (n-1)*s2 + (kappa0*n / (kappa0 + n)) * (ybar - mu0)^2
# What is s2? Variance of data?
sigmaN <- vNsigmaN / vN

# Simulations - Gibbs Sampling
sims <- data.frame(mu=0, sigma2=0)

muN <- (kappa0 / (kappa0 + n)) * mu0 + (n / (kappa0 + n)) * ybar 
vNsigmaN <- v0*sigma0 + (n-1)*s2 + (kappa0*n / (kappa0 + n)) * (ybar - muN)^2
sigmaN <- vNsigmaN / vN
X <- rchisq(1, vN)
sigma2 <- (vN * sigmaN / X)
mu <- rnorm(1, muN, sqrt(sigma2/kappaN))  
sims[1,1] <-  mu
sims[1,2] <-  sigma2 
for (i in 2:1000){
  # Byter ut mu0 mot [i-1,1] i sims
  # Byter ut sigma0 mot [i-1,2] i sims
  muN <- (kappa0 / (kappa0 + n)) * sims[i-1,1] + (n / (kappa0 + n)) * ybar 
  # Byter ut mu0 mot muN
  vNsigmaN <- v0*sims[i-1,2] + (n-1)*s2 + (kappa0*n / (kappa0 + n)) * (ybar - muN)^2
  sigmaN <- vNsigmaN / vN
  
  # Below, should it be vN or n-1? (See lecture 3)
  X <- rchisq(1, n-1)
  sigma2 <- (vN * sigmaN / X)
  mu <- rnorm(1, muN, sqrt(sigma2/kappaN))  
  sims[i,1] <-  mu
  sims[i,2] <-  sigma2 
}
#trace plot - mu, remove burn-in!
tr_w_burn <- ggplot(sims, aes(x=1:nrow(sims), y=mu)) + geom_line() + theme_bw() + xlab("1:1000") +
  ylab("Mu") + ggtitle("Gibbs Sampling - Trace Plot - Mu")
ggplot(sims,aes(mu))+geom_histogram(aes(y = ..density..),alpha=0.65,
  fill="royalblue", binwidth=0.2)+geom_density(size=1.05, binwidth=0.2) +
  theme_bw() + ggtitle("Gibbs Sampling - Density Plot + Histogram")
# trace plot - sigma2
tr_w_burn2 <- ggplot(sims, aes(x=1:nrow(sims), y=sigma2)) + geom_line() +
  theme_bw() +ggtitle("Gibbs Sampling - Trace Plot - Sigma^2")+ xlab("1:1000")
grid.arrange(tr_w_burn, tr_w_burn2, ncol=2)

# Burn-in is first 100 ? If so:
ggplot(sims[101:1000,], aes(x=101:nrow(sims), y=mu)) + geom_line() + theme_bw() + xlab("1:1000") +
  ylab("Mu") + ggtitle("Gibbs Sampling - Trace Plot - Mu")
ggplot(sims[101:1000,],aes(mu))+geom_histogram(aes(y = ..density..),alpha=0.65,
                                    fill="royalblue", binwidth=0.2)+geom_density(size=1.05, binwidth=0.2) +
  theme_bw() + ggtitle("Gibbs Sampling - Density Plot + Histogram")
# trace plot - sigma2
ggplot(sims[101:1000,], aes(x=101:nrow(sims), y=sigma2)) + geom_line() + theme_bw()
# Look at efficieny in terms of auto-correlation
acf_m <- acf(sims[101:1000,1], lag.max = 30, type = c("correlation"), plot=FALSE)
acf_s <- acf(sims[101:1000,2], lag.max = 30, type = c("correlation"), plot=FALSE)
acf_m <- data.frame(ACF=as.numeric(acf_m$acf), Lag=0:30)
acf_s <- data.frame(ACF=as.numeric(acf_s$acf), Lag=0:30)
ggplot(acf_m, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw()
ggplot(acf_s, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw()

## b)
## Se r-script - upg1.b)

## c)
colMeans(sims[101:1000,])
a1 <- data.frame(x=rnorm(1000, 32.27564, sqrt(1546.53868)))
b1 <- data.frame(y=mixDensMean, x=xGrid)
ggplot(rainfall, aes(X136)) + geom_histogram(aes(y = ..density..),alpha=0.9,
  fill="black") + theme_bw() +
  geom_density(data=a1, aes(x),col="royalblue", size=1.05) +
  geom_line(data=b1, aes(x=x, y=y), col="red", size=1.05) + 
  ggtitle("Density for Gibbs Sampling and Mixture of normals\n Blue = 1.a) - Red = 1.b)") + xlab("")

#### Assignment 2 ####

# b-c)
library(msm)
?rtnorm

tau <- 10
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara)
nPara <- dim(X)[2]

rtnorm(16, mean=as.vector(rep(0,nPara)), sqrt(10))

library(mvtnorm)
initVal <- as.vector(rep(0,dim(X)[2]))

set.seed(1234)
initval <- as.vector(rmvnorm(n=1, mean = initVal, sigma = diag(x=tau^2, 16, 16)))

mean_p <- t(as.matrix(as.vector(rep(0,dim(X)[2])), ncol=1))
sigma_p <- diag(x=tau^2, 16, 16)
sigma_p2 <- ( as.matrix(diag(sigma_p)))

#B_n <- solve(t(X)%*%X + sigma_p2%*%mean_p) %*% t(X)%*%y
#mean_p <- t(as.matrix(B_n))
#sigma_p <- solve(t(X)%*%X + sigma_p)
#sigma_p2 <- ( as.matrix(diag(sigma_p)))

rmvnorm(1, mean_p, sigma_p)

emptyB <- data.frame(matrix(vector(), 101, 16))
u <- data.frame(matrix(vector(), 4601, 101))
set.seed(311015)
u[,1] <- rtnorm(4601, X%*%t(mean_p), sd = rep(1, 16)) 

for (i in 1:101){
  B_n <- solve(t(X)%*%X + sigma_p2%*%mean_p) %*% t(X)%*%u[,i]
  mean_p <- t(as.matrix(B_n))
  sigma_p <- solve(t(X)%*%X + sigma_p)
  sigma_p2 <- ( as.matrix(diag(sigma_p)))
  
  emptyB[i,] <-  rmvnorm(1, mean_p, sigma_p)
  newB <- t(matrix(as.numeric(emptyB[i,])))
  
  for(j in 1:4601){
    if(Data$spam[j] == 1){
      u[j,i+1] <- rtnorm(1, t(X[i,] %*% newB),sd = rep(1, 16), lower=0, upper=Inf)
    }else{
      u[j,i+1] <- rtnorm(1, t(X[i,] %*% newB),sd = rep(1, 16), lower=-Inf, upper=0)
    }
  }
}


Betas <- matrix(as.numeric(emptyB[20,]))

y_fit <- (X) %*% Betas
for(i in 1:4601){
  if(y_fit[i] > 0){
    y_fit[i] = 1
  }else{
    y_fit[i] = 0
  }
}
table(y,y_fit)

betasOptim <- matrix(as.numeric(OptimResults$par))
y_fitOptim <- (X) %*% betasOptim
for(i in 1:4601){
  if(y_fitOptim[i] > 0){
    y_fitOptim[i] = 1
  }else{
    y_fitOptim[i] = 0
  }
}
table(y,y_fitOptim)



