#### Assignment 1 ####
rainfall <- read.delim("C:/Users/Gustav/Documents/Bayesian-Learning/Lab3/rainfall.dat", sep="", header = TRUE)
library(ggplot2)
library(gridExtra)
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
  X <- rchisq(1, n-1)
  sigma2 <- (vN * sigmaN / X)
  mu <- rnorm(1, muN, sqrt(sigma2/kappaN))  
  sims[i,1] <-  mu
  sims[i,2] <-  sigma2 
}
#trace plots - with burn-in!
tr_w_burn <- ggplot(sims, aes(x=1:nrow(sims), y=mu)) + geom_line() + theme_bw() + xlab("1:1000") + ylab("Mu") + ggtitle("Gibbs Sampling - Trace Plot - Mu") + geom_vline(xintercept=100, col="red", size=1.05)
tr_w_burn2 <- ggplot(sims, aes(x=1:nrow(sims), y=sigma2)) + geom_line() +
  theme_bw() +ggtitle("Gibbs Sampling - Trace Plot - Sigma^2")+ xlab("1:1000") + geom_vline(xintercept=100, col="red", size=1.05)
grid.arrange(tr_w_burn, tr_w_burn2, ncol=2)
tr_wh_burn <- ggplot(sims[101:1000,], aes(x=101:nrow(sims), y=mu)) + geom_line() + theme_bw() + xlab("101:1000") + ylab("Mu") + ggtitle("Gibbs Sampling - Trace Plot - Mu\n Without burn-in")
tr_wh_burn2 <- ggplot(sims[101:1000,], aes(x=101:nrow(sims), y=sigma2)) + geom_line() + theme_bw() + ggtitle("Gibbs Sampling - Trace Plot - Sigma^2\n Without burn-in") + xlab("101:1000")
grid.arrange(tr_wh_burn, tr_wh_burn2, ncol=2)
# Looks at efficieny in terms of auto-correlation
acf_m <- acf(sims[101:1000,1], lag.max = 25, type = c("correlation"), plot=FALSE)
acf_s <- acf(sims[101:1000,2], lag.max = 25, type = c("correlation"), plot=FALSE)
acf_m <- data.frame(ACF=as.numeric(acf_m$acf), Lag=0:25)
acf_s <- data.frame(ACF=as.numeric(acf_s$acf), Lag=0:25)
acf_mu <- ggplot(acf_m, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw() + ggtitle("Auto-correlation for chain - Mu")
acf_sigma <- ggplot(acf_s, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw() + ggtitle("Auto-correlation for chain - Sigma^2")
grid.arrange(acf_mu, acf_sigma, ncol=2)
rainfall <- read.delim("C:/Users/Gustav/Documents/Bayesian-Learning/Lab3/rainfall.dat",
                       sep="", header = TRUE)
x <- as.matrix(rainfall['X136'])

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- c(10, 10)
muPrior <- c(32.2681, 32.2681)
tau2Prior <- rep(10,nComp) # Prior std theta
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2


# MCMC options
nIter <- 1000 # Number of Gibbs sampling draws


# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  thetaDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    thetaDraws[j] <- rgamma(1,param[j],1)
  }
  thetaDraws = thetaDraws/sum(thetaDraws) # Diving every column of ThetaDraws by the sum of the elements in that column.
  return(thetaDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
theta <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x,plot = FALSE)$density))

simulations <- data.frame(w_1 = 0, w_2 = 0, mu_1 = 0, mu_2 = 0, sigma1=0, sigma2=0)

for (k in 1:nIter){
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  # Update components probabilities
  w <- rDirichlet(alpha + nAlloc)
simulations[k,1] <- w[1]
simulations[k,2] <- w[2]  

  # Update theta's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    theta[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
simulations[k,3] <- theta[1]
simulations[k,4] <- theta[2]  
simulations$Expected[k] <- sum(simulations$w_1[k] * simulations$mu_1[k] + 
                                 simulations$w_2[k] * simulations$mu_2[k])
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - theta[j])^2))/(nu0[j] + nAlloc[j]))
  }
simulations[k,5] <- sigma2[1]
simulations[k,6] <- sigma2[2]
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- w[j]*dnorm(x[i], mean = theta[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    #hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,theta[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + w[j]*compDens
      #lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    #lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    #legend("topleft", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
          # col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
  }
}
Wh_b <- ggplot(simulations, aes(x=1:nrow(simulations), y=Expected)) + geom_line() + theme_bw() + geom_vline(xintercept=100, col="red", size=1.05) + xlab("1:1000") + ylab("Mu") + ggtitle("Trace plot for chain of values from the\n mixture model")
Wh_o_b <- ggplot(simulations[101:nrow(simulations),], aes(x=101:nrow(simulations), y=Expected)) + geom_line() + theme_bw() + xlab("101:1000") + ylab("Mu") + ggtitle("Trace plot for chain of values from the\n mixture model (without burn-in)")
grid.arrange(Wh_b, Wh_o_b, ncol=2)
# Looks at efficieny in terms of auto-correlation
acf_mix <- acf(simulations[101:nrow(simulations),7], lag.max = 25, type = c("correlation"), plot=FALSE)
acf_mix <- data.frame(ACF=as.numeric(acf_mix$acf), Lag=0:25)
ggplot(acf_mix, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw() + ggtitle("Auto-correlation for values from mixture model")

## trace plots for the other parameters
mu1_b <- ggplot(simulations, aes(x=1:1000, y=mu_1)) + geom_line() + theme_bw() + ggtitle("Trace plot for mu_1")
mu2_b <- ggplot(simulations, aes(x=1:1000, y=mu_2)) + geom_line() + theme_bw()+ ggtitle("Trace plot for mu_2")
w1_b <- ggplot(simulations, aes(x=1:1000, y=w_1)) + geom_line() + theme_bw()+ ggtitle("Trace plot for w_1")
w2_b <- ggplot(simulations, aes(x=1:1000, y=w_2)) + geom_line() + theme_bw()+ ggtitle("Trace plot for w_2")
sigma1_b <- ggplot(simulations, aes(x=1:1000, y=sigma1)) + geom_line() + theme_bw()+ ggtitle("Trace plot for sigma_1")
sigma2_b <- ggplot(simulations, aes(x=1:1000, y=sigma2)) + geom_line() + theme_bw()+ ggtitle("Trace plot for sigma_2")
grid.arrange(mu1_b, mu2_b,w1_b,w2_b,sigma1_b,sigma2_b, ncol=2)
# Autocorrelation plots for the other parameters
acf_mu1 <- acf(simulations[101:1000,3], lag.max = 25, type = c("correlation"), plot=FALSE)
acf_mu2 <- acf(simulations[101:1000,4], lag.max = 25, type = c("correlation"), plot=FALSE)
acf_w1 <- acf(simulations[101:1000,1], lag.max = 25, type = c("correlation"), plot=FALSE)
acf_w2 <- acf(simulations[101:1000,2], lag.max = 25, type = c("correlation"), plot=FALSE)
acf_s1 <- acf(simulations[101:1000,5], lag.max = 25, type = c("correlation"), plot=FALSE)
acf_s2 <- acf(simulations[101:1000,6], lag.max = 25, type = c("correlation"), plot=FALSE)
acf_mu1 <- data.frame(ACF=as.numeric(acf_mu1$acf), Lag=0:25)
acf_mu2 <- data.frame(ACF=as.numeric(acf_mu2$acf), Lag=0:25)
acf_w1 <- data.frame(ACF=as.numeric(acf_w1$acf), Lag=0:25)
acf_w2 <- data.frame(ACF=as.numeric(acf_w2$acf), Lag=0:25)
acf_s1 <- data.frame(ACF=as.numeric(acf_s1$acf), Lag=0:25)
acf_s2 <- data.frame(ACF=as.numeric(acf_s2$acf), Lag=0:25)

acf_mu1 <- ggplot(acf_mu1, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw() + ggtitle("Auto-correlation for chain - Mu_1")
acf_mu2 <- ggplot(acf_mu2, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw()+ ggtitle("Auto-correlation for chain - Mu_2")
acf_w1 <- ggplot(acf_w1, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw()+ ggtitle("Auto-correlation for chain - W_1")
acf_w2 <- ggplot(acf_w2, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw()+ ggtitle("Auto-correlation for chain - W_2")
acf_s1 <- ggplot(acf_s1, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw()+ ggtitle("Auto-correlation for chain - Sigma_1")
acf_s2 <- ggplot(acf_s2, aes(x=Lag, y=ACF))+geom_bar(stat="identity", fill="black")+theme_bw() + ggtitle("Auto-correlation for chain - Sigma_2")
grid.arrange(acf_mu1, acf_mu2,acf_w1,acf_w2,acf_s1,acf_s2,ncol=2)

## c)
a1 <- data.frame(x=rnorm(1000, 32.27564, sqrt(1546.53868)))
b1 <- data.frame(y=mixDensMean, x=xGrid)
ggplot(rainfall, aes(X136)) + geom_histogram(aes(y = ..density..),alpha=0.9,
  fill="black") + theme_bw() +
  geom_density(data=a1, aes(x),col="royalblue", size=1.05) +
  geom_line(data=b1, aes(x=x, y=y), col="red", size=1.05) + 
  ggtitle("Density for normal model and mixture of normals\n Blue = 1.a) - Red = 1.b)") + xlab("")
###########   BEGIN USER INPUTS   ################
Probit <- 1 # If Probit <-0, then logistic model is used.
chooseCov <- c(1:16) # Here we choose which covariates to include in the model
tau <- 10; # Prior scaling factor such that Prior Covariance = (tau^2)*I
###########     END USER INPUT    ################

# install.packages("mvtnorm") # Loading a package that contains the multivariate normal pdf
library("mvtnorm") # This command reads the mvtnorm package into R's memory. NOW we can use dmvnorm function.

# Loading data from file
Data<-read.table("C:/Users/Gustav/Documents/Bayesian-Learning/Lab3/SpamReduced.dat",header=TRUE)  # Spam data from Hastie et al.
y <- as.vector(Data[,1]); # Data from the read.table function is a data frame. Let's convert y and X to vector and matrix.
X <- as.matrix(Data[,2:17]);
covNames <- names(Data)[2:length(names(Data))];
X <- X[,chooseCov]; # Here we pick out the chosen covariates.
covNames <- covNames[chooseCov];
nPara <- dim(X)[2];

# Setting up the prior
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara);

# Defining the functions that returns the log posterior (Logistic and Probit models). Note that the first input argument of

# this function must be the one that we optimize on, i.e. the regression coefficients.

LogPostProbit <- function(betaVect,y,X,mu,Sigma){
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  # MQ change:
  # instead of logLik <- sum(y*log(pnorm(linPred)) + (1-y)*log(1-pnorm(linPred)) ) type in the equivalent and
  # much more numerically stable; 
  logLik <- sum(y*pnorm(linPred, log.p = TRUE) + (1-y)*pnorm(linPred, log.p = TRUE, lower.tail = FALSE))
  # The old expression: logLik2 <- sum(y*log(pnorm(linPred)) + (1-y)*log(1-pnorm(linPred)) )
  abs(logLik) == Inf
  #print('-----------------')
  #print(logLik)
  #print(logLik2)
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  return(logLik + logPrior)
}

# Calling the optimization routine Optim. Note the auxilliary arguments that are passed to the function logPost
# Note how I pass all other arguments of the function logPost (i.e. all arguments except betaVect which is the one that we are trying to optimize over) to the R optimizer.
# The argument control is a list of options to the optimizer. Here I am telling the optimizer to multiply the objective function (i.e. logPost) by -1. This is because
# Optim finds a minimum, and I want to find a maximum. By reversing the sign of logPost I can use Optim for my maximization problem.

# Different starting values. Ideally, any random starting value gives you the same optimum (i.e. optimum is unique)
initVal <- as.vector(rep(0,dim(X)[2])); 
# Or a random starting vector: as.vector(rnorm(dim(X)[2]))
# Set as OLS estimate: as.vector(solve(crossprod(X,X))%*%t(X)%*%y); # Initial values by OLS

if (Probit==1){
  logPost = LogPostProbit;
} else{
  logPost = LogPostLogistic;
}

OptimResults<-optim(initVal,logPost,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

library(msm)
library(mvtnorm)
library(coda)

tau <- 10
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara)
nPara <- dim(X)[2]

mean_p <- t(as.matrix(as.vector(rep(0,dim(X)[2])), ncol=1))
sigma_p <- diag(x=tau^2, 16, 16)
sigma_p2 <- ( as.matrix(diag(sigma_p)))

emptyB <- data.frame(matrix(vector(), 250, 16))
u <- data.frame(matrix(vector(), 4601, 20))
set.seed(311015)
u[,1] <- rtnorm(4601, X%*%t(mean_p), sd = rep(1, 16)) 

for (i in 1:250){
  B_n <- solve(t(X)%*%X + sigma_p2%*%mean_p) %*% t(X)%*%u[,i]
  mean_p <- t(as.matrix(B_n))
  sigma_p <- solve(t(X)%*%X + sigma_p) 
  sigma_p2 <- ( as.matrix(diag(sigma_p))) #same as sigma_p, just modified format
  
  emptyB[i,] <-  rmvnorm(1, mean_p, sigma_p)
  newB <- t(matrix(as.numeric(emptyB[i,])))
  
  for(j in 1:4601){
    if(Data$spam[j] == 1){
      u[j,i+1] <- rtnorm(1, t(X[i,] %*% newB),sd = rep(1, 16), lower=0, upper=Inf)
    }else{
      u[j,i+1] <- rtnorm(1, t(X[i,] %*% newB),sd = rep(1, 16), lower=-Inf, upper=0)
    }
  }
  print(i+1)
}

effectiveSize(emptyB)
effectiveSize(u)

Betas <- matrix(as.numeric(emptyB[250,]))
options(scipen = 999)
Res <- data.frame(covs=covNames,OptimBeta = as.numeric(OptimResults$par), GibbsBeta=as.numeric(emptyB[250,]), OptimStd = sqrt(diag(-solve(OptimResults$hessian))), GibbsStd=sqrt(sigma_p2),row.names = NULL)
Res

par(mfrow=c(4,2))
for(i in 1:16){
  plot(emptyB[,i], type="l")
}
for(i in 9:16){
  plot(emptyB[,i], type="l")
}


y_fitGibbs <- (X) %*% Betas
for(i in 1:4601){
  if(y_fitGibbs[i] > 0){
    y_fitGibbs[i] = 1
  }else{
    y_fitGibbs[i] = 0
  }
}
table(y,y_fitGibbs)

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
## NA
