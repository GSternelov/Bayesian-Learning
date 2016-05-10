# Lab 2 - Bayeisian analysis

### Assignment 1 ###
age18_29 <- c(M=208, C=45, FP=46, KD=35, MP=110,
              S=189, V=34, SD=53, Others=88)
age30_49 <- c(M=403, C=58, FP=74, KD=42, MP=146,
              S=413, V=127, SD=93, Others=57)
age50_64 <- c(M=370, C=51, FP=60, KD=47, MP=67,
              S=401, V=59, SD=61, Others=15)
age65 <- c(M=383, C=89, FP=86, KD=65, MP=45,
           S=567, V=74, SD=79, Others=17)
Votes <- data.frame(cbind(age18_29,age30_49,age50_64,age65))
Votes$Total <- rowSums(Votes)

priorDir <- c(30,6,7,6,7,30,6,6,2)
Posteriors <- data.frame(matrix(vector(), 9, 4))
for(i in 1:4){
for(j in 1:9){
    Posteriors[j,i] <- priorDir[j] + Votes[j,i] 
  }
}

library(gtools)
library(ggplot2)
library(plyr)
library(dplyr)

set.seed(311015)
age1 <- data.frame(y=rdirichlet(1000, c(Posteriors[,1])),age=1 )
age2 <- data.frame(y=rdirichlet(1000, c(Posteriors[,2])),age=2 )
age3 <- data.frame(y=rdirichlet(1000, c(Posteriors[,3])),age=3 )
age4 <- data.frame(y=rdirichlet(1000, c(Posteriors[,4])),age=4 )
ageframe <- data.frame(rbind(age1,age2,age3,age4))

ggplot(ageframe, aes(x=age, y=y.1/1000)) + geom_bar(stat="identity")
ggplot(ageframe, aes(x=age, y=y.2/1000)) + geom_bar(stat="identity", aes(fill=age))+
  theme_bw()

# c)
ageframe$Winner <- 0
for(i in 1:4000){
  if( sum(ageframe[i, 5:7]) > sum(ageframe[i, 1:4]) ){
    ageframe$Winner[i] = 1    
  }else{
    ageframe$Winner[i] = 0
  }
}
mean(ageframe[,11] > 0)
mean(ageframe[1:1000,11] > 0)
mean(ageframe[1001:2000,11] > 0)
mean(ageframe[2001:3000,11] > 0)
mean(ageframe[3001:4000,11] > 0)

ggplot(ageframe, aes(x=age, y=Winner/1000, group=age, fill=age))+
  geom_bar(stat="identity") + theme_bw()

# d)
## already have draws for theta, se ageframe
set.seed(311015)
dframe <- data.frame(t(x=rmultinom(n=1000, 6300000, c(0.2,0.3,0.2,0.3))))

  Posteriors2 <- data.frame(matrix(vector(), 1, 9))
  h <- 0
  for(i in 1:4){
  for(j in 1:1000){
      h <- h+1
      Posteriors2[h,1:9] <- dframe[j,i] * ageframe[h, 1:9]  
    }
  }

Posteriors2$age <- ageframe$age
Winner <- 0
for(i in 1:1000){
  index <- c(i, i+1000, i+2000,i+3000)
  if( sum(Posteriors2[index, 5:7]) > sum(Posteriors2[index, 1:4]) ){
    Winner[i] <- 1    
  }else{
    Winner[i] <- 0
  }
}
mean(Winner > 0)
mean(Posteriors2[1:1000,11] > 0)
mean(Posteriors2[1001:2000,11] > 0)
mean(Posteriors2[2001:3000,11] > 0)
mean(Posteriors2[3001:4000,11] > 0)

## Assignment 2
JapanTemp <- read.delim("C:/Users/Gustav/Documents/Machine-Learning/Lab 6/JapanTemp.dat", 
                        sep="", header = TRUE)
summary(lm(temp ~ time+ I(time^2), data=JapanTemp))
ClassicLM <- lm(temp ~ time+ I(time^2), data=JapanTemp)

summary(lm(temp ~ time+ I(time^2)+ I(time^3)+ I(time^4)+ I(time^5)+I(time^6)+I(time^7),
   data=JapanTemp))

sum((11.58 - JapanTemp$time)^2) / 365

# b)
library(geoR)
library(mvtnorm)
regLine <- data.frame(matrix(vector(), 365, 100)) 
set.seed(311015)
for(i in 1:100){
  sigma0 <- rinvchisq(1, df = 5, scale = 122.82)
  priorCoef <- rmvnorm(n=1, mean = c(11.58,58,-50), sigma = diag(x=sigma0/30, 3, 3))
  regLine[,i] <- priorCoef[1] + priorCoef[2] * JapanTemp$time + priorCoef[3] * JapanTemp$time^2
}
#plot(regLine[,1], ylim=c(5, 30), type="l")
#for(i in 1:99){
#  points(regLine[,i+1], type="l")
#}
#points(ClassicLM$fitted.values, type="l", col="red", lwd=3)

require(reshape2)
regLine_m <- melt(regLine)
regLine_m$x <- rep(1:365, 100)

ggplot()+geom_line(data=regLine_m,aes(x=x,y=value,group=variable)) +
  geom_line(data=data.frame(y=ClassicLM$fitted.values), aes(y=y,x=1:365),
            col="indianred", size=1.25) + theme_bw() +
  geom_point(data=JapanTemp,aes(x=1:365, y=temp),col="red", size=3)
  
  
# c)
X <- as.matrix(data.frame(int=rep(1, 365), x=JapanTemp$time, x2=JapanTemp$time^2))
omega0 <- diag(x=30, 3,3)
beta0 <- (as.matrix(c(11.58,58,-50)))
#beta0 <- (as.matrix(c(0,0,0)))
v0 <- 5
s0 <- 122.82

betaHat <- solve(t(X)%*%X) %*%
  t(X) %*% JapanTemp[,2]
omegaNew <- t(X)%*%X + omega0
betaNew <- (solve(t(X)%*%X + omega0)) %*% 
  ((t(X)%*%X%*%betaHat)+(omega0%*%beta0)) 
vNew <- v0 + nrow(JapanTemp)
vNew_sNew <- v0*s0 + t(JapanTemp[,2])%*%JapanTemp[,2] +
  t(beta0)%*%omega0%*%(beta0)- t(betaNew)%*%omegaNew%*%(betaNew)
sNew <- vNew_sNew/vNew

PosteriorLine <- data.frame(matrix(vector(), 365, 250))
PosteriorCoef <- data.frame(matrix(vector(), 250, 3)) 
set.seed(311015)
for(i in 1:250){
  sigma0 <- rinvchisq(1, df = vNew, scale = sNew)
  PosteriorCoef[i,] <- rmvnorm(n=1, mean = betaNew, sigma = as.numeric(sigma0) * solve(omegaNew))
  PosteriorLine[,i] <- PosteriorCoef[i,1]+PosteriorCoef[i,2]*JapanTemp$time+PosteriorCoef[i,3]*JapanTemp$time^2
}

colMeans(PosteriorCoef)
PosteriorLine_m <- melt(PosteriorLine)
PosteriorLine_m$x <- rep(1:365, 250)

ggplot()+geom_line(data=PosteriorLine_m,aes(x=x,y=value,group=variable),col="royalblue") +
  theme_bw()+geom_point(data=JapanTemp,aes(x=1:365, y=temp),col="red", size=3) +
  ggtitle("Posterior Distribution") + xlab("Day of year") + ylab("Temperature")


## d)
Posterior <- colMeans(PosteriorCoef)
maxDay <- -(Posterior[2] / (2*Posterior[3]))
DayMax <- data.frame(x=-(PosteriorCoef[,2] / (2*PosteriorCoef[,3])))
ggplot(DayMax, aes(x)) + geom_histogram(binwidth=0.002) + theme_bw() + 
  scale_x_continuous(breaks=c(0.56,0.565,0.57,0.575,0.58,0.585))

PosteriorMaxTemp <- data.frame(y=PosteriorCoef[,1]+PosteriorCoef[,2]*
                                 maxDay+PosteriorCoef[,3]*maxDay^2, x=maxDay)
ggplot(PosteriorMaxTemp, aes(y)) + geom_histogram(binwidth=0.075) + theme_bw() 

PosteriorLine1 <- data.frame(y=Posterior[1]+Posterior[2]*JapanTemp$time+
                               Posterior[3]*JapanTemp$time^2, x=JapanTemp$time)
ggplot(PosteriorLine1, aes(y=y, x=x)) + geom_line() + theme_bw() + 
  geom_line(data=PosteriorMaxTemp, aes(y=y, x=x)) + 
  geom_point(data=JapanTemp,aes(x=time, y=temp), col="red")


maxTemp <- ddply(PosteriorLine_m,.(x), summarize, EstTemp=mean(value))
day209 <- filter(PosteriorLine_m, x == 209)
ggplot(day209, aes(x=value)) + geom_histogram(binwidth=0.04, fill="indianred") +
  theme_bw()
  ggplot(day209, aes(x=1:250, y=value)) + geom_point()
ggplot(maxTemp, aes(x=x, y=EstTemp)) + geom_line(size=1.25, col="royalblue") +
  theme_bw() + geom_point(data=JapanTemp,aes(x=1:365, y=temp), col="red") +
  geom_point(data=day209, aes(y=value, x=x),col="springgreen3")

mday <- mean(day209$value)
mvar <- var(day209$value)
posteriorDay <- data.frame(x=rnorm(1000, mean=mday, sd=sqrt(mvar)))

# Använder varje värde från d) och dag 209, sätter som mean i normalfördelning?

ggplot(maxTemp, aes(x=x, y=EstTemp)) + geom_line(size=1.25, col="royalblue") +
  theme_bw() + geom_point(data=JapanTemp,aes(x=1:365, y=temp), col="red") +
  geom_point(data=day209, aes(y=value, x=x),col="darkorange") + 
  geom_point(data=posteriorDay, aes(y=x, x=209),col="springgreen3") +
  ylim(25,29)



time209 <- filter(JapanTemp, time==round(209/365,4))[,1]
betaTemp <- as.matrix(c(betaNew[1], betaNew[2]*time209, betaNew[3]*time209^2),nrow=3)
vNew_sNewTemp <- v0*s0 + t(JapanTemp[,2])%*%JapanTemp[,2] +
  t(beta0)%*%omega0%*%(beta0)- t(betaTemp)%*%omegaNew%*%(betaTemp)
sNewTemp <- vNew_sNewTemp/vNew

PosteriorLine2 <- data.frame(matrix(vector(), 1, 250)) 
set.seed(311015)
for(i in 1:250){
  sigma0 <- rinvchisq(1, df = vNew, scale = sNewTemp)
  CoefPost2 <- rmvnorm(n=1, mean = betaTemp, sigma = as.numeric(sigma0) * solve(omegaNew))
  PosteriorLine2[,i] <- CoefPost2[1] +CoefPost2[2]*JapanTemp[209,1] +CoefPost2[3]*JapanTemp[209,1]^2
}

PosteriorLine2_m <- melt(PosteriorLine2)
PosteriorLine2_m$x <- rep(209, 250)

ggplot()+geom_line(data=PosteriorLine2_m,aes(x=x,y=value),col="royalblue") +
  theme_bw()+geom_point(data=JapanTemp,aes(x=1:365, y=temp),col="red", size=3) +
  ggtitle("Posterior Distribution") + xlab("Day of year") + ylab("Temperature")

ggplot()+geom_line(data=PosteriorLine2_m,aes(x=1:250,y=value),col="royalblue") +
  theme_bw()+geom_hline(yintercept=JapanTemp[209,2],col="red", size=1.25) +
  ylab("Temperature")

ggplot(PosteriorLine2_m, aes(x=value)) + geom_histogram(binwidth=0.1) +
  geom_vline(xintercept=JapanTemp[209, 2], size=4, col="indianred") +
  theme_bw()

ggplot(maxTemp, aes(x=x, y=EstTemp)) + geom_line(size=1.25, col="royalblue") +
  theme_bw() + geom_point(data=JapanTemp,aes(x=1:365, y=temp), col="red") +
  geom_point(data=posteriorDay, aes(y=x, x=209),col="springgreen3") +
  geom_point(data=PosteriorLine2_m, aes(y=value, x=x),col="darkorange") +
  ylim(25,29)
