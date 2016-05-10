# Reading data
age18_29 <- c(M=208, C=45, FP=46, KD=35, MP=110,
  S=189, V=34, SD=53, Others=88)

age30_49 <- c(M=403, C=58, FP=74, KD=42, MP=146,
              S=413, V=127, SD=93, Others=57)

age50_64 <- c(M=370, C=51, FP=60, KD=47, MP=67,
              S=401, V=59, SD=61, Others=15)

age65 <- c(M=383, C=89, FP=86, KD=65, MP=45,
              S=567, V=74, SD=79, Others=17)

df <- cbind(age18_29,age30_49,age50_64,age65) # Put everything in one matrix
prior <- c(alpha1=30, alpha2=6, alpha3=7, alpha4=6, alpha5=7, # Reading prior
alpha6=30, alpha7=6, alpha8=6, alpha9=2)  
df <- as.data.frame(cbind(df, prior)) # Make a data.frame of all data and given priors



#install.packages("gtools") # Needed for simulating from Drichlet.
library(gtools)

mylist <- list(age18_29=c(),age30_49 =c(),age50_64=c(),age65=c()) # Create a list to fill with data for every age_group


# Creates a list with 1000 simulations of probabbilites of the partys for each age group.
for(i in 1:4){
  set.seed(311015)
  temp <- as.data.frame(rdirichlet(1000, df[,i] + df$prior))
  colnames(temp) <- c("M", "C", "FP", "KD", "MP", "S", "V", "SD", "Others")
  mylist[[i]] <- list(Prob=temp)
}
library(reshape2)
library(ggplot2)

# Setting up the required data.frame for ggplot.
gg_groupedBar <- data.frame(NA)
for(i in 1:4){
  for(j in 1:1000){
gg_groupedBar[i,j] <- c(colMeans(mylist[[i]]$Prob), agegroup=i)[j]
}}
colnames(gg_groupedBar) <- c("M", "C", "FP", "KD", "MP", "S", "V", "SD", "Others", "agegroup")
gg_groupedBar <- melt(gg_groupedBar,id.vars = "agegroup")

# Plottings
ggplot(data=gg_groupedBar[1:36,], aes(x=variable, y=value, fill=factor(agegroup))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  labs(title="Predicted party votes", x="", y="Votes in %") +
  scale_fill_manual(values=c("royalblue","seagreen","darkorange","indianred"),
                    name="Agegroups",
                    labels=c("18-29", "30-49", "50-64", "65<")) +
  theme_bw() + theme(legend.position="bottom")
# Redgreen vs Alliansen. Caluclating amount of simulations in favor of redgreen wins.
redgreen_alliansen <- c()
for(i in 1:4){
  for(j in 1:1000){
    redgreen_alliansen[j] <- sum(mylist[[i]]$Prob[j,c(1:4)]) < sum(mylist[[i]]$Prob[j,c(5:7)])
  }
  mylist[[i]]$redgreen_alliansen <- mean(redgreen_alliansen)
}

# Setting up the ggplot required dataframe.
gg_redgreen_alliansen <- data.frame(agegroup=1:4, variable=rep("redgreen_alliansen", 4),
                                    value=c(mylist$age18_29$redgreen_alliansen, mylist$age30_49$redgreen_alliansen,
                                      mylist$age50_64$redgreen_alliansen, mylist$age65$redgreen_alliansen))

# Plotting
ggplot(data=gg_redgreen_alliansen, aes(x=variable, y=value, fill=factor(agegroup))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  labs(title="Predicted wins for left parties", x="", y="Wins in %") +
  scale_fill_manual(values=c("royalblue","seagreen","darkorange","indianred"),
                    name="Agegroups",
                    labels=c("18-29", "30-49", "50-64", "65<")) +
  theme_bw() + theme(legend.position="bottom")
# Simulate amount of voters for each age group.
set.seed(311015)
numVoter <- as.data.frame(t(rmultinom(1000,6300000, c(0.2, 0.3, 0.3, 0.2))))
colnames(numVoter) <- c("age18_29", "age30_49", "age50_64", "age65")

# Calculate amount of voters for each party and age group using the
# simulated number of voters from numVoter.
Prob_numVoter <- as.data.frame(matrix(NA, nrow = 1000, ncol=9)) 
colnames(Prob_numVoter) <- c("M", "C", "FP", "KD", "MP", "S", "V", "SD", "Others")
for(i in 1:4){
  for(j in 1:1000) {
    Prob_numVoter[j,] <- mylist[[i]][[1]][j,]*numVoter[j,i]    
  }
  mylist[[i]]$Prob_numVoter <- Prob_numVoter
}


total_win <- c()
# Calculate how many times rödgröna wins over allisansen.
for(j in 1:1000){
  temp <- mylist[[1]]$Prob_numVoter[j,] + mylist[[2]]$Prob_numVoter[j,] +
  mylist[[3]]$Prob_numVoter[j,] + mylist[[4]]$Prob_numVoter[j,]
total_win[j] <- sum(temp[1:4]) < sum(temp[5:7])
}

mylist$Total <- mean(total_win)
library(gtools)
library(ggplot2)
library(plyr)
library(dplyr)
JapanTemp <- read.delim("C:/Users/Gustav/Documents/Machine-Learning/Lab 6/JapanTemp.dat", sep="", header = TRUE)
ClassicLM <- lm(temp ~ time+ I(time^2), data=JapanTemp)
summary(ClassicLM)
library(geoR)
library(mvtnorm)
regLine <- data.frame(matrix(vector(), 365, 100)) 
set.seed(311015)
for(i in 1:100){
  sigma0 <- rinvchisq(1, df = 1, scale = 122.82)
  priorCoef <- rmvnorm(n=1, mean = c(11.58,58,-50), sigma = diag(x=sigma0/5, 3, 3))
  regLine[,i] <- priorCoef[1] + priorCoef[2] * JapanTemp$time + priorCoef[3] * JapanTemp$time^2
}
require(reshape2)
regLine_m <- melt(regLine)
regLine_m$x <- rep(1:365, 100)
ggplot()+geom_line(data=regLine_m,aes(x=x,y=value,group=variable),col="royalblue") +  geom_line(data=data.frame(y=ClassicLM$fitted.values), aes(y=y,x=1:365),
            col="indianred", size=1.25) + theme_bw() +
  geom_point(data=JapanTemp,aes(x=1:365, y=temp),col="indianred", size=3)+
  ggtitle("Regression curves for the prior distribution, \n red dots are observed values") + xlab("Day of year") + ylab("Temperature")
regLine <- data.frame(matrix(vector(), 365, 100)) 
set.seed(311015)
for(i in 1:100){
  sigma0 <- rinvchisq(1, df = 10, scale = 122.82)
  priorCoef <- rmvnorm(n=1, mean = c(11.58,58,-50), sigma = diag(x=sigma0/30, 3, 3))
  regLine[,i] <- priorCoef[1] + priorCoef[2] * JapanTemp$time + priorCoef[3] * JapanTemp$time^2
}
regLine_m <- melt(regLine)
regLine_m$x <- rep(1:365, 100)
ggplot()+geom_line(data=regLine_m,aes(x=x,y=value,group=variable),col="royalblue") + geom_line(data=data.frame(y=ClassicLM$fitted.values), aes(y=y,x=1:365),
            col="indianred", size=1.25) + theme_bw() +
  geom_point(data=JapanTemp,aes(x=1:365, y=temp),col="indianred", size=3)+
  ggtitle("Regression curves for the prior distribution, \n red dots are observed values") + xlab("Day of year") + ylab("Temperature")
X <- as.matrix(data.frame(int=rep(1, 365), x=JapanTemp$time,x2=JapanTemp$time^2))
omega0 <- diag(x=30, 3,3)
beta0 <- (as.matrix(c(11.58,58,-50)))
v0 <- 10
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
PosteriorLine_m <- melt(PosteriorLine)
PosteriorLine_m$x <- rep(1:365, 250)
Posterior <- colMeans(PosteriorCoef)
PosteriorLine1 <- data.frame(y=Posterior[1]+Posterior[2]*JapanTemp$time+
                               Posterior[3]*JapanTemp$time^2, x=JapanTemp$time)
ggplot()+geom_line(data=PosteriorLine_m,aes(x=x,y=value,group=variable),col="royalblue")  +
  theme_bw()+geom_point(data=JapanTemp,aes(x=1:365, y=temp),col="indianred", size=3) +
  ggtitle("Regression curves for the posterior distribution, \n red dots are observed values") + xlab("Day of year") + ylab("Temperature") +
  geom_line(data=PosteriorLine1,aes(x=1:365,y=y),col="black")
DayMax <- data.frame(x=-(PosteriorCoef[,2] / (2*PosteriorCoef[,3])))
ggplot(DayMax, aes(x)) + geom_histogram(binwidth=0.002) + theme_bw() + 
  scale_x_continuous(breaks=c(0.56,0.565,0.57,0.575,0.58,0.585))
## NA
