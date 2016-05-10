
# Lab 1 - Bayesian Learning
 # 1
# a)
set.seed(12345)
beta10 <- rbeta(10, 16, 8)
beta100 <- rbeta(100, 16, 8)
beta1000 <- rbeta(1000, 16, 8)
sd(beta10); sd(beta100) ; sd(beta1000)
sqrt((16*8)/((16+8)^2*(16+8+1)))
mean(beta10); mean(beta100) ; mean(beta1000)
16/(16+8)

# b)
length(subset(beta10, beta10 < 0.4)) / 10
length(subset(beta100, beta100 < 0.4)) / 100
length(subset(beta1000, beta1000 < 0.4)) / 1000
pbeta(0.4, 16, 8)

# c)
log_odds <- log(beta1000 / (1-beta1000))
hist(log_odds)
density(log_odds)


# 2
 # b)
y <- c(14,25,45,25,30,33,19,50,34,67)
sigma2 <- 0
set.seed(12345)
for(i in 1:10000){
  chisq <- rchisq(1, 10)
  sigma2[i] <- sum((log(y)-3.5)^2) / chisq
}
hist(sigma2, breaks=100)
mean(sigma2)
var(sigma2)
# theoretical mean
(sum((log(y)-3.5)^2))/(10-2)
# theoretical variance
(2*10^2 * ((sum((log(y)-3.5)^2))/10)^2 ) / (8^2 * 6)


 # c)
G <- 2*dnorm((sqrt(sigma2)/sqrt(2))-1)
hist(G, breaks=75, col="blue")

# 3
?besselI
radians <- c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)
# a)
kappa <- seq(from=0, to=5, by=0.001)
distKappa <- (1/besselI(kappa, 0))^10 * exp(kappa * sum(cos(radians-2.39))-kappa)
kappaFrame <- data.frame(x=kappa, y=distKappa)
plot(y=distKappa,x=kappa)
mean(distKappa)
which.max(distKappa)
kappa[2126]

test <- sort(table(distKappa))
head(test)

plot(density(distKappa))

kappaVals <- data.frame(x=distKappa)
ggplot(kappaVals, aes(x)) + geom_density()

# 4
hist(rgamma(10000, shape = 4, rate = 2))
mean(rgamma(10000, shape = 4, rate = 2))

bseq <- seq(1,3, 0.01)
prob <- 0
j <- 0
for(i in bseq){
  j <- j+1
  prob[j] <-  pgamma(q = 5, 4*i, i) - pgamma(q = 3, 4*i, i) 
}
prob <- data.frame(x=bseq, y=prob)
ggplot(prob, aes(x=x, y=y)) + geom_line() + ggtitle("Percentage of values in interval 3-5\n for different values of beta")

hist(rgamma(10000, shape = 4*1.8, rate = 1.8))
pgamma(q = 5, 4*1.8, 1.8) - pgamma(q = 3, 4*1.8, 1.8)

# c
ysum <- sum(120342 + 235967 + 243745 + 197452 + 276935 + 157222) / 100000
 # prior
prior <- data.frame(y=rgamma(10000, shape = 4*1.81, rate = 1.81))
plot(density(rgamma(10000, shape = 4*1.81, rate = 1.81)), ylim=c(0, 1.15))
plot(x=seq(0,10, 0.01), y=dgamma(seq(0, 10, by=0.01), 4*1.81, 1.81), type="l",
     ylim=c(0, 1.15), col="royalblue")
 # posterior
hist(rgamma(10000, shape = 8+19, rate = 2+12.31663))
lines(x=seq(0,10, 0.01),y=dgamma(seq(0, 10, by=0.01), 4*1.81 + 19, 1.81 + 12.31663), 
      type="l", col="darkorange")
lines(density(rgamma(10000, shape = 4*1.81+19, rate = 4*1.81+12.31663)), ylim=c(0, 1.15))


posterior <- data.frame(y=rgamma(10000, shape = 8+19, rate = 2+12.31663))
library(ggplot2)
library(scales)
ggplot(posterior, aes(y)) + scale_y_continuous(labels = percent_format())+
  geom_freqpoly(binwidth=0.05)  + geom_freqpoly(data=prior, aes(y), binwidth=0.05) 

lines(density(rgamma(10000, shape = 8+19, rate = 2+12.31663)))
pgamma(q = 5, 4*1.8+19, 1.8+12.31663) - pgamma(q = 3, 4*1.8+19, 1.8+12.31663)





