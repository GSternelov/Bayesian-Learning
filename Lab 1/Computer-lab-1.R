library(ggplot2)
# 1.a)
set.seed(12345)
beta10 <- rbeta(10, 16, 8) 
beta100 <- rbeta(100, 16, 8)
beta1000 <- rbeta(1000, 16, 8)
#sd(beta10); sd(beta100) ; sd(beta1000)
true_sd <- sqrt((16*8)/((16+8)^2*(16+8+1)))
#mean(beta10); mean(beta100) ; mean(beta1000)
true_mean <- 16/(16+8) # True value, (alpha + s) / ((alpha + s) + (beta + f))
# 1.b)
beta_10_val <-length(subset(beta10, beta10 < 0.4)) / 10
beta_100_val <- length(subset(beta100, beta100 < 0.4)) / 100
beta_1000_val <- length(subset(beta1000, beta1000 < 0.4)) / 1000
true_prob <- pbeta(0.4, 16, 8)
log_odds <- log(beta1000 / (1-beta1000))
log_odds_frame <- data.frame(log_odds=log_odds)
ggplot(log_odds_frame, aes(log_odds)) + geom_histogram(binwidth=0.1, fill="royalblue") + theme_bw() + ggtitle("Histogram for 1000 values from \n the posterior distribution of the log-odds")
density(log_odds)
y <- c(14,25,45,25,30,33,19,50,34,67)
sigma2 <- 0
set.seed(12345)
for(i in 1:10000){
  chisq <- rchisq(1, 10)
  sigma2[i] <- (sum((log(y)-3.5)^2)) / chisq
}
sigma2_frame <- data.frame(sigma2=sigma2)
ggplot(sigma2_frame, aes(sigma2)) + geom_histogram(binwidth=0.1, fill="royalblue") + theme_bw() + ggtitle("Histogram for 10 000 values from \n the posterior distribution of sigma^2")
max(sigma2)
#mean(sigma2)
#var(sigma2)
# theoretical mean
#c(theoretical_mean=(sum((log(y)-3.5)^2))/(10-2))
# theoretical variance
#c(theoretical_variance=(2*10^2 * ((sum((log(y)-3.5)^2))/10)^2 ) / (8^2 * 6))
c(Mean=mean(sigma2))
c(Variance=var(sigma2))
c(theoretical_mean=(sum((log(y)-3.5)^2))/(10-2))
c(theoretical_variance=(2*10^2 * ((sum((log(y)-3.5)^2))/10)^2 ) / (8^2 * 6))
G <- 2*pnorm((sqrt(sigma2)/sqrt(2)))-1
G_frame <- data.frame(G=G)
ggplot(G_frame, aes(G)) + geom_histogram(binwidth=0.008, fill="royalblue") + theme_bw() + ggtitle("Histogram for the posterior distribution of the Gini coefficient G")
# a)
radians <- c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)
kappa <- seq(from=0, to=5, by=0.001)
distKappa <- (1/besselI(kappa, 0))^10 * exp(kappa * sum(cos(radians-2.39))-kappa)
kappaFrame <- data.frame(x=kappa, y=distKappa)
maxkappa <- which.max(distKappa)
modeKappa <- kappa[2126]
ggplot(kappaFrame, aes(x=x,y=y)) + geom_line(col="royalblue", size=1.5) + geom_vline(xintercept=modeKappa,col="royalblue",size=1.25) +
  ggtitle("Posterior distribution of k \n mode = 2.125") + theme_bw() +ylab("frequency") + xlab("kappa")
bseq <- seq(1,3, 0.01)
prob <- 0
j <- 0
for(i in bseq){
  j <- j+1
  prob[j] <-  (pgamma(q = 5, 4*i, i) - pgamma(q = 3, 4*i, i) ) *100
}
prob <- data.frame(x=bseq, y=prob)
ggplot(prob, aes(x=x, y=y)) + geom_line(col="royalblue", size=1.35) + ggtitle("Percentage of values in interval 3-5\n for different values of beta") + ylab("Percentage between 3-5") + xlab("Beta value") + geom_point(data=data.frame(y=0.5004194*100,x=1.8), size=4, col="darkorange") + theme_bw()
# c
ysum <- sum(120342 + 235967 + 243745 + 197452 + 276935 + 157222) / 100000
 # prior
prior <- data.frame(x=seq(0,10, 0.01), y=dgamma(seq(0, 10, by=0.01), 4*1.81, 1.81))
posterior <- data.frame(x=seq(0,10, 0.01),y=dgamma(seq(0, 10, by=0.01), 4*1.81 + 19, 1.81 + 12.31663))
ggplot(prior, aes(x=x, y=y)) + geom_line(col="royalblue", size=1.3) + ylim(0,1.2) + geom_line(data=posterior, aes(x=x, y=y), col="darkorange", size=1.3) + theme_bw() + ggtitle("Probability density functions \n Prior(blue) compared to Posterior(orange)") + xlab("lambda")
## NA
