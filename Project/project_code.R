# 1. Read in libraries
library(plyr); library(dplyr)
library(tidyr)
library(rjags)
library(coda)
library(stringr)
library(ggplot2)
library(gridExtra)

# 2. Read in data
res2015 <- read.csv("C:\\Users\\Gustav\\Documents\\Bayesian-Learning\\Project\\Allsvenskan2015.csv", sep=",", header = TRUE)
res2015 <- res2015[,c(21,5,14,23)]
names(res2015) <- c("Home.Team",  "Away.Team", "Home.Goals", "Away.Goals")
res2016 <- read.csv("C:\\Users\\Gustav\\Documents\\Bayesian-Learning\\Project\\Allsv2016.csv", sep=";", header = TRUE)

res1516 <- (bind_rows(res2015, res2016))

# 3. Specify model
## First, a list with information from the data frame
teams <- unique(c(res1516$Home.Team, res1516$Away.Team))

Allsv_list <- list(HomeGoals=res1516$Home.Goals,AwayGoals=res1516$Away.Goals,
                   HomeTeam=as.numeric(factor(res1516$Home.Team, levels=teams)), AwayTeam=as.numeric(factor(res1516$Away.Team, levels=teams)),
                teams = length(teams), games = nrow(res1516))

## Second, write a string with the model for the JAGS function
MyModelString <- "model {
  for (i in 1:games){
    HomeGoals[i] ~ dpois(lambdaHome[HomeTeam[i],AwayTeam[i]])
    AwayGoals[i] ~ dpois(lambdaAway[HomeTeam[i],AwayTeam[i]])
  }
  for(home in 1:teams){
    for(away in 1:teams){
      lambdaHome[home, away] <- exp(homeInt + skill[home] - skill[away])
      lambdaAway[home, away] <- exp(awayInt + skill[away] - skill[home]) 
    }
  }
skill[1] <- 0
for(j in 2:teams) {
  skill[j] ~ dnorm(group_skill, group_precision)
}  
group_skill ~ dnorm(0, 0.0625)
group_precision <- 1 / pow(group_sigma, 2)
group_sigma ~ dunif(0, 3)
homeInt ~ dnorm(0, 0.0625)
awayInt ~ dnorm(0, 0.0625)
} "

# 4. Compile model and generate MCMC samples
## Compiling model
AllsvModel <- jags.model(textConnection(MyModelString), data=Allsv_list, n.chains=3, n.adapt=10000)
## Burning some samples
update(AllsvModel, 10000)
## Generate MCMC samples
s1 <- coda.samples(AllsvModel, variable.names=c("homeInt", "awayInt", "skill", "group_skill", "group_sigma"), n.iter=20000, thin=2)
## Merge the MCMC chains into one matrix
ms1 <- as.matrix(s1)
ms1df <- as.data.frame(ms1)

# 5. Analyze convergence of the parameters
### Think i have got three chains with 10000 iterations where half of each sample is classed as burn-in
## awayInt = 1, groupSigma=2, groupSkill=3, homeInt=4
## For the skill parameters:
##  [1] "Hammarby IF"      "Kalmar FF"        "Falkenbergs FF"   "IFK Göteborg"     "Djurgårdens IF"   "IFK Norrköping"   "AIK"             
#[8] "GIF Sundsvall"    "Gefle IF"         "Helsingborgs IF"  "Örebro SK"        "IF Elfsborg"      "BK Häcken"        "Halmstads BK"    
#[15] "Åtvidabergs FF"   "Malmö FF"         "Jönköpings Södra" "Östersunds FK"   

conv_plots <- function(indexParam){
  trace <- ggplot(ms1df, aes(x=1:30000, y=ms1df[,indexParam])) + geom_line(col="royalblue") + theme_bw()
  HPD_lower <- HPDinterval(as.mcmc(ms1[,indexParam]), prob=0.95)[1]
  HPD_upper <- HPDinterval(as.mcmc(ms1[,indexParam]), prob=0.95)[2]
  density <- ggplot(ms1df, aes(x=ms1df[,indexParam])) + geom_density() + theme_bw() + geom_vline(xintercept=c(HPD_lower, HPD_upper), linetype=2)
  hist <- ggplot(ms1df, aes(x=ms1df[,indexParam])) + geom_histogram(aes(y=..density..),fill="royalblue", binwidth=0.01) + theme_bw() + 
    geom_vline(xintercept=c(HPD_lower, HPD_upper), linetype=2)
  ACF <- data.frame(acf=as.numeric(acf(ms1df[,indexParam], plot = FALSE, lag.max = 40)$acf), lag=0:40)
  Efficiency <- ggplot(ACF, aes(x=lag, y=acf)) + geom_bar(fill="royalblue", stat="identity") + theme_bw()
  plots <- grid.arrange(trace, hist, Efficiency, ncol=2)
  return(plots)
}

awInt <- conv_plots(1)
sigma <- conv_plots(2)
skill <- conv_plots(3)
hmInt <- conv_plots(4)
Peking <- conv_plots(10)

# 6. Compare estimates of the skill parameters for the teams, home and away
## Home
teamSkill <- ms1[,5:22]
teamSkillHome <- (teamSkill - rowMeans(teamSkill)) + ms1[, "homeInt"]
teamSkillHome <- exp(teamSkillHome)
colnames(teamSkillHome) <- teams
teamSkillHome <- teamSkillHome[,order(colMeans(teamSkillHome), decreasing=T)]

teamSkillHome_frame <- data.frame(HPDinterval(as.mcmc(teamSkillHome), prob=0.95), median = as.numeric(apply(teamSkillHome, 2, FUN = median)) )
teamSkillHome_frame$team <- row.names(teamSkillHome_frame)
teamSkillHome_frame[,5:6] <-  HPDinterval(as.mcmc(teamSkillHome), prob=0.65)
ggplot(teamSkillHome_frame, aes(y=median, x=18:1)) + theme_bw() + geom_pointrange(aes(ymin = lower, ymax = upper),col="royalblue") +
  coord_flip() + geom_text(aes(label=team),nudge_x = 0.5) + geom_pointrange(aes(ymin = V5, ymax =  V6), size=1.05,col="royalblue")

## Away
teamSkillAway <- (teamSkill - rowMeans(teamSkill)) + ms1[, "awayInt"]
teamSkillAway <- exp(teamSkillAway)
colnames(teamSkillAway) <- teams
teamSkillAway <- teamSkillAway[,order(colMeans(teamSkillAway), decreasing=T)]

teamSkillAway_frame <- data.frame(HPDinterval(as.mcmc(teamSkillAway), prob=0.95), median = as.numeric(apply(teamSkillAway, 2, FUN = median)) )
teamSkillAway_frame$team <- row.names(teamSkillAway_frame)
teamSkillAway_frame[,5:6] <-  HPDinterval(as.mcmc(teamSkillAway), prob=0.65)
ggplot(teamSkillAway_frame, aes(y=median, x=18:1)) + theme_bw() + geom_pointrange(aes(ymin = lower, ymax = upper),col="royalblue") +
  coord_flip() + geom_text(aes(label=team),nudge_x = 0.5) + geom_pointrange(aes(ymin = V5, ymax =  V6), size=1.05,col="royalblue")


# 7. Predictions





