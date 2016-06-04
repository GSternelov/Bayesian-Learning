## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
# 1. Read in libraries
library(plyr); library(dplyr)
library(rjags)
library(coda)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
# 2. Read in data
res2015 <- read.csv("C:\\Users\\Gustav\\Documents\\Bayesian-Learning\\Project\\Allsvenskan2015.csv", sep=",", header = TRUE)
res2015 <- res2015[,c(21,5,14,23)]
names(res2015) <- c("Home.Team",  "Away.Team", "Home.Goals", "Away.Goals")
res2016 <- read.csv("C:\\Users\\Gustav\\Documents\\Bayesian-Learning\\Project\\Allsv2016.csv", sep=";", header = TRUE)
res1516 <- data.frame(bind_rows(res2015, res2016[1:96,]))
head(res1516)
teams <- unique(c(res1516$Home.Team, res1516$Away.Team))

## ---- echo=FALSE,fig.width=4, fig.height=4, fig.align='center', message=FALSE, warning=FALSE----
library(png)
img <- readPNG("C:\\Users\\Gustav\\Documents\\Bayesian-Learning\\Project\\DAG_final.PNG")
grid.raster(img)

## ---- echo=FALSE, warning=FALSE, message=FALSE, fig.height=2.7, fig.width=6----
load("C:/Users/Gustav/Documents/Bayesian-Learning/samples.Rda")
load("C:/Users/Gustav/Documents/Bayesian-Learning/samplesMatrix.Rda")
conv_plots <- function(indexParam, nameParam){
  trace <- ggplot(ms1df, aes(x=1:30000, y=ms1df[,indexParam])) + geom_line(col="royalblue") + theme_bw() + ylab(nameParam)
  HPD_lower <- HPDinterval(as.mcmc(ms1[,indexParam]), prob=0.95)[1]
  HPD_upper <- HPDinterval(as.mcmc(ms1[,indexParam]), prob=0.95)[2]
  density <- ggplot(ms1df, aes(x=ms1df[,indexParam])) + geom_density() + theme_bw() + geom_vline(xintercept=c(HPD_lower, HPD_upper), linetype=2)
  hist <- ggplot(ms1df, aes(x=ms1df[,indexParam])) + geom_histogram(aes(y=..density..),fill="royalblue", binwidth=0.01) + theme_bw() + 
    geom_vline(xintercept=c(HPD_lower, HPD_upper), linetype=2) + xlab(nameParam)
  ACF <- data.frame(acf=as.numeric(acf(ms1df[,indexParam], plot = FALSE, lag.max = 40)$acf), lag=0:40)
  Efficiency <- ggplot(ACF, aes(x=lag, y=acf)) + geom_bar(fill="royalblue", stat="identity") + theme_bw()
  text <- nameParam
  plots <- grid.arrange(trace, hist, Efficiency, ncol=2, top=textGrob(text,gp=gpar(fontsize=16,font=3)))
  return(plots)
}
hmInt <- conv_plots(4, "Home intercept")

## ---- echo=FALSE, fig.height=2.7, fig.width=6, warning=FALSE, message=FALSE----
awInt <- conv_plots(1, "Away intercept") 

## ---- echo=FALSE, fig.height=2.7, fig.width=6, warning=FALSE, message=FALSE----
skill <- conv_plots(3, "Mu - Teams")

## ---- echo=FALSE, fig.height=2.7, fig.width=6, warning=FALSE, message=FALSE----
sigma <- conv_plots(2, "Sigma - Teams")

## ---- echo=FALSE, fig.height=2.7, fig.width=6, warning=FALSE, message=FALSE----
PekingSkill <- conv_plots(10, "Skill - IFK Norrköping")

## ---- echo=FALSE, fig.height=2.7, fig.width=6, warning=FALSE, message=FALSE----
DIFSkill <- conv_plots(9, "Skill - Djurgårdens IF")

## ---- echo=FALSE, fig.height=5.5, fig.width=6----------------------------
teamSkill <- ms1[,5:22]
teamSkillHome <- (teamSkill - rowMeans(teamSkill)) + ms1[, "homeInt"]
teamSkillHome <- exp(teamSkillHome)
colnames(teamSkillHome) <- teams
teamSkillHome <- teamSkillHome[,order(colMeans(teamSkillHome), decreasing=T)]

teamSkillHome_frame <- data.frame(HPDinterval(as.mcmc(teamSkillHome), prob=0.95), median = as.numeric(apply(teamSkillHome, 2, FUN = median)) )
teamSkillHome_frame$team <- row.names(teamSkillHome_frame)
teamSkillHome_frame[,5:6] <-  HPDinterval(as.mcmc(teamSkillHome), prob=0.65)
ggplot(teamSkillHome_frame, aes(y=median, x=1:18)) + theme_bw() + geom_pointrange(aes(ymin = lower, ymax = upper),col="black", size=0.5) + coord_flip() + geom_text(aes(label=team),nudge_x = 0.5) + geom_pointrange(aes(ymin = V5, ymax =  V6), size=0.95,col="royalblue") +
  scale_x_reverse(breaks=c(1,4,8,12,16)) + ylab("Expected number of goals") + xlab("") + ggtitle("Expected number of goals per home game \nBlack line = 95% HPD, Blue line = 65% HPD")

## ---- echo=FALSE, fig.height=5.5, fig.width=6----------------------------
teamSkillAway <- (teamSkill - rowMeans(teamSkill)) + ms1[, "awayInt"]
teamSkillAway <- exp(teamSkillAway)
colnames(teamSkillAway) <- teams
teamSkillAway <- teamSkillAway[,order(colMeans(teamSkillAway), decreasing=T)]

teamSkillAway_frame <- data.frame(HPDinterval(as.mcmc(teamSkillAway), prob=0.95), median = as.numeric(apply(teamSkillAway, 2, FUN = median)) )
teamSkillAway_frame$team <- row.names(teamSkillAway_frame)
teamSkillAway_frame[,5:6] <-  HPDinterval(as.mcmc(teamSkillAway), prob=0.65)
ggplot(teamSkillAway_frame, aes(y=median, x=1:18)) + theme_bw() + geom_pointrange(aes(ymin = lower, ymax = upper), col="black", size=0.5) +  coord_flip() + geom_text(aes(label=team),nudge_x = 0.5) + geom_pointrange(aes(ymin = V5, ymax =  V6), size=0.95,col="royalblue") +
  scale_x_reverse(breaks=c(1,4,8,12,16)) + ylab("Expected number of goals") + xlab("")+ ggtitle("Expected number of goals per away game \nBlack line = 95% HPD, Blue line = 65% HPD")

## ----  echo=FALSE, fig.height=5.5, fig.width=6, warning=FALSE------------
load("C:/Users/Gustav/Documents/Bayesian-Learning/SeasonGoals.Rda")
load("C:/Users/Gustav/Documents/Bayesian-Learning/SeasonPred.Rda")
load("C:/Users/Gustav/Documents/Bayesian-Learning/SimGoalsHome.Rda")
load("C:/Users/Gustav/Documents/Bayesian-Learning/SimGoalsAway.Rda")
res1516$MatchResult <- sign(res1516$Home.Goals - res1516$Away.Goals) 

MatchResults <- data.frame(matrix(vector(), 30000, 480))
for(i in 1:480){
  MatchResults[,i] <- sign(PredGoalsHome[,i] - PredGoalsAway[,i]) 
}
MatchResults <- data.frame(t(MatchResults))

Allres1516 <- (bind_rows(res2015, res2016))
AllRes <- cbind(Allres1516[1:336,1:2],MatchResults[1:336,])

seasonInt <- data.frame(matrix(vector(), 18, 6))
for(i in 1:18){
  seasonInt[i,2:3] <- data.frame(HPDinterval(as.mcmc(as.numeric(seasonP[,i])), prob=0.95) )
  seasonInt[i,4:6] <- data.frame(HPDinterval(as.mcmc(as.numeric(seasonP[,i])), prob=0.65), 
                                 median=as.numeric(apply(data.frame(as.numeric(seasonP[,i])), 2, FUN = median)) )
}
seasonInt[,1] <- data.frame(table(AllRes$Home.Team))[,1]
names(seasonInt) <- c("Team", "Lower95", "Upper95", "Lower65","Upper65", "Median")

seasonInt$Actual <- as.numeric(table(res1516$Home.Team,res1516$MatchResult)[,3] *3) + as.numeric(table(res1516$Home.Team, res1516$MatchResult)[,2]) +
  as.numeric(table(res1516$Away.Team, res1516$MatchResult)[,1] *3) + as.numeric(table(res1516$Away.Team, res1516$MatchResult)[,2]) 

seasonInt <-seasonInt[order(-seasonInt$Median),]
cols <- c("Median"="royalblue","Actual"="darkorange")
ggplot(seasonInt, aes(y=Median,x=1:18))+scale_x_reverse(breaks=c(1,4,8,12,16))+theme_bw()+coord_flip()+geom_text(aes(label=Team),nudge_x=0.5)+ geom_pointrange(aes(ymin=Lower95,ymax=Upper95),col="black")+ geom_pointrange(aes(ymin=Lower65,ymax=Upper65,col="Median"), size=0.8) + geom_point(aes(y=Actual,col="Actual"), size=3.5) + scale_colour_manual(name="",values=cols) + theme(legend.position=c(.9, .1),legend.background = element_rect(color = "black",   fill = "white", size = 1, linetype = "solid"), legend.text = element_text(size = 14)) + ggtitle("Simulated number of points vs actual \nBlack line = 95% HPD, Blue line = 65% HPD") + ylab("Simulated number of points") + xlab("")

## ----  echo=FALSE, fig.height=5.5, fig.width=6, warning=FALSE------------
seasonGoalInt <- data.frame(matrix(vector(), 18, 6))
for(i in 1:18){
  seasonGoalInt[i,2:3] <- data.frame(HPDinterval(as.mcmc(as.numeric(SeasonGoals[,i])), prob=0.95) )
  seasonGoalInt[i,4:6] <- data.frame(HPDinterval(as.mcmc(as.numeric(SeasonGoals[,i])), prob=0.65), 
                                 median=as.numeric(apply(data.frame(as.numeric(SeasonGoals[,i])), 2, FUN = mean)) )
}
seasonGoalInt[,1] <- data.frame(table(AllRes$Home.Team))[,1]
names(seasonGoalInt) <- c("Team", "Lower95", "Upper95", "Lower65","Upper65", "Median")

seasonGoalInt$Actual <- aggregate(res1516$Home.Goals, by=list(Team=res1516$Home.Team), FUN=sum)[,2] + 
  aggregate(res1516$Away.Goals, by=list(Team=res1516$Away.Team), FUN=sum)[,2]
  
seasonGoalInt <-seasonGoalInt[order(-seasonGoalInt$Median),]
cols <- c("Median"="royalblue","Actual"="darkorange")
ggplot(seasonGoalInt, aes(y=Median,x=1:18))+scale_x_reverse(breaks=c(1,4,8,12,16))+theme_bw()+coord_flip()+geom_text(aes(label=Team),nudge_x=0.5)+ geom_pointrange(aes(ymin=Lower95,ymax=Upper95),col="black")+ geom_pointrange(aes(ymin=Lower65,ymax=Upper65,col="Median"), size=0.8) + geom_point(aes(y=Actual,col="Actual"), size=3.5) + scale_colour_manual(name="",values=cols) + theme(legend.position=c(.9, .1),legend.background = element_rect(color = "black",   fill = "white", size = 1, linetype = "solid"),legend.text = element_text(size = 14)) +  ggtitle("Simulated number of goals vs actual \nBlack line = 95% HPD, Blue line = 65% HPD")+ ylab("Simulated number of goals") + xlab("")

## ---- echo=FALSE---------------------------------------------------------
res2016$MatchResult <- sign(res2016$Home.Goals - res2016$Away.Goals) 
AllFut <- cbind(Allres1516[337:480,1:2],MatchResults[337:480,])
load("C:/Users/Gustav/Documents/Bayesian-Learning/futurePred.Rda")
seasonPred <- data.frame(matrix(vector(), 16, 6))
for(i in 1:16){
  seasonPred[i,2:3] <- data.frame(HPDinterval(as.mcmc(as.numeric(seasonFut[,i])), prob=0.95) )
  seasonPred[i,4:6] <- data.frame(HPDinterval(as.mcmc(as.numeric(seasonFut[,i])), prob=0.65), 
                                 median=as.numeric(apply(data.frame(as.numeric(seasonFut[,i])), 2, FUN = median)) )
}
seasonPred[,1] <- data.frame(table(AllFut$Home.Team))[,1]
names(seasonPred) <- c("Team", "Lower95", "Upper95", "Lower65","Upper65", "Median")
seasonPred$Current <- as.numeric(table(res2016$Home.Team,res2016$MatchResult)[,3] *3) + as.numeric(table(res2016$Home.Team, res2016$MatchResult)[,2]) +
  as.numeric(table(res2016$Away.Team, res2016$MatchResult)[,1] *3) + as.numeric(table(res2016$Away.Team, res2016$MatchResult)[,2]) 
seasonPred$Total <- seasonPred$Median + seasonPred$Current
Current = seasonPred$Current
seasonPred <-data.frame(Team=data.frame(table(AllFut$Home.Team))[,1],Current=seasonPred$Current, Median=seasonPred$Median, Total=seasonPred$Total, Lower95=seasonPred$Lower95 + Current,  Upper95=seasonPred$Upper95 + Current,Lower65=seasonPred$Lower65 + Current,Upper65=seasonPred$Upper65 + Current)
seasonPred <-seasonPred[order(-seasonPred$Total),]
row.names(seasonPred) <- NULL
seasonPred

## ---- echo=FALSE---------------------------------------------------------
load("C:/Users/Gustav/Documents/Bayesian-Learning/SeasonPos.Rda")
prob <- data.frame(Team=data.frame(table(AllFut$Home.Team))[,1],Winner=0, Top4=0, Top8=0, Rel_Play_off=0, Relegation=0)
Win <- 1
Top4 <- seq(1,4,1)
Top8 <- seq(1,8,1)
Rel_Play_off <- 14
Relegation <- seq(15,16,1)
for(i in 1:16){
  prob[i,2] <- round(mean(as.numeric(SeasonPos[,i]) %in% Win) * 100,2)
  prob[i,3]  <-round(mean(as.numeric(SeasonPos[,i]) %in% Top4) * 100,2)
  prob[i,4] <- round(mean(as.numeric(SeasonPos[,i]) %in% Top8) * 100,2)
  prob[i,5]  <- round(mean(as.numeric(SeasonPos[,i]) %in% Rel_Play_off) * 100,2)
  prob[i,6]  <- round(mean(as.numeric(SeasonPos[,i]) %in% Relegation) * 100,2)
}
prob <-prob[order(-prob$Winner,-prob$Top4,-prob$Top8,-prob$Rel_Play_off,-prob$Relegation),]
row.names(prob) <- NULL
prob

## ---- echo=FALSE, warning=FALSE, message=FALSE, fig.height=2.7, fig.width=6----
skill <- conv_plots(6, "Skill - Kalmar FF")
skill <- conv_plots(7, "Skill - Falkenbergs FF")
skill <- conv_plots(8, "Skill - IFK Göteborg")
skill <- conv_plots(11, "Skill - AIK")
skill <- conv_plots(12, "Skill - GIF Sundsvall")
skill <- conv_plots(13, "Skill - Gefle IF")
skill <- conv_plots(14, "Skill - Helsingborgs IF")
skill <- conv_plots(15, "Skill - Örebro SK")
skill <- conv_plots(16, "Skill - IF Elfsborg")
skill <- conv_plots(17, "Skill - BK Häcken")
skill <- conv_plots(18, "Skill - Halmstads BK")
skill <- conv_plots(19, "Skill - Åtvidabergs FF")
skill <- conv_plots(20, "Skill - Malmö FF")
skill <- conv_plots(21, "Skill - Jönköpings Södra")
skill <- conv_plots(22, "Skill - Östersunds FK")

## ----code=readLines(knitr::purl("C:\\Users\\Gustav\\Documents\\Bayesian-Learning\\Project\\The report.Rmd",documentation = 1)), eval = FALSE----
## 
## 

