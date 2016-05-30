library(plyr); library(dplyr)
library(tidyr)

res2015 <- read.csv("C:\\Users\\Gustav\\Documents\\Bayesian-Learning\\Project\\Allsvenskan2015.csv", sep=",", header = TRUE)
res2015 <- res2015[,c(21,5,14,23)]
names(res2015) <- c("Home.Team",  "Away.Team", "Home.Goals", "Away.Goals")
res2016 <- read.csv("C:\\Users\\Gustav\\Documents\\Bayesian-Learning\\Project\\Allsv2016.csv", sep=";", header = TRUE)

res1516 <- (bind_rows(res2015, res2016))


# La Liga example
library(rjags)
library(coda)
library(mcmcplots)
library(stringr)


load("C:\\Users\\Gustav\\Documents\\Bayesian-Learning\\Project\\rasmus_baath_user_13_data_analysis_contest\\laliga.RData")

laliga$MatchResult <- sign(laliga$HomeGoals - laliga$AwayGoals) 

d <- na.omit(laliga)
teams <- unique(c(d$HomeTeam, d$AwayTeam))
seasons <- unique(d$Season)

data_list <- list(HomeGoals = d$HomeGoals, AwayGoals = d$AwayGoals, 
                  HomeTeam = as.numeric(factor(d$HomeTeam, levels=teams)),
                  AwayTeam = as.numeric(factor(d$AwayTeam, levels=teams)),
                  Season = as.numeric(factor(d$Season, levels=seasons)),
                  n_teams = length(teams), n_games = nrow(d), 
                  n_seasons = length(seasons))

old_par <- par(mfcol=c(2,1), mar=rep(2.2, 4))
hist(c(d$AwayGoals, d$HomeGoals), xlim=c(-0.5, 8), breaks = -1:9 + 0.5)
mean_goals <- mean(c(d$AwayGoals, d$HomeGoals))
hist(rpois(9999, mean_goals), xlim=c(-0.5, 8), breaks = -1:9 + 0.5)
# Allsvenskan data
hist(c(res1516$Away.Goals, res1516$Home.Goals), xlim=c(-0.5, 8), breaks = -1:9 + 0.5)
mean_goals <- mean(c(res1516$Away.Goals, res1516$Home.Goals))
hist(rpois(9999, mean_goals), xlim=c(-0.5, 8), breaks = -1:9 + 0.5)
par(old_par)

# Convenience function to generate the type of column names Jags outputs.
col_name <- function(name, ...) {
  paste0(name, "[", paste(..., sep=",") , "]")
}

m1_string <- "model {
for(i in 1:n_games) {
HomeGoals[i] ~ dpois(lambda_home[HomeTeam[i],AwayTeam[i]])
AwayGoals[i] ~ dpois(lambda_away[HomeTeam[i],AwayTeam[i]])
}

for(home_i in 1:n_teams) {
for(away_i in 1:n_teams) {
lambda_home[home_i, away_i] <- exp(baseline + skill[home_i] - skill[away_i])
lambda_away[home_i, away_i] <- exp(baseline + skill[away_i] - skill[home_i])
}
}

skill[1] <- 0
for(j in 2:n_teams) {
skill[j] ~ dnorm(group_skill, group_tau)
}  

group_skill ~ dnorm(0, 0.0625)
group_tau <- 1 / pow(group_sigma, 2)
group_sigma ~ dunif(0, 3)
baseline ~ dnorm(0, 0.0625)
}"

# Compiling model 1
m1 <- jags.model(textConnection(m1_string), data=data_list, n.chains=3, n.adapt=5000)
# Burning some samples on the altar of the MCMC god
update(m1, 5000)
# Generating MCMC samples
s1 <- coda.samples(m1, variable.names=c("baseline", "skill", "group_skill", "group_sigma"), n.iter=10000, thin=2)
# Merging the three MCMC chains into one matrix
ms1 <- as.matrix(s1)

plot(s1[,col_name("skill", which(teams == "FC Sevilla"))])
plot(s1[,col_name("skill", which(teams == "FC Valencia"))])

## @knitr fig.height=5
# Plots histograms over home_goals, away_goals, the difference in goals and a barplot over the credible match results.
plot_goals <- function(home_goals, away_goals) {
  n_matches <- length(home_goals)
  goal_diff <- home_goals - away_goals
  match_result <- ifelse(goal_diff < 0, "away_win", ifelse(goal_diff > 0, "home_win", "equal"))
  hist(home_goals, xlim=c(-0.5, 10), breaks=(0:100) - 0.5)
  hist(away_goals, xlim=c(-0.5, 10), breaks=(0:100) - 0.5)
  hist(goal_diff, xlim=c(-6, 6), breaks=(-100:100) - 0.5 )
  barplot(table(match_result) / n_matches , ylim=c(0, 1))
}

# Simulates game goals scores using the MCMC samples from the m1 model.
plot_pred_comp1 <- function(home_team, away_team, ms) {
  old_par <- par(mfrow = c(2, 4))
  baseline <- ms[, "baseline"]
  home_skill <- ms[, col_name("skill", which(teams == home_team))]
  away_skill <- ms[, col_name("skill", which(teams == away_team))]
  home_goals <- rpois(nrow(ms),  exp(baseline +  home_skill - away_skill))
  away_goals <- rpois(nrow(ms),  exp(baseline +  away_skill - home_skill))
  plot_goals(home_goals, away_goals)
  home_goals <- d$HomeGoals[ d$HomeTeam == home_team & d$AwayTeam == away_team]
  away_goals <- d$AwayGoals[ d$HomeTeam == home_team & d$AwayTeam == away_team]
  plot_goals(home_goals, away_goals)
  par(old_par)
}

plot_pred_comp1("FC Valencia", "FC Sevilla", ms1)

plot_pred_comp1("FC Sevilla", "FC Valencia",ms1)

# model 2
m2_string <- "model {
for(i in 1:n_games) {
HomeGoals[i] ~ dpois(lambda_home[HomeTeam[i],AwayTeam[i]])
AwayGoals[i] ~ dpois(lambda_away[HomeTeam[i],AwayTeam[i]])
}

for(home_i in 1:n_teams) {
for(away_i in 1:n_teams) {
lambda_home[home_i, away_i] <- exp( home_baseline + skill[home_i] - skill[away_i])
lambda_away[home_i, away_i] <- exp( away_baseline + skill[away_i] - skill[home_i])
}
}

skill[1] <- 0 
for(j in 2:n_teams) {
skill[j] ~ dnorm(group_skill, group_tau)
}

group_skill ~ dnorm(0, 0.0625)
group_tau <- 1/pow(group_sigma, 2)
group_sigma ~ dgamma(0.1,0.1) 

home_baseline ~ dnorm(0, 0.0625)
away_baseline ~ dnorm(0, 0.0625)
}"


## @knitr results='hide'
m2 <- jags.model(textConnection(m2_string), data=data_list, n.chains=3, n.adapt=5000)
update(m2, 5000)
s2 <- coda.samples(m2, variable.names=c("home_baseline", "away_baseline","skill", "group_sigma", "group_skill"), n.iter=10000, thin=2)
ms2 <- as.matrix(s2)

plot(s2[,"home_baseline"])
plot(s2[,"away_baseline"])

plotPost(exp(ms2[,"home_baseline"]) - exp(ms2[,"away_baseline"]), compVal=0, 
         xlab="Home advantage in number of goals")

Gm2 <- jags.model(textConnection(m2_string), data=data_list, n.chains=3, n.adapt=5000)
update(Gm2, 5000)
Gs2 <- coda.samples(Gm2, variable.names=c("home_baseline", "away_baseline","skill", "group_sigma", "group_skill"), n.iter=10000, thin=2)
Gms2 <- as.matrix(Gs2)


