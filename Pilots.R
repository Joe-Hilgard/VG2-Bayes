# Load packages
library(BayesFactor)
library(magrittr)

# Function for pilot test w/ null-equivalent region.
equivTest = function(n1, n2, t, nullInterval=NULL, mu=0, rscale=.5) {
  x = rnorm(n1); x = x - mean(x); x = x / sd(x) # center & scale @ 0
  y = rnorm(n2); y = y - mean(y); y = y / sd(y) # center and scale at 0
  # Assume pooled variance of 1 (as we have scaled it as such)
  # Thus, SE depends only on sample size
  SE = sqrt(1/n1 + 1/n2)
  y = y + t*SE # center at t*SE
  return(ttestBF(x, y, nullInterval=nullInterval, paired=F, mu=mu, rscale=rscale))
}
equivTestPaired = function(N, t, nullInterval=NULL, rscale=.5) {
  x = rnorm(N)
  diff = rnorm(N, sd = sqrt(N))
  diff = diff / sd(diff) * sqrt(N) # rescale 
  diff = diff - mean(diff) + t # recenter
  y = x + diff
  # check that the t-value reproduces appropriately
  stopifnot(round(t.test(y, x, paired=T)$statistic, 2) == round(t, 2)) 
  return(ttestBF(x, y, nullInterval=nullInterval, paired=T, rscale=rscale))
}
#welch's t
welch.t = function(m1, m2, sd1, sd2, n1, n2) {
  sp = sqrt( ((n1-1)*sd1^2 + (n2-1) * sd2^2) / (n1 + n2 -2) )
  se_diff = sp * sqrt(1/n1 + 1/n2)
  mean_diff = m1-m2
  return(mean_diff/se_diff)
}

###
# Arriaga et al., 2008 and the insufficient pilot test
###
equivTestPaired(20, 2.63, NULL) #Difficulty
equivTestPaired(20, 2.27, NULL) #Competence
1/equivTestPaired(20, 1.67, NULL) #Discomfort
1/equivTestPaired(20, 1.56, NULL) #Realism
1/equivTestPaired(20, 1.32, NULL) #Frustration
1/equivTestPaired(20, 1.29, NULL) #Pleasure
1/equivTestPaired(20, 1.24, NULL) #Action
1/equivTestPaired(20, 1.14, NULL) #Disorientation
1/equivTestPaired(20, .89, NULL) #Excitement
1/equivTestPaired(20, .86, NULL) #Identification
1/equivTestPaired(20, .83, NULL) #Satisfaction
1/equivTestPaired(20, .79, NULL) #Boredom
1/equivTestPaired(20, .53, NULL) #Presence
1/equivTestPaired(20, .48, NULL) #Involvement

#####
# Valadez & Ferguson and the insufficient pilot test
#####
# t-tests for each pairwise comparison
# All values as reported in Table 1 (means, SDs, Ns)

# Cell sizes: RDR-v, 15; RDR-nv, 10; FIFA, 15.
# RDR "violent" vs RDR "nonviolent"
  # Difficulty
(t = welch.t(11.5, 9.6, 3.31, 1.89, 15, 10))
equivTest(15, 10, t, NULL)
  # Pace
(t = welch.t(12.0, 10.4, 3.38, 2.71, 15, 10))
1/equivTest(15, 10, t, NULL) 
  # Competitveness
(t = welch.t(8.53, 5.60, 3.50, 1.17, 15, 10))
equivTest(15, 10, t, NULL) 
# RDR "violent" vs FIFA
  # Difficulty
(t = welch.t(11.5, 13.2, 3.31, 3.01, 15, 15))
1/equivTest(15, 10, t, NULL) 
  # Pace
(t = welch.t(12.0, 14.3, 3.38, 2.89, 15, 15))
equivTest(15, 10, t, NULL) 
  # Competitveness
(t = welch.t(8.53, 8.47, 3.50, 3.42, 15, 15))
1/equivTest(15, 10, t, NULL) 
# RDR "nonviolent" vs FIFA 
  # Difficulty
(t = welch.t(9.6, 13.2, 1.89, 3.01, 10, 15))
equivTest(15, 10, t, NULL) 
  # Pace
(t = welch.t(10.4, 14.3, 2.71, 2.89, 10, 15))
equivTest(15, 10, t, NULL) 
  # Competitveness
(t = welch.t(5.60, 8.47, 1.17, 3.42, 10, 15))
equivTest(15, 10, t, NULL) 

###
# Adachi & Willoughby 2010 and the insufficient pilot test
###
# pilot 1, 14 subjects (paired t-test)
# F-tests provided in Table 1. Square-root taken to convert to t statistic
equivTestPaired(14, sqrt(19.31), NULL)  # Violence
1/equivTestPaired(14, -sqrt(0.21), NULL)  # Competition
1/equivTestPaired(14, sqrt(0.35), NULL)   # Difficulty
1/equivTestPaired(14, sqrt(0.64), NULL)   # Pace

# experiment 1, 21 subs in 2 cells each
equivTest(21, 21, sqrt(61.75), NULL)  # Violence
1/equivTest(21, 21, -sqrt(0.15), NULL)  # Competition
1/equivTest(21, 21, -sqrt(2.54), NULL)  # Difficulty
1/equivTest(21, 21, -sqrt(2.56), NULL)  # Pace 
# Note that using means&SDs gets you a different answer 
  # presumably b/c of the Gender variable.

###
# Anderson et al., 2004 and the insufficient pilot test
###

# In which it is argued that Glider Pro and Marathon 2 are equivalent
# 120 subjects across 10 games = 12 participants per game
# Violence difference between Marathon 2 and Glider Pro is r = .842 [.392, .854]
# They report means and mean squared errors, e.g. the average variance around each mean
# We'll have to assume that MSE as the variance within each cell.
# All values pulled from Table 1

# Generate t-values, then plug those ts into equivTest:
# Difficulty
t.dif = welch.t(4.25, 3.62, sqrt(2.37), sqrt(2.37), 12, 12)
# Enjoyment
t.enj = welch.t(3.69, 3.94, sqrt(2.35), sqrt(2.35), 12, 12)
# Action
t.act = welch.t(3.67, 2.31, sqrt(2.01), sqrt(2.01), 12, 12)
# Frustration
t.fru = welch.t(4.25, 4.75, sqrt(2.38), sqrt(2.38), 12, 12)
# Violence
t.vio = welch.t(4.86, 1.41, sqrt(2.38), sqrt(2.38), 12, 12)
# Bayes factors: 
equivTest(12, 12, t.act, nullInterval=NULL) # BF here favors the alternative 2.41-to-1
1/equivTest(12, 12, t.dif, nullInterval=NULL)
1/equivTest(12, 12, t.fru, nullInterval=NULL)
1/equivTest(12, 12, t.enj, nullInterval=NULL)
equivTest(12, 12, t.vio, nullInterval=NULL)

# Glider Pro & Marathon 2 might be equivalent on difficulty, enjoyment, frustration, but might differ in pace of action
# Experiments 2 and 3 of this manuscript have larger sample sizes and use similar manipulations.
#   Experiment 2 does not report any check of these matching variables.
#   Experiment 3 reports an absence of statistical significance, but the manipulation is a little different,
#     and the necessary summary statistics are not reported.

### CUT MATERIALS ###
# # I have decided not to attempt evaluating Adachi's second pilot study
# # It is too difficult to conduct the full 2x2 at this time.
# # pilot study 2, N = 19 within 
# dat = data.frame("Condition" = rep(c("C_NV", "C_V", "LC_NV", "LC_V"), each=4),
#                  "Violence" = rep(c(0, 1, 0, 1), each=4),
#                  "Compete" = rep(c(1, 1, 0, 0), each=4),
#                  "Outcome" = rep(c("Violence", "Competitiveness", "Difficulty", "Pace"),4),
#                  "Mean" = c(1.52, 5.86, 3.63, 4.74,
#                             5.37, 6.32, 4.47, 5.42,
#                             1, 1.36, 4.68, 4.57,
#                             6.42, 3.18, 4.42, 5.11),
#                  "SD" = c(.7, .76, 1.38, 1.05,
#                           .9, .58, 1.58, .90,
#                           0, .44, 1.16, 1.26,
#                           .84, 1.67, 1.74, 1.29)
# )
# # generate data
# # Rhis is challenging because it's repeated-measures.
# N = 19
# genData = c()
# for (i in dat$Condition) {
#   for (j in dat$Outcome) {
#     activeRow = dat[dat$Condition == i & dat$Outcome == j,]
#     cellMean = activeRow$Mean
#     cellSD = activeRow$SD
#     cellDat = rnorm(N, cellMean, cellSD)
#     cellDat = cellDat - mean(cellDat) + cellMean # center precisely
#     if (cellSD > 0) cellDat = cellDat / sd(cellDat) # scale precisely
#     newData = data.frame("Condition" = i, "Outcome" = j, "Value" = cellDat)
#     rbind(genData, newData)
#   }
# }
# 
# # experiment 2 it would be necessary to have the raw data to perform anovaBF()
# # No, I could simulate this still... I need to pull from a normal,
#   # Then standardize & center to appropriate means, that's all.
# # Doing the full model-comparison BF could get messy, but would be very interesting.
# dat = data.frame("x" = rep(c("A", "B", "C"), each=15))
# dat$y = rnorm(45)
# dat$y = dat$y + .3*(dat$x=="A") - .5*(dat$x == "C")
# anovaBF(y ~ x, data=dat)
