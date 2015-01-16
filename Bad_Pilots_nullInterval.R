# Function for pilot test w/ null-equivalent region.
equivTestPaired = function(N, t, nullInterval=NULL, rscale=.5) {
  x = rnorm(N)
  diff = rnorm(N, sd = sqrt(N))
  diff = diff / sd(diff) * sqrt(N) # rescale 
  diff = diff - mean(diff) + t # recenter
  y = x + diff
  return(ttestBF(x, y, nullInterval=nullInterval, paired=T, rscale=rscale))
}
equivTest = function(n1, n2, t, nullInterval=NULL, mu=0, rscale=.5) {
  x = rnorm(n1); x = x - mean(x); x = x / sd(x) # center & scale @ 0
  y = rnorm(n2)
  y = y/sd(y); y = y - mean(y) # center and scale at 0
  SE = sqrt(var(x)/n1 + var(y)/n2)
  y = y + t*SE # center at t*SE
  return(ttestBF(x, y, nullInterval=nullInterval, paired=F, mu=mu, rscale=rscale))
}
#welch's t
welch.t = function(m1, m2, sd1, sd2, n1, n2) {
  se_diff = sqrt(sd1^2/n1 + sd2^2/n2)
  mean_diff = m1-m2
  return(mean_diff/se_diff)
}
invertBF = function(model) {
  return(1/exp(model@bayesFactor[['bf']]))
}
# Load in BayesFactor package 
require(BayesFactor)

###
# Arriaga et al., 2008 and the insufficient pilot test
###
equivTestPaired(20, 2.63, c(-.1, .1)) #Difficulty
equivTestPaired(20, 2.27, c(-.1, .1)) #Competence
equivTestPaired(20, 1.67, c(-.1, .1)) #Discomfort
equivTestPaired(20, 1.56, c(-.1, .1)) #Realism
equivTestPaired(20, 1.32, c(-.1, .1)) #Frustration
equivTestPaired(20, 1.29, c(-.1, .1)) #Pleasure
equivTestPaired(20, 1.24, c(-.1, .1)) #Action
equivTestPaired(20, 1.14, c(-.1, .1)) #Disorientation
equivTestPaired(20, .89, c(-.1, .1)) #Excitement
equivTestPaired(20, .86, c(-.1, .1)) #Identification
equivTestPaired(20, .83, c(-.1, .1)) #Satisfaction
equivTestPaired(20, .79, c(-.1, .1)) #Boredom
equivTestPaired(20, .53, c(-.1, .1)) #Presence
equivTestPaired(20, .48, c(-.1, .1)) #Involvement

#####
# Valadez & Ferguson and the insufficient pilot test
#####
# t-tests for each pairwise comparison

# Cell sizes: RDR-v, 15; RDR-nv, 10; FIFA, 15.
# RDR "violent" vs RDR "nonviolent"
equivTest(15, 10, 1.82, c(-.1, .1)) # Difficulty
equivTest(15, 10, 1.31, c(-.1, .1)) # Pace
equivTest(15, 10, 3.00, c(-.1, .1)) # Competitveness
# RDR "violent" vs FIFA
equivTest(15, 10, 1.47, c(-.1, .1)) # Difficulty
equivTest(15, 10, 2.00, c(-.1, .1)) # Pace
equivTest(15, 10, .047, c(-.1, .1)) # Competitveness
# RDR "nonviolent" vs FIFA 
equivTest(15, 10, 3.45, c(-.1, .1)) # Difficulty
equivTest(15, 10, 3.43, c(-.1, .1)) # Pace
equivTest(15, 10, 3.00, c(-.1, .1)) # Competitveness

###
# Adachi & Willoughby 2010 and the insufficient pilot test
###
# pilot 1, 14 subjects (paired t-test)
equivTestPaired(14, 4.39, c(-.1, .1)) # Violence?
equivTestPaired(14, -.46, c(-.1, .1)) # Competition
equivTestPaired(14, .59, c(-.1, .1))  # Difficulty
equivTestPaired(14, .8, c(-.1, .1))   # Pace

# experiment 1, 21 subs in 2 cells each
equivTest(21, 21, 7.858, c(-.1, .1)) # Violence?
equivTest(21, 21, -.387, c(-.1, .1)) # Competition
equivTest(21, 21, -1.59, c(-.1, .1))  # Difficulty
equivTest(21, 21, .89, c(-.1, .1))   # Pace


###
# Anderson et al., 2004 and the insufficient pilot test
###

# In which it is argued that Glider Pro and Marathon 2 are equivalent
# 120 subjects across 10 games = 12 participants per game
# Violence difference between Marathon 2 and Glider Pro is r = .842 [.392, .854]
# They report means and mean squared errors, e.g. the average variance around each mean
# We'll have to assume that MSE as the variance within each cell.

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
# Bayes factors: (inverted so they give BF01 instead of BF10)
invertBF(equivTest(12, 12, t.dif, nullInterval=c(-.1, .1)))
invertBF(equivTest(12, 12, t.enj, nullInterval=c(-.1, .1)))
invertBF(equivTest(12, 12, t.act, nullInterval=c(-.1, .1))) # BF here favors the alternative 2.61-to-1
invertBF(equivTest(12, 12, t.fru, nullInterval=c(-.1, .1)))
invertBF(equivTest(12, 12, t.vio, nullInterval=c(-.1, .1)))

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
# # Shit... this is a huge pain in the ass because it's repeated-measures.
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