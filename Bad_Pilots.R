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
# Load in BayesFactor package 
require(BayesFactor)

# Arriaga et al., 2008 and the no good rotten pilot test
equivTestPaired(20, .48, c(-.1, .1))
equivTestPaired(20, .53, c(-.1, .1))
equivTestPaired(20, .79, c(-.1, .1))
equivTestPaired(20, .83, c(-.1, .1))
equivTestPaired(20, .86, c(-.1, .1))
equivTestPaired(20, .89, c(-.1, .1))
equivTestPaired(20, 1.14, c(-.1, .1))
equivTestPaired(20, 1.24, c(-.1, .1))
equivTestPaired(20, 1.29, c(-.1, .1))
equivTestPaired(20, 1.32, c(-.1, .1))
equivTestPaired(20, 1.56, c(-.1, .1))
equivTestPaired(20, 1.67, c(-.1, .1))
equivTestPaired(20, 2.27, c(-.1, .1))
equivTestPaired(20, 2.63, c(-.1, .1))

## Update the section below:

# Valadez & Ferguson and the no good very bad pilot test
# what was I doing here?
# I want to make t-tests for each pairwise comparison
# To do that, I need the t-values of the pairwise comparisons...

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

# Adachi & Willoughby -- Are they matched?
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

# pilot study 2, N = 19 within 
dat = data.frame("Condition" = rep(c("C_NV", "C_V", "LC_NV", "LC_V"), each=4),
                 "Violence" = rep(c(0, 1, 0, 1), each=4),
                 "Compete" = rep(c(1, 1, 0, 0), each=4),
                 "Outcome" = rep(c("Violence", "Competitiveness", "Difficulty", "Pace"),4),
                 "Mean" = c(1.52, 5.86, 3.63, 4.74,
                            5.37, 6.32, 4.47, 5.42,
                            1, 1.36, 4.68, 4.57,
                            6.42, 3.18, 4.42, 5.11),
                 "SD" = c(.7, .76, 1.38, 1.05,
                          .9, .58, 1.58, .90,
                          0, .44, 1.16, 1.26,
                          .84, 1.67, 1.74, 1.29)
)
# generate data
# Shit... this is a huge pain in the ass because it's repeated-measures.
N = 19
genData = c()
for (i in dat$Condition) {
  for (j in dat$Outcome) {
    activeRow = dat[dat$Condition == i & dat$Outcome == j,]
    cellMean = activeRow$Mean
    cellSD = activeRow$SD
    cellDat = rnorm(N, cellMean, cellSD)
    cellDat = cellDat - mean(cellDat) + cellMean # center precisely
    if (cellSD > 0) cellDat = cellDat / sd(cellDat) # scale precisely
    newData = data.frame("Condition" = i, "Outcome" = j, "Value" = cellDat)
    rbind(genData, newData)
  }
}

# experiment 2 it would be necessary to have the raw data to perform anovaBF()
# No, I could simulate this still... I need to pull from a normal,
  # Then standardize & center to appropriate means, that's all.
# Doing the full model-comparison BF could get messy, but would be very interesting.
dat = data.frame("x" = rep(c("A", "B", "C"), each=15))
dat$y = rnorm(45)
dat$y = dat$y + .3*(dat$x=="A") - .5*(dat$x == "C")
anovaBF(y ~ x, data=dat)

