# Function for pilot test w/ null-equivalent region.
equivTestPaired = function(N, t, lo, hi) {
  x = rnorm(N)
  diff = rnorm(N, sd = sqrt(N))
  diff = diff / sd(diff) * sqrt(N) # rescale 
  diff = diff - mean(diff) + t # recenter
  y = x + diff
  return(ttestBF(x, y, nullInterval=c(lo, hi), paired=T, rscale=.5))
}
equivTest = function(n1, n2, t, lo, hi) {
  x = rnorm(n1); x = x - mean(x); x = x / sd(x) # center & scale @ 0
  y = rnorm(n2)
  y = y/sd(y); y = y - mean(y) # center and scale at 0
  SE = sqrt(var(x)/n1 + var(y)/n2)
  y = y + t*SE # center at t*SE
  return(ttestBF(x, y, nullInterval=c(lo, hi), paired=F, rscale=.5))
}
# Load in BayesFactor package 
require(BayesFactor)

# Arriaga et al., 2008 and the no good rotten pilot test
equivTestPaired(20, .48, -.1, .1)
equivTestPaired(20, .53, -.1, .1)
equivTestPaired(20, .79, -.1, .1)
equivTestPaired(20, .83, -.1, .1)
equivTestPaired(20, .86, -.1, .1)
equivTestPaired(20, .89, -.1, .1)
equivTestPaired(20, 1.14, -.1, .1)
equivTestPaired(20, 1.24, -.1, .1)
equivTestPaired(20, 1.29, -.1, .1)
equivTestPaired(20, 1.32, -.1, .1)
equivTestPaired(20, 1.56, -.1, .1)
equivTestPaired(20, 1.67, -.1, .1)
equivTestPaired(20, 2.27, -.1, .1)
equivTestPaired(20, 2.63, -.1, .1)

## Update the section below:

# Valadez & Ferguson and the no good very bad pilot test
# what was I doing here?
# I want to make t-tests for each pairwise comparison
# To do that, I need the t-values of the pairwise comparisons...

# Cell sizes: RDR-v, 15; RDR-nv, 10; FIFA, 15.
# RDR "violent" vs RDR "nonviolent"
equivTest(15, 10, 1.82, -.1, .1) # Difficulty
equivTest(15, 10, 1.31, -.1, .1) # Pace
equivTest(15, 10, 3.00, -.1, .1) # Competitveness
# RDR "violent" vs FIFA
equivTest(15, 10, 1.47, -.1, .1) # Difficulty
equivTest(15, 10, 2.00, -.1, .1) # Pace
equivTest(15, 10, .047, -.1, .1) # Competitveness
# RDR "nonviolent" vs FIFA 
equivTest(15, 10, 3.45, -.1, .1) # Difficulty
equivTest(15, 10, 3.43, -.1, .1) # Pace
equivTest(15, 10, 3.00, -.1, .1) # Competitveness

# Adachi & Willoughby -- Are they matched?
# pilot 1
# inverse so that we get BF01 not BF10
1/exp(ttest.tstat(t=4.39, n1=14, rscale = 0.707)[['bf']])
1/exp(ttest.tstat(t=-.46, n1=14, rscale = 0.707)[['bf']])
1/exp(ttest.tstat(t=.59, n1=14, rscale = 0.707)[['bf']])
1/exp(ttest.tstat(t=0.8, n1=14, rscale = 0.707)[['bf']])
# experiment 1, inverse so we get BF01
1/exp(ttest.tstat(t=7.858, n1=21, n2=21, rscale = 0.707)[['bf']])
1/exp(ttest.tstat(t=-.387, n1=21, n2=21, rscale = 0.707)[['bf']])
1/exp(ttest.tstat(t=-1.59, n1=21, n2=21, rscale = 0.707)[['bf']])
1/exp(ttest.tstat(t=.89, n1=21, n2=21, rscale = 0.707)[['bf']])

# experiment 2 it would be necessary to have the raw data to perform anovaBF()
# No, I could simulate this still... I need to pull from a normal,
  # Then standardize & center to appropriate means, that's all.
# Doing the full model-comparison BF could get messy, but would be very interesting.
dat = data.frame("x" = rep(c("A", "B", "C"), each=15))
dat$y = rnorm(45)
dat$y = dat$y + .3*(dat$x=="A") - .5*(dat$x == "C")
anovaBF(y ~ x, data=dat)

