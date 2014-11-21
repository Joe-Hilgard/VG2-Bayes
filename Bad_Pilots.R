# Function for pilot test w/ null-equivalent region.
equivTest = function(N, t, lo, hi) {
  x = rnorm(N)
  diff = rnorm(N, sd = sqrt(N))
  diff = diff / sd(diff) * sqrt(N) # rescale 
  diff = diff - mean(diff) + t # recenter
  y = x + diff
  return(ttestBF(x, y, nullInterval=c(lo, hi), paired=T, rscale=.5))
}
# Load in BayesFactor package 
require(BayesFactor)

# Arriaga et al., 2008 and the no good rotten pilot test
equivTest(20, .48, -.1, .1)
equivTest(20, .53, -.1, .1)
equivTest(20, .79, -.1, .1)
equivTest(20, .83, -.1, .1)
equivTest(20, .86, -.1, .1)
equivTest(20, .89, -.1, .1)
equivTest(20, 1.14, -.1, .1)
equivTest(20, 1.24, -.1, .1)
equivTest(20, 1.29, -.1, .1)
equivTest(20, 1.32, -.1, .1)
equivTest(20, 1.56, -.1, .1)
equivTest(20, 1.67, -.1, .1)
equivTest(20, 2.27, -.1, .1)
equivTest(20, 2.63, -.1, .1)

## Update the section below:

# Valadez & Ferguson and the no good very bad pilot test
exp(linearReg.R2stat(N=25, p=2, R2=0.5, rscale = 0.353553390593274)[['bf']])

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
dat = data.frame("x" = rep(c("A", "B", "C"), each=15))
dat$y = rnorm(45)
dat$y = dat$y + .3*(dat$x=="A") - .5*(dat$x == "C")
anovaBF(y ~ x, data=dat)

# Elson CRTT
# inverse so we get BF01
Flist = c(3.28,1.46,4.14,.95,.19,1.74,2.78,2.77,2.17,16.01,.17,.08,.01,10.57)
tlist = sqrt(Flist)
# 21 in each of 4 cells. main effect of violence?
sink(file="Elson-output-BF.txt")
for (i in 1:length(tlist)) {
  BF01 = 1/exp(ttest.tstat(t=tlist[i], n1=42, n2=42, rscale = 0.707)[['bf']])
  print(BF01)
}
sink()
# or so we can use Dienes:
sink(file="Elson-output-ESCI.txt")
for (f in Flist) {
  print(F2R_noZ(f, 84))
}
sink()
