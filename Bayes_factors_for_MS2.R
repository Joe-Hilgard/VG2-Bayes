# load Christie, Kaye, & Baguley's version of Dienes' BF calculator
source("http://www.lifesci.sussex.ac.uk/home/Zoltan_Dienes/inference/bayesFactorCalc2.R") 
# load Rouder' & Morey 's code for one-tailed Cauchy
source("MetaBF.R")
# Load BayesFactor
require(BayesFactor) 
#Bf(sd, obtained, uniform = FALSE, lower=0, upper=1, meanoftheory=0, sdtheory=1, tails=2)

# Load helper functions
# Note that equivTest() and esciTest() give BF (alt/null) but metaBF and BF03 give BF (null / alt)
equivTest = function(n1, n2, t, nullInterval=NULL, mu=0, rscale=.4) {
  x = rnorm(n1); x = x - mean(x); x = x / sd(x) # center & scale @ 0
  y = rnorm(n2); y = y - mean(y); y = y / sd(y) # center and scale at 0
  # Assume pooled variance of 1 (as we have scaled it as such)
  # Thus, SE depends only on sample size
  SE = sqrt(1/n1 + 1/n2)
  y = y + t*SE # center at t*SE
  return(ttestBF(x, y, nullInterval=nullInterval, paired=F, mu=mu, rscale=rscale))
}
equivTestPaired = function(N, t, nullInterval=NULL, rscale=.4) {
  x = rnorm(N)
  diff = rnorm(N, sd = sqrt(N))
  diff = diff / sd(diff) * sqrt(N) # rescale 
  diff = diff - mean(diff) + t # recenter
  y = x + diff
  return(ttestBF(x, y, nullInterval=nullInterval, paired=T, rscale=rscale))
}
# esciTest() for use when effect size is available but F, t, M(SD) stats aren't.
# It converts effect size r into a t-test, then runs equivTest().
esciTest = function(r, n, rscale=.4, paired=F) {
  se.r = sqrt((1-r^2)/(n-2)) # from http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
  # Maybe it only holds for paired designs?
  t = r / se.r 
  if (paired==F) return(equivTest(n1=ceiling(n/2), n2=floor(n/2), t=t, rscale=rscale))
  else return(equivTestPaired(N=n, t=t, rscale=rscale))
} 
# Why do esciTest(.4, 40, .4) and esciTest(.4, 40, .4, paired=T)
  # return different values? Theoretically/mathematically speaking.
  # Proportion of error variance changes for between vs within designs?
# welch's t
welch.t = function(m1, m2, sd1, sd2, n1, n2) {
  sp = sqrt( ((n1-1)*sd1^2 + (n2-1) * sd2^2) / (n1 + n2 -2) )
  se_diff = sp * sqrt(1/n1 + 1/n2)
  mean_diff = m1-m2
  return(mean_diff/se_diff)
}
pool.sd = function(sds, ns) {
  SSlist = sds^2 %*% (ns-1)
  pool.var = sum(SSlist) / sum(ns-1)
  return(sqrt(pool.var))
}
# for between-subjects tests only
t2d = function(t, n1, n2) {
  d = 2*t / sqrt(n1 + n2 - 2) # equation from http://www.uccs.edu/~lbecker/
  # variance of d from http://stats.stackexchange.com/questions/8487/how-do-you-calculate-confidence-intervals-for-cohens-d
  d.var = ( (n1+n2)/(n1*n2) + d^2 / (2*(n1+n2-2)) ) * ( (n1 + n2)/(n1 + n2 - 2) )
  d.se = sqrt(d.var)
  CI = c(d - 1.96*d.se, d + 1.96*d.se)
  return(list("d"=d, "d.se" = d.se, "d.CI"=CI))
}
# convert effect size r to d
# from http://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf
r2d = function(r, n1, n2) {
  d = 2*r / sqrt(1 - r^2) 
  d.var = ( (n1+n2)/(n1*n2) + d^2 / (2*(n1+n2-2)) ) * ( (n1 + n2)/(n1 + n2 - 2) )
  d.se = sqrt(d.var)
  CI = c(d - 1.96*d.se, d + 1.96*d.se)
  return(list("d"=d, "d.se" = d.se, "d.CI"=CI))
}
# This function just takes the result from the D-C calculator and inverts it
# could also use my equivTest function but that requires the exact per-cell N
BF03 = function(mean, sd, lower, meanoftheory, sdtheory) {
  return(1/Bf(obtained=mean, sd=sd, lower=lower, meanoftheory=meanoftheory, sdtheory=sdtheory)$BayesFactor)
}

# Aggressive affect
# This is complicated by the pre-post design, as means and SDs for pre-post difference scores
  # were not reported.
# Thus, I have to use the effect size r as reported in the manuscript or in personal correspondence
# Some error may be introduced via assumption of equal sample sizes.
# This is the best we get given the numbers made available, people. Share your data!

# Valadez & Ferguson, interaction effect as originally reported (F-value on p. 614)
equivTest(18+15, 14+14+21+18, sqrt(3.11))
1/equivTest(18+15, 14+14+21+18, sqrt(3.11))
d = t2d(sqrt(3.11), 18+15, 14+14+21+18)
1/metaBF(sqrt(3.11),100,half=T,r=.4)
1/BF03(d$d, d$d.se, lower=-1, meanoftheory=.61, sdtheory=.056)


# Valadez & Ferguson, interaction effect violent (RDR either version) vs nonviolent (FIFA)
# F-value received in personal communication 5/10/14
# "For the 2x2 (time x game) interaction (no covariate), I get F (1,97) = 4.884, sig = .029 
#  For the 2x2 (time x game) interaction (gender covariate), I get F (1,97) = 4.884, sig = .029"

equivTest(18+15, 14+14+21+18, sqrt(4.884))
1/equivTest(18+15, 14+14+21+18, sqrt(4.884))
1/metaBF(sqrt(4.884),100,half=T,r=.4)
d = r2d(.22, 18+14+15+21, 14+18)
1/BF03(d$d, d$d.se, lower=-1, meanoftheory=.61, sdtheory=.056)

# Przybylski et al. 2014.
# Numbers retrieved in personal correspondence 5/13/2014
# (There was a very minor copy-paste error)
# Przybylski et al., Study 1
esciTest(.004, 99, .4)
1/esciTest(.004, 99, .4)
t = qt(.515, 98)
equivTest(54, 55, t, rscale=.4)
metaBF(t,100,half=T,r=.4)
d = r2d(.004, 50, 49) # assuming nearly-equal samples
BF03(d$d, d$d.se, lower=-1, meanoftheory=.61, sdtheory=.056)

# Przybylski et al., Study 2
esciTest(-.08, 101, .4) 
1/esciTest(-.08, 101, .4)
t = qt(.5-(1-.41)/2, 101-2)
equivTest(50, 51, t, rscale=.4)
metaBF(t,100,half=T,r=.4)
d = r2d(-.08, 50, 51) # assuming nearly-equal samples
BF03(d$d, d$d.se, lower=-1, meanoftheory=.61, sdtheory=.056)
# Pryzyblski et al., Study 5
  # using standardized regression weight as effect size r
  # reported on p450, left column, last paragraph 
esciTest(.03, 109, .4)
1/esciTest(.03, 109, .4)
t = qt(.5+(1-.74)/2, 109-2)
1/equivTest(54, 55, t)
metaBF(t, 109, half=T, r=.4)
d = r2d(.03, 55, 54)
BF03(d$d, d$d.se, lower=-1, meanoftheory=.61, sdtheory=.056)

# Ivory & Kalyanaraman, 2007
# using F-value reported in Table 1, assuming equal sample sizes
equivTest(60, 60, sqrt(3.83), rscale = .4)
1/equivTest(60, 60, sqrt(3.83), rscale = .4)
1/metaBF(sqrt(3.83), 120, half=T, r=.4)
d = t2d(sqrt(3.83), 60, 60)
1/BF03(d$d, d$d.se, lower=-1, meanoftheory=.61, sdtheory=.056)

# Aggressive behavior
# Elson et al. 2013 noise intensity
equivTest(84/2, 84/2, sqrt(3.28), rscale=.4)
1/equivTest(84/2, 84/2, sqrt(3.28), rscale=.4)
1/metaBF(sqrt(3.28), 84, half=T, r=.4)
d = t2d(sqrt(3.28), 42, 42)
1/BF03(d$d, d$d.se, lower=-1, meanoftheory=.43, sdtheory=.046)
# Elson et al. 2013 noise duration
equivTest(84/2, 84/2, sqrt(.95), rscale=.4)
1/equivTest(84/2, 84/2, sqrt(.95), rscale=.4)
metaBF(sqrt(.95), 84, half=T, r=.4)
d=t2d(sqrt(.95), 42, 42)
1/BF03(d$d, d$d.se, lower=-1, meanoftheory=.43, sdtheory=.046)

# Ferguson et al. 2008
# As I had originally conducted it based on the p-value reported in Table 1:
tval = qt(.55, 50-2) # get t-value according to p=.90, two-tailed
equivTest(26, 24, tval, rscale=.4)
1/equivTest(26, 24, tval, rscale=.4)
metaBF(tval, 50, half=T, r=.4)
d = t2d(tval, 26, 24)
BF03(d$d, d$d.se, lower=-1, meanoftheory=.43, sdtheory=.046)

# As Ferguson has corrected me, saying to use the means and SDs instead,
  # and providing appropriate SPSS output in personal correspondence:
equivTest(26, 24, -.722, rscale=.4)
1/equivTest(26, 24, -.722, rscale=.4)
metaBF(-.722, 50, half=T, r=.4)
d = t2d(-.722, 26, 24)
BF03(d$d, d$d.se, lower=-1, meanoftheory=.43, sdtheory=.046)

# Ferguson & Rueda, 2010, Violent (Hitman, Call of Duty) vs Nonviolent (Madden 07) game
# Using means and SDs retrieved from Table 1
t = welch.t(mean(c(6.03, 6.02)), 5.89, 
            pool.sd(c(1.95, 2.05), c(26, 26)), 2.03,
            26+26, 25)
equivTest(26+26, 25, t, rscale=.4)
1/equivTest(26+26, 25, t, rscale=.4)
metaBF(t, 26+26+25, half=T, r=.4)
d = t2d(t, 26+26, 25)
BF03(d$d, d$d.se, lower=-1, meanoftheory=.43, sdtheory=.046)
# Adachi & Willoughby 2011 exp 1
equivTest(21, 21, 0, rscale=.4)
1/equivTest(21, 21, 0, rscale=.4)
metaBF(0, 42, half=T, r=.4)
d = t2d(0, 21, 21)
BF03(d$d, d$d.se, lower=-1, meanoftheory=.43, sdtheory=.046)
# Adachi & Willoughby 2011 exp 2
# Means and SDs received in personal correspondence from Adachi, 4/29/14
t = welch.t(mean(c(-.776, .904)), mean(c(.785, -.913)),
            pool.sd(c(1.418, 1.347), c(15, 15)), pool.sd(c(1.542, 1.160), c(15, 15)),
            15+15, 15+15)
equivTest(30, 30, t, rscale=.4)
1/equivTest(30, 30, t, rscale=.4)
metaBF(t, 60, half=T, r=.4)
d = t2d(t, 30, 30)
BF03(d$d, d$d.se, lower=-1, meanoftheory=.43, sdtheory=.046)

# Tear & Nielsen, 2014, using values reported in Table 1
# Combine violent & ultraviolent cells, generate t-value of pairwise contrast vs nonviolent
t = welch.t(3.35, mean(3.3, 3.35), 1.97, pool.sd(c(1.52, 1.21), c(40, 40)), 40, 80)
equivTest(40, 80, t, rscale=.4)
1/equivTest(40, 80, t, rscale=.4)
metaBF(t, 120, half=T, r=.4)
d = t2d(t, 40, 80)
BF03(d$d, d$d.se, lower=-1, meanoftheory=.43, sdtheory=.046)

# Aggressive cognition
# Ivory & Kalyaraman, 2007
# using F-value reported in Table 1, assuming equal cell sizes
equivTest(60, 60, sqrt(.17), rscale=.4)
1/equivTest(60, 60, sqrt(.17), rscale=.4)
metaBF(sqrt(.17), 120, half=T, r=.4)
d = t2d(sqrt(.17), 60, 60)
BF03(d$d, d$d.se, lower=-1, meanoftheory=.43, sdtheory=.046)



## BF01 / BF10 for study 2 / table 2 of Elson's CRTT complaint
# inverse so we get BF01
# These are F(1, 80) from Elson's Table 2, plus the reported count of low volume settings
Flist = c(3.28, 1.46, 4.14, .95, .19, 1.74, 2.78, 2.77, 2.17, 16.01, .17, .08, .01,
          10.57) # this last 10.57 is in text but not table, has neg. r.
tlist = sqrt(Flist)
tail(tlist, 1) = -tail(tlist, 1) # turn that one back to negative
# Function to turn these to effect size r
F2R = function(Fstat, N, width=.95, neg=F) {
  r.equiv = sqrt(Fstat/(Fstat + N - 2)) 
  if(neg==T) {r.equiv=-(r.equiv)}
  zScore = 1/2 * log((1+r.equiv)/(1-r.equiv))
  z.se = 1/sqrt(N-3)
  margin = -qnorm((1-width)/2)
  z.low = r.equiv-margin*z.se
  z.hi = r.equiv+margin*z.se
  r.equiv.low = (exp(2*z.low)-1)/(exp(2*z.low)+1)
  r.equiv.hi = (exp(2*z.hi)-1)/(exp(2*z.hi)+1)
  return(c(r.equiv.low, r.equiv, r.equiv.hi))
}
# 21 in each of 4 cells. main effect of violence?
rList = c()
for (f in Flist) {
  r = F2R(f, 84)[2]
  rList = c(rList, r)
  rList = round(rList, 2)
}
tail(rList, 1) = -tail(rList, 1) # flipping that last one
write(rList, "Elson-output-r.txt", ncolumns=1)

# Bayes factors
# BF01
bf01list = c()
for (i in 1:length(tlist)) {
  BF01 = 1/exp(ttest.tstat(t=tlist[i], n1=42, n2=42, rscale = 0.4)[['bf']])
  bf01list = c(bf01list, BF01)
}
write(bf01list, "Elson-output-BF01.txt", ncolumns=1)
#BF02
bf02list = c()
for (i in 1:length(tlist)) {
  BF02 = metaBF(tlist[i], 84, half=T, r=.4)
  bf02list = c(bf02list, BF02)
}
write(bf02list, "Elson-output-BF02.txt", ncolumns=1)
# or use Dienes for BF03:
bf03list = c(); N = 84
for (f in Flist) {
  r.equiv = sqrt(f/(f + N - 2))
  z = atanh(r.equiv)
  if (f == 10.57) z = -z # turn that last one negative, sorry about the kludge.
  z.se = 1/sqrt(N-3)
  bf03 = BF03(z, z.se, lower=-1, meanoftheory=.213171, sdtheory=.026252)
  bf03list = c(bf03list, bf03)
  #zlist = c(zlist, z)
  #zlist.se = c(zlist.se, z.se)
}
write(bf03list, file="Elson-output-BF03.txt", ncolumns=1)


# Other re-analysis of Tear & Nielsen, 2014, using values reported in Table 1
# prosocial stuff irrelevant to current MS, but could bear a mention
# Helper functions to combine violent & ultraviolent cells, generate t-value of pairwise contrast

# helping behavior
help = welch.t(3.2, mean(3.2, 3.02), 1.92, pool.sd(c(1.34, 1.03), c(40, 40)), 40, 80)
# hurting behavior
hurt = welch.t(3.35, mean(3.3, 3.35), 1.97, pool.sd(c(1.52, 1.21), c(40, 40)), 40, 80)
# Charity donation
charity = welch.t(2, mean(2.33, 3.0), 1.94, pool.sd(c(2.03, 2.10), c(40, 40)), 40, 80)
## effect sizes
r.help = t2R(help, 120)$r
r.hurt = t2R(hurt, 120)$r
r.charity = t2R(charity, 120)$r

# Now make BF01s:
1/exp(ttest.tstat(t=help, n1=40, n2=80, rscale = 0.5)[['bf']])
1/exp(ttest.tstat(t=hurt, n1=40, n2=80, rscale = 0.5)[['bf']])
1/exp(ttest.tstat(t=charity, n1=40, n2=80, rscale = 0.5)[['bf']])
# and BF03s:
#BF03(atanh(0), 1/sqrt(120-3), lower=-1, meanoftheory=.213171, sdtheory=.026252) # different mean/sd of theory for prosoc
BF03(atanh(r.hurt), 1/sqrt(120-3), lower=-1, meanoftheory=.213171, sdtheory=.026252)
#BF03(atanh(-.079), 1/sqrt(120-3), lower=-1, meanoftheory=.213171, sdtheory=.026252) # different mean/sd of theory for prosoc
