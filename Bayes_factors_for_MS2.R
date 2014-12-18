# load Christie, Kaye, & Baguley's version of Dienes' BF calculator
source("http://www.lifesci.sussex.ac.uk/home/Zoltan_Dienes/inference/bayesFactorCalc2.R") 
# Load BayesFactor
require(BayesFactor) 
#Bf(sd, obtained, uniform = FALSE, lower=0, upper=1, meanoftheory=0, sdtheory=1, tails=2)

# Load helper functions
equivTest = function(n1, n2, t, nullInterval=NULL, mu=0, rscale=.5) {
  x = rnorm(n1); x = x - mean(x); x = x / sd(x) # center & scale @ 0
  y = rnorm(n2); y = y - mean(y); y = y / sd(y)  # center and scale at 0
  SE = sqrt(var(x)/n1 + var(y)/n2)
  y = y + t*SE # center y at t*SE
  return(ttestBF(x, y, nullInterval=nullInterval, paired=F, mu=mu, rscale=rscale))
}
equivTestPaired = function(N, t, nullInterval=NULL, rscale=.5) {
  x = rnorm(N)
  diff = rnorm(N, sd = sqrt(N))
  diff = diff / sd(diff) * sqrt(N) # rescale 
  diff = diff - mean(diff) + t # recenter
  y = x + diff
  return(ttestBF(x, y, nullInterval=nullInterval, paired=T, rscale=rscale))
}
esciTest = function(r, n, rscale, paired=F) {
  se.r = sqrt((1-r^2)/(n-2)) # from http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
  # Maybe it only holds for paired designs?
  t = r / se.r 
  if (paired==F) return(equivTest(n1=ceiling(n/2), n2=floor(n/2), t=t, rscale=rscale))
  else return(equivTestPaired(N=n, t=t, rscale=rscale))
} 
invertBF = function(model) {
  return(1/exp(model@bayesFactor[['bf']]))
}
# Why do esciTest(.4, 40, .4) and esciTest(.4, 40, .4, paired=T)
  # return different values? Theoretically/mathematically speaking.
r2t = function(r, n) {
  se.r = sqrt((1-r^2)/(n-2)) # from http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
  t = r / se.r
  return(t)
}

# This function just takes the result from the D-C calculator and inverts it
BF02 = function(mean, sd, lower, meanoftheory, sdtheory) {
  return(1/Bf(obtained=mean, sd=sd, lower=lower, meanoftheory=meanoftheory, sdtheory=sdtheory)$BayesFactor)
}
# could also use my equivTest function but that requires the exact per-cell N
# List of BF01's from JZS Bayes, scale=.4
# Remember slotNames() to see items w/in S4 environment

# Aggressive affect
# Valadez & Ferguson, interaction effect violent (RDR either version) vs nonviolent (FIFA)
esciTest(.22, 100, .4) # MAY BE INVALIDATED BY UNEQUAL CELL SIZES?
invertBF(esciTest(.22, 100, .4))
# Valdez & Ferguson, interaction effect as originally reported
esciTest(.17, 100, .4)
invertBF(esciTest(.17, 100, .4))
# Przybylski et al., Study 1
esciTest(.004, 99, .4)
invertBF(esciTest(.004, 99, .4))
# Przybylski et al., Study 2
esciTest(.08, 101, .4) # or is it -.08?
invertBF(esciTest(.08, 101, .4))
# Pryzyblski et al., Study 5
esciTest(.03, 109, .4)
invertBF(esciTest(.03, 109, .4))
# Ivory & Kalyaraman, 2007
esciTest(.13, 120, .4)
invertBF(esciTest(.13, 120, .4))
# Aggressive behavior
# Elson et al. 2013 noise intensity
esciTest(.2, 84, .4)
invertBF(esciTest(.2, 84, .4))
# Elson et al. 2013 noise duration
esciTest(.11, 84, .4)
invertBF(esciTest(.11, 84, .4))
# Ferguson et al. 2008
esciTest(.02, 50, .4)
invertBF(esciTest(.02, 50, .4))
# Ferguson & Rueda, Violent vs Nonviolent game
esciTest(.01, 77, .4)
invertBF(esciTest(.01, 77, .4))
# Adachi & Willoughby 2011 exp 1
esciTest(0, 42, .4)
invertBF(esciTest(0, 42, .4))
# Adachi & Willoughby 2011 exp 2
esciTest(.03, 60, .4)
invertBF(esciTest(.03, 60, .4))

# Aggressive cognition
# Ivory & Kalyaraman, 2007
esciTest(-.08, 120, .4)
invertBF(esciTest(-.08, 120, .4))


# List of BF02's from Dienes-Christie calculator
  
# All effect sizes are in r-to-z units s.t. normal distribution applies 
# Aggressive affect
# Valdez & Ferguson, interaction effect violent (RDR either version) vs nonviolent (FIFA)
BF02(sd=.101535, .223656, lower=-1, meanoftheory=.298566, sdtheory=.01996) # where am I getting all these sig digs from?
#BF02(sd=.10, .22, lower=-1, meanoftheory=.30, sdtheory=.02) # rounding matters very little
# Valdez & Ferguson, interaction effect as originally reported
BF02(sd=.101535, .171667, lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Przybylski et al., Study 1
BF02(sd=.102062, 0, lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Przybylski et al., Study 2
BF02(sd=.101015, .080171, lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Przybylski et al., Study 5
BF02(.03, 1/sqrt(109-3), lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Ivory & Kalyraman, 2007
BF02(sd=.09245, mean=.13074, lower=-1, meanoftheory=.298566, sdtheory=.01996)

# Aggressive behavior
# Elson et al. 2013 noise intensity
BF02(sd=.111111, .202733, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Elson et al. 2013 noise duration
BF02(sd=.111111, .110447, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Ferguson et al. 2008
BF02(sd=.145865, .020003, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Ferguson & Rueda, Violent vs Nonviolent game
BF02(sd=.116248, .01, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Adachi & Willoughby 2011 exp 1
BF02(sd=.160128, .0, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Adachi & Willoughby 2011 exp 2
BF02(sd=.132453, .030009, lower=-1, meanoftheory=.213171, sdtheory=.026252)

# Aggressive cognition
# Ivory & Kalyaraman, 2007
BF02(sd=.09245, -.08017, lower=-1, meanoftheory=.223656, sdtheory=.018621)


## BF01 / BF10 for study 2 / table 2 of Elson's CRTT complaint
# inverse so we get BF01
# These are F(1, 80) from Elson's Table 2, plus the reported count of low volume settings
Flist = c(3.28, 1.46, 4.14, .95, .19, 1.74, 2.78, 2.77, 2.17, 16.01, .17, .08, .01,
          10.57) # this last 10.57 is in text but not table, has neg. r.
tlist = sqrt(Flist)
# Function to turn these to effect size r
F2R = function(Fstat, N, width=.95, neg=F) {
  r.equiv = sqrt(Fstat/(Fstat + N - 2)) # is this where the loss of fidelity happens?
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
rList[length(rList)] = -rList[length(rList)] # flipping that last one
write(rList, "Elson-output-r.txt", ncolumns=1)

bf01list = c()
for (i in 1:length(tlist)) {
  BF01 = 1/exp(ttest.tstat(t=tlist[i], n1=42, n2=42, rscale = 0.4)[['bf']])
  bf01list = c(bf01list, BF01)
}
write(bf01list, "Elson-output-BF01.txt", ncolumns=1)
# or use Dienes for BF02:
bf02list = c(); N = 84
for (f in Flist) {
  r.equiv = sqrt(f/(f + N - 2))
  z = atanh(r.equiv)
  if (f == 10.57) z = -z # turn that last one negative, sorry about the kludge.
  z.se = 1/sqrt(N-3)
  bf02 = BF02(z, z.se, lower=-1, meanoftheory=.213171, sdtheory=.026252)
  bf02list = c(bf02list, bf02)
  #zlist = c(zlist, z)
  #zlist.se = c(zlist.se, z.se)
}
write(bf02list, file="Elson-output-BF02.txt", ncolumns=1)

BF02(sd=.111111, .202733, lower=-1, meanoftheory=.213171, sdtheory=.026252)


# sink(file="Elson-output-ESCI.txt")
# for (f in Flist) {
#   print(F2R_noZ(f, 84))
# }
# sink()

# Re-analysis of Tear & Nielsen, 2014, using values reported in Table 1
# Helper functions to combine violent & ultraviolent cells, generate t-value of pairwise contrast
pool.sd = function(sds, ns) {
  num = sum((ns-1)*sds^2)
  denom = sum(ns-1)
  pool.var = num/denom
  pool.sd = sqrt(pool.var)
  return(pool.sd)
}
welch.t = function(m1, m2, sd1, sd2, n1, n2) {
  se_diff = sqrt(sd1^2/n1 + sd2^2/n2)
  mean_diff = m1-m2
  return(mean_diff/se_diff)
}
t2R = function(tstat, N, digits=2) {
  neg = tstat<0
  Fstat = tstat^2
  r.equiv = sqrt(Fstat/(Fstat + N - 2))
  if(neg==T) {r.equiv=-(r.equiv)}
  zScore = 1/2 * log((1+r.equiv)/(1-r.equiv))
  z.se = 1/sqrt(N-3)
  z.low = r.equiv-1.96*z.se
  z.hi = r.equiv+1.96*z.se
  r.equiv.low = (exp(2*z.low)-1)/(exp(2*z.low)+1)
  r.equiv.hi = (exp(2*z.hi)-1)/(exp(2*z.hi)+1)
  print(paste("Point estimate:", r.equiv))
  print(paste("95% CI: [", r.equiv.low, ", ", r.equiv.hi, "]", sep=""))
  return(c(r.equiv.low, r.equiv, r.equiv.hi))
}
n = 120/3
# helping behavior
help = welch.t(3.2, mean(3.2, 3.02), 1.92, pool.sd(c(1.34, 1.03), c(40, 40)), 40, 80)
# hurting behavior
hurt = welch.t(3.35, mean(3.3, 3.35), 1.97, pool.sd(c(1.52, 1.21), c(40, 40)), 40, 80)
# Charity donation
charity = welch.t(2, mean(2.33, 3.0), 1.94, pool.sd(c(2.03, 2.10), c(40, 40)), 40, 80)
## effect sizes
t2R(help, 120)
t2R(hurt, 120)
t2R(charity, 120)

# Now make BF01s:
1/exp(ttest.tstat(t=help, n1=40, n2=80, rscale = 0.5)[['bf']])
1/exp(ttest.tstat(t=hurt, n1=40, n2=80, rscale = 0.5)[['bf']])
1/exp(ttest.tstat(t=charity, n1=40, n2=80, rscale = 0.5)[['bf']])
# and BF02s:
#BF02(atanh(0), 1/sqrt(120-3), lower=-1, meanoftheory=.213171, sdtheory=.026252) # different mean/sd of theory for prosoc
BF02(atanh(.013), 1/sqrt(120-3), lower=-1, meanoftheory=.213171, sdtheory=.026252)
#BF02(atanh(-.079), 1/sqrt(120-3), lower=-1, meanoftheory=.213171, sdtheory=.026252) # different mean/sd of theory for prosoc
