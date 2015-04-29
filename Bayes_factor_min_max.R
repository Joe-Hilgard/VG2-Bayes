# Maximum and minimum BF for the null in small samples.
require(BayesFactor)

# Independent groups t-test
# BF is maximal at t = 0
nVector = seq(4, 100, 2)
# BF is maximal at t = 0
BF01_twoSamp_max = NA
for (i in 1:length(nVector)) {
  nCell = nVector[i]/2
  BF01_twoSamp_max[i] = 
    1/exp(ttest.tstat(t=0, n1=nCell, n2=nCell, rscale = 0.5)[['bf']])
}
# BF is minimal at p = .05
BF01_twoSamp_min = NA
for (i in 1:length(nVector)) {
  nCell = nVector[i]/2
  tVal = qt(.975, nVector[i]-2)
  BF01_twoSamp_min[i] = 
    1/exp(ttest.tstat(t=tVal, n1=nCell, n2=nCell, rscale = 0.5)[['bf']])  
}

plot(nVector, BF01_twoSamp_min, typ='l', log="y", ylim=c(0.5,10))
lines(nVector, BF01_twoSamp_max, typ='l')

# One-sample t-test
nVector1 = seq(2, 100, 1)
# BF is maximal at t = 0
BF01_oneSamp_max = NA
for (i in 1:length(nVector1)) {
  BF01_oneSamp_max[i] = 
    1/exp(ttest.tstat(t=0, n1=nVector1[i], rscale = 0.5)[['bf']])  
}
# BF is minimal at p = .05
BF01_oneSamp_min = NA
for (i in 1:length(nVector1)) {
  tVal = qt(.975, nVector1[i]-1)
  BF01_oneSamp_min[i] = 
    1/exp(ttest.tstat(t=tVal, n1=nVector1[i], rscale = 0.5)[['bf']])  
}

plot(nVector1, BF01_oneSamp_min, typ='l', log="y", ylim=c(0.5,10))
lines(nVector1, BF01_oneSamp_max, typ='l')
