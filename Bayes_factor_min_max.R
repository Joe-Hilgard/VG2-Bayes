# Maximum and minimum BF for the null in small samples.
require(BayesFactor)

# Prepare plot space
pdf('BFminmax.pdf', width=8, height=4) # attempting to scale approximately to A4 dimensions #width=12,height=24) 
par(mfrow=c(1,2),cex=1,mar=c(4,4,.5,1),mgp=c(2.2,1,0))

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
# Another line for p = .10, just for contours
BF01_twoSamp_med = NA
for (i in 1:length(nVector)) {
  nCell = nVector[i]/2
  tVal = qt(.95, nVector[i]-2)
  BF01_twoSamp_med[i] = 
    1/exp(ttest.tstat(t=tVal, n1=nCell, n2=nCell, rscale = 0.5)[['bf']])  
}

plot(nVector, BF01_twoSamp_min, typ='l', log="y", ylim=c(3^-1,6),
     xlab = "Sample size (Total N across 2 cells)",
     ylab = "Bayes factor (null / alternative)", 
     yaxt = 'n',
     lty=3, lwd=2)
axis(2, at=c(10^-1, 6^-1, 3^-1, 1, 3, 6, 10)
     , labels=c("1:10", "1:6", "1:3", "1:1", "3:1", "6:1", "10:1"),
)
abline(h=1, col='grey')
lines(nVector, BF01_twoSamp_max, typ='l')
lines(nVector, BF01_twoSamp_med, typ='l', lty=5)
mtext(side=3,adj=.05,cex=1.5,'A.',line=-1.5)

# One-sample t-test
nVector1 = seq(3, 100, 1)
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
# BF at p = .10 just for contour
BF01_oneSamp_med = NA
for (i in 1:length(nVector1)) {
  tVal = qt(.95, nVector1[i]-1)
  BF01_oneSamp_med[i] = 
    1/exp(ttest.tstat(t=tVal, n1=nVector1[i], rscale = 0.5)[['bf']])  
}


plot(nVector1, BF01_oneSamp_max, typ='l', log="y", ylim=c(3^-1,6),
     xlab = "Sample size (Total N, repeated measures)",
     ylab = "Bayes factor (null / alternative)", 
     yaxt = 'n')
axis(2, at=c(10^-1, 6^-1, 3^-1, 1, 3, 6, 10)
     , labels=c("1:10", "1:6", "1:3", "1:1", "3:1", "6:1", "10:1"),
)
abline(h=1, col='grey')
lines(nVector1, BF01_oneSamp_min, typ='l', lty=3, lwd=2)
lines(nVector1, BF01_oneSamp_med, typ='l', lty=5)
mtext(side=3,adj=.05,cex=1.5,'B.',line=-1.5)

# end
dev.off()
