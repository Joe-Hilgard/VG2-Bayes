#library('msm')
#one and two-sample implemented

setwd("~/VG2-Bayes/")
N=40  #two-sample
bigN=400 # two-sample
# N=20 #one-sample

z=seq(-2, 2,.005) # over what interval do we plot?
I=length(z)

# Generate estimates
#Cauchy
s=.4 # set scale
f=dcauchy(z,0,s) # get density function over interval Z for cauchy w/ scale s

# Epred gives prediction of the null model:
EDens=function(obs,N) dt(sqrt(N)/2*obs,(N-2))
# EDens=function(obs,N) dt(sqrt(N)*obs,N-1) # pdf of expected data (will treat N as constant) 
# effect size above is t-value * sqrt(N) so that it is scaled to standard dev. not standard error.
Epred=EDens(z,N) # run this over z for the expected data given H0
bigEpred=EDens(z,bigN)

# Ppred gives prediction of the point-alternative effect model:
PDens = function(obs, N) dt(sqrt(N)/2*obs, N-2, ncp = sqrt(N)/2*.2)
Ppred = PDens(z, N)
bigPpred = PDens(z, bigN)
# Cpred gives prediction of the Cauchy-distributed effect model:
Cintgrand=function(delta,obs,N) dt(sqrt(N/2)*obs,2*(N-1),sqrt(N/2)*delta)*dcauchy(delta,0,s)
#Cintgrand=function(delta,obs,N) {
#  dt(sqrt(N)*obs, df = N-1, ncp = sqrt(N) * delta) * dcauchy(delta, 0 ,s) # convolved t with cauchy
#  # This is product of prob of observed effect size at some effect size delta TIMES
#  # the prior probability of that effect size delta given the model (Cauchy with scale s)
#}
# integrate Cintgrand over obs=[-20, 20]
CDens=function(obs,N) integrate(Cintgrand,lower=-20,upper=20,obs=obs,N=N)$value 
# Now integrate over every possible effect size on z
Cpred=1:I; bigCpred=1:I
for (i in 1:I) Cpred[i]=CDens(z[i],N) # but obs is only two points?
for (i in 1:I) bigCpred[i]=CDens(z[i],bigN)

# Npred gives prediction of Anderson's point hypothesis
Nintgrand = function(delta, obs, N) { # function of delta, treating obs and N as constants
  dt(sqrt(N/2)*obs, df = 2*(N-1), ncp = sqrt(N/2)*delta) * dnorm(delta, mean = 0.43, sd = .1)
}
NDens = function(obs, N) integrate(Nintgrand, lower=-20, upper=20, obs=obs, N=N, 
                                   subdivisions=1000, rel.tol = .Machine$double.eps^0.5)$value
Npred = 1:I; bigNpred = 1:I
for (i in 1:I) Npred[i] = NDens(z[i], N)
for (i in 1:I) bigNpred[i] = NDens(z[i], bigN)

# generate data for figure 2A, 2B
BF = function(obs, N) EDens(obs, N)/CDens(obs, N)
effect = seq(-1, 1, .005)
BFList = 1:length(effect); bigBFList = BFList
for (i in 1:length(effect)) BFList[i] = EDens(effect[i],N)/CDens(effect[i],N)
for (i in 1:length(effect)) bigBFList[i] = EDens(effect[i],bigN)/CDens(effect[i],bigN)
############

# Prepare plot space
pdf('BFFigure.pdf', width=8, height=11) # attempting to scale approximately to A4 dimensions #width=12,height=24) 
par(mfrow=c(3,2),cex=1,mar=c(4,4,.5,1),mgp=c(2.2,1,0)) # 

# Plotting, row 1
# Priors:
plot(z,f,typ='n', # type 'n' for no plotting, just scales
     xlab=expression(paste("True Effect Size, ",delta)), xlim = c(-1, 1),
     ylab="Density",ylim=c(0,1)) 
#lines(z,f,col="darkgreen",lty=2,lwd=2) # draw Cauchy density function
arrows(0,0,0,1) # Draw Dirac for null hypothesis
arrows(.2,0,.2,1) # Draw Dirac for point-alternative hypothesis
# arrows(.21, 0, .21, 1, col="red") # Draw Dirac for Anderson's hypothesis (not truly Dirac but might as well be)
# text labels
# Add text for red dirac arrow?
text(-.2,.8,expression(paste("Model ",M[0])),adj=1)
text(.3,.25,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
mtext(side=3,adj=.05,cex=1.5,'A.',line=-1.2)

# Probabilities given data:
plot(z,Epred,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
lines(z,Ppred,lty=3,lwd=2,col='darkred') # data given point-alternative
obs = .4
abline(v=obs,col='grey')
points(obs,EDens(obs,N),cex=1.3,pch=19)
points(obs,PDens(obs,N),cex=1.3,pch=19)
PDens(obs,N)/EDens(obs,N)
#points(obs,CDens(obs[1],N),cex=1.3,pch=19,col='darkgreen')
#points(obs,NDens(obs, N), cex=1.3, pch=19, col='darkred')
#lines(z, Cpred,lty=2,lwd=2,col='darkgreen') # data given JZS prior
#lines(z, Npred, lty=3, lwd=2, col='darkred') # data given Anderson's hypothesis
text(-.35,.3,expression(paste("Model ",M[0])),adj=1)
text(.9,.15,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
# points(obs[2],EDens(obs[2],N),cex=1.3,pch=21,bg='white',lwd=2)
# points(obs[2],CDens(obs[2],N),cex=1.3,pch=21,col='darkgreen',bg='white',lwd=2)
mtext(side=3,adj=.05,cex=1.5,'B.',line=-1.2)

# Plotting, row 2
# Priors
plot(z,f,typ='n', # type 'n' for no plotting, just scales
     xlab=expression(paste("True Effect Size, ",delta)), xlim = c(-1, 1),
     ylab="Density",ylim=c(0,1)) 
lines(z,f,col="darkgreen",lty=2,lwd=2) # draw Cauchy density function
arrows(0,0,0,1) # Draw Dirac for null hypothesis
# text labels
# Add text for red dirac arrow?
text(-.2,.9,expression(paste("Model ",M[0])),adj=1)
text(.4,.5,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
mtext(side=3,adj=.05,cex=1.5,'C.',line=-1.2)
# Probabilities given data:
plot(z,Epred,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
lines(z,Cpred,lty=3,lwd=2,col='darkred') # data given JZS-prior alternative
obs = .4
abline(v=obs,col='grey')
points(obs,EDens(obs,N),cex=1.3,pch=19)
points(obs,CDens(obs,N),cex=1.3,pch=19)
text(-.4,.2,expression(paste("Model ",M[0])),adj=1)
text(.4,.1,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
mtext(side=3,adj=.05,cex=1.5,'D.',line=-1.2)
EDens(obs,N)/CDens(obs,N)

# Plotting, row 3
# H0 vs Craig meta with big N
plot(z,bigEpred,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
lines(z,bigPpred,lty=3,lwd=2,col='darkred') # data given point-alternative
obs = .4
abline(v=obs,col='grey')
points(obs,EDens(obs,bigN),cex=1.3,pch=19)
points(obs,PDens(obs,bigN),cex=1.3,pch=19)
PDens(obs,bigN)/EDens(obs,bigN)
text(-.35,.3,expression(paste("Model ",M[0])),adj=1)
text(.9,.15,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
mtext(side=3,adj=.05,cex=1.5,'E.',line=-1.2)
# H0 vs vague alternative with big N
plot(z,bigEpred,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
lines(z,bigCpred,lty=3,lwd=2,col='darkred') # data given JZS-prior alternative
obs = .4
abline(v=obs,col='grey')
points(obs,EDens(obs,bigN),cex=1.3,pch=19)
points(obs,CDens(obs,bigN),cex=1.3,pch=19)
text(-.4,.2,expression(paste("Model ",M[0])),adj=1)
text(.4,.1,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
mtext(side=3,adj=.05,cex=1.5,'F.',line=-1.2)
EDens(obs,bigN)/CDens(obs,bigN)

dev.off()

# Figure 2:
# BF10 vs observed effect size:
plot(x=effect, xlab = expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),
     y=BFList, log="y", ylim=c(10, .00001), ylab="Bayes factor (alternative/null)",
     typ='l',
     yaxt='n'
     )
axis(2, at=c(.0001, .001, .01, .1, 1, 10), 
     labels=c(".0001", ".001", ".01", ".1", "1", "10"),
     )
lines(x=effect, y=bigBFList, lty=2, col="darkred")
#mtext(side=3,adj=.05,cex=1.5,'A.',line=-1.2)

dev.off()

# #################################################
# myRate=4
# myShape=3
# mySpec=function(z,rate,shape) .5*dgamma(z,shape=shape,rate=rate)+.5*dgamma(-z,shape=shape,rate=rate)
# f=mySpec(z,rate=myRate,shape=myShape)
# plot(z,f,typ='n',xlab=expression(paste("True Effect Size, ",delta)),ylab="Density",ylim=c(0,1))
# lines(z,f,col="darkgreen",lty=2,lwd=2)
# arrows(0,0,0,1)
# text(-.2,.9,expression(paste("Model ",M[0])),adj=1)
# text(1,.4,expression(paste("Model ",M[b])),col="darkgreen",adj=0)
# mtext(side=3,adj=.05,cex=1.7,'C.',line=-1.1)
# 
# Uintgrand=function(delta,obs,N)
#   #dt(sqrt(N/2)*obs,2*(N-1),sqrt(N/2)*delta)*mySpec(delta,rate=myRate,shape=myShape)
#   dt(sqrt(N)*obs,N-1,sqrt(N)*delta)*mySpec(delta,rate=myRate,shape=myShape)
# 
# UDens=function(obs,N) integrate(Uintgrand,lower=-6,upper=6,obs=obs,N=N)$value
# 
# Upred=1:I
# for (i in 1:I) Upred[i]=UDens(z[i],N)
# plot(z,Epred,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
# abline(v=obs,col='grey')
# lines(z,Upred,lty=2,lwd=2,col='darkgreen')
# text(-.35,.3,expression(paste("Model ",M[0])),adj=1)
# text(1.05,.08,expression(paste("Model ",M[b])),col="darkgreen",adj=0)
# points(obs[1],EDens(obs[1],N),cex=1.3,pch=19)
# points(obs[1],UDens(obs[1],N),cex=1.3,pch=19,col='darkgreen')
# points(obs[2],EDens(obs[2],N),cex=1.3,pch=21,bg='white',lwd=2)
# points(obs[2],UDens(obs[2],N),cex=1.3,pch=21,col='darkgreen',bg='white',lwd=2)
# mtext(side=3,adj=.05,cex=1.7,'D.',line=-1.1)
# 
# ####################
# 
# f=dtnorm(z,0,1,0,Inf)
# plot(z,f,typ='n',xlab=expression(paste("True Effect Size, ",delta)),ylab="Density",ylim=c(0,1))
# lines(z,f,col="darkgreen",lty=2,lwd=2)
# #lines(-z,f,typ='l',col="darkred",lwd=2,lty=2)
# arrows(0,0,0,1)
# text(-.2,.9,expression(paste("Model ",M[0])),adj=1)
# text(1.1,.5,expression(paste("Model ",M[p])),col="darkgreen",adj=0)
# mtext(side=3,adj=0,cex=1.5,'E.')
# 
# WMintgrand=function(delta,obs,N)
#   #dt(sqrt(N/2)*obs,2*(N-1),sqrt(N/2)*delta)*dtnorm(delta,0,1,0,Inf)
#   dt(sqrt(N)*obs,(N-1),sqrt(N)*delta)*dtnorm(delta,0,1,0,Inf)
# 
# WMDens=function(obs,N) integrate(WMintgrand,lower=0,upper=10,obs=obs,N=N)$value
# 
# WMpred=1:I
# for (i in 1:I) WMpred[i]=WMDens(z[i],N)
# plot(z,Epred,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
# abline(v=obs,col='grey')
# lines(z,WMpred,lty=2,lwd=2,col='darkgreen')
# lines(-z,WMpred,lty=3,lwd=2,col='darkred')
# text(-.35,.3,expression(paste("Model ",M[0])),adj=1)
# text(1.0,.13,expression(paste("Model ",M[p])),col="darkgreen",adj=0)
# text(-1.0,.13,expression(paste("Model ",M[n])),col="darkred",adj=1)
# points(obs[1],EDens(obs[1],N),cex=1.3,pch=19)
# points(obs[1],WMDens(obs[1],N),cex=1.3,pch=19,col='darkgreen')
# points(obs[2],EDens(obs[2],N),cex=1.3,pch=21,bg='white',lwd=2)
# points(obs[2],WMDens(obs[2],N),cex=1.3,pch=21,col='darkgreen',bg='white',lwd=2)
# points(obs[1],WMDens(-obs[1],N),cex=1.3,pch=19,col='darkred')
# points(obs[2],WMDens(-obs[2],N),cex=1.3,pch=21,col='darkred',bg='white',lwd=2)
# mtext(side=3,adj=.05,cex=1.7,'F.',line=-1.1)
# 
# 
# ######################
# plot(z,f,typ='n',xlab=expression(paste("True Effect Size, ",delta)),ylab="Density",ylim=c(0,1))
# arrows(0,0,0,1)
# arrows(.2,0,.2,1,lty=2,col='darkgreen',lwd=2)
# text(-.2,.9,expression(paste("Model ",M[0])),adj=1)
# text(.4,.9,expression(paste("Model ",M[d])),col="darkgreen",adj=0)
# mtext(side=3,adj=0,cex=1.5,'E.')
# mtext(side=3,adj=.05,cex=1.7,'G.',line=-1.1)
# 
# #ChDens=function(obs,N) dt(sqrt(N/2)*obs,2*(N-1),sqrt(N/2)*.2)
# ChDens=function(obs,N) dt(sqrt(N)*obs,N-1,sqrt(N)*.2)
# Chpred=ChDens(z,N)
# plot(z,Epred,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
# abline(v=obs,col='grey')
# lines(z,Chpred,lty=2,lwd=2,col='darkgreen')
# text(-.35,.3,expression(paste("Model ",M[0])),adj=1)
# text(.8,.1,expression(paste("Model ",M[d])),adj=0,col='darkgreen')
# points(obs[1],EDens(obs[1],N),cex=1.3,pch=19)
# points(obs[1],ChDens(obs[1],N),cex=1.3,pch=19,col='darkgreen')
# points(obs[2],EDens(obs[2],N),cex=1.3,pch=21,bg='white',lwd=2)
# points(obs[2],ChDens(obs[2],N),cex=1.3,pch=21,col='darkgreen',bg='white',lwd=2)
# mtext(side=3,adj=.05,cex=1.7,'H.',line=-1.1)
# 
# dev.off()
# 
# 
# bf=matrix(ncol=4,c(
#   CDens(obs[1],N)/EDens(obs[1],N),
#   CDens(obs[2],N)/EDens(obs[2],N),
#   UDens(obs[1],N)/EDens(obs[1],N),
#   UDens(obs[2],N)/EDens(obs[2],N),
#   WMDens(obs[1],N)/EDens(obs[1],N),
#   WMDens(obs[2],N)/EDens(obs[2],N),
#   ChDens(obs[1],N)/EDens(obs[1],N),
#   ChDens(obs[2],N)/EDens(obs[2],N)))
# 
# print(c(CDens(obs[1],N),EDens(obs[1],N)))
# print(c(CDens(obs[2],N),EDens(obs[2],N)))
