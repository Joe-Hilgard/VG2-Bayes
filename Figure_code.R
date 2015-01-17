#library('msm')
#one and two-sample implemented

setwd("~/GitHub/VG2-Bayes/")
N.s = 40  #two-sample, small N
N.m = 100 # med N
N.l = 200 # large N
N.xl = 400 # extra-large N
# N=20 #one-sample

z=seq(-2, 2,.005) # over what interval do we plot?
I=length(z)

# Generate estimates
#Cauchy
s=.43 # set scale
f=dcauchy(z,0,s) # get density function over interval Z for cauchy w/ scale s

# Epred gives prediction of the null model:
EDens=function(obs,N) dt(sqrt(N)/2*obs,(N-2))
# EDens=function(obs,N) dt(sqrt(N)*obs,N-1) # pdf of expected data (will treat N as constant) 
# effect size above is t-value * sqrt(N) so that it is scaled to standard dev. not standard error.
Epred.s=EDens(z,N.s) # run this over z for the expected data given H0
#Epred.m=EDens(z,N.m)
Epred.l=EDens(z,N.l)
#Epred.xl=EDens(z,N.xl)

# Ppred gives prediction of the point-alternative effect model:
PDens = function(obs, N) dt(sqrt(N)/2*obs, N-2, ncp = sqrt(N)/2*.43)
Ppred.s = PDens(z, N.s)
#Ppred.m = PDens(z,N.m)
Ppred.l = PDens(z,N.l)
#Ppred.xl = PDens(z,N.xl)

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
Cpred.s = Cpred.m = Cpred.l = Cpred.xl = 1:I
for (i in 1:I) Cpred.s[i]=CDens(z[i],N.s) 
#for (i in 1:I) Cpred.m[i]=CDens(z[i],N.m)
for (i in 1:I) Cpred.l[i]=CDens(z[i],N.l) 
#for (i in 1:I) Cpred.xl[i]=CDens(z[i],N.xl)

# Npred gives prediction of Anderson's point hypothesis
Nintgrand = function(delta, obs, N) { # function of delta, treating obs and N as constants
  dt(sqrt(N/2)*obs, df = 2*(N-1), ncp = sqrt(N/2)*delta) * dnorm(delta, mean = 0.43, sd = .1)
}
NDens = function(obs, N) integrate(Nintgrand, lower=-20, upper=20, obs=obs, N=N, 
                                   subdivisions=1000, rel.tol = .Machine$double.eps^0.5)$value
Npred.s = Npred.m = Npred.l = Npred.xl = 1:I
for (i in 1:I) Npred.s[i]=NDens(z[i],N.s) 
#for (i in 1:I) Npred.m[i]=NDens(z[i],N.m)
for (i in 1:I) Npred.l[i]=NDens(z[i],N.l) 
#for (i in 1:I) Npred.xl[i]=NDens(z[i],N.xl)

# generate data for figure 2A
BF = function(obs, N) EDens(obs, N)/CDens(obs, N)
BF2 = function(obs, N) EDens(obs, N)/NDens(obs, N)
effect = seq(-1, 1, .005)
BFList.s = BFList.m = BFList.l = BFList.xl = 1:length(effect)
for (i in 1:length(effect)) BFList.s[i] = EDens(effect[i],N.s)/CDens(effect[i],N.s)
for (i in 1:length(effect)) BFList.m[i] = EDens(effect[i],N.m)/CDens(effect[i],N.m)
for (i in 1:length(effect)) BFList.l[i] = EDens(effect[i],N.l)/CDens(effect[i],N.l)
for (i in 1:length(effect)) BFList.xl[i] = EDens(effect[i],N.xl)/CDens(effect[i],N.xl)
# generate data for figure 2B
# Shrinking scale on x-axis for lots of reasons (visibility, but also integration was borked)
effectNarrow = seq(-.2, 1, .005)
BFList2.s = BFList2.m = BFList2.l = BFList2.xl = 1:length(effectNarrow)
for (i in 1:length(effectNarrow)) BFList2.s[i] = EDens(effectNarrow[i],N.s)/NDens(effectNarrow[i],N.s)
for (i in 1:length(effectNarrow)) BFList2.m[i] = EDens(effectNarrow[i],N.m)/NDens(effectNarrow[i],N.m)
for (i in 1:length(effectNarrow)) BFList2.l[i] = EDens(effectNarrow[i],N.l)/NDens(effectNarrow[i],N.l)
for (i in 1:length(effectNarrow)) BFList2.xl[i] = EDens(effectNarrow[i],N.xl)/NDens(effectNarrow[i],N.xl)############

## Figure 1
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
arrows(.43,0,.43,1) # Draw Dirac for point-alternative hypothesis
# arrows(.21, 0, .21, 1, col="red") # Draw Dirac for Anderson's hypothesis (not truly Dirac but might as well be)
# text labels
# Add text for red dirac arrow?
text(-.2,.8,expression(paste("Model ",M[0])),adj=1)
text(.3,.25,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
mtext(side=3,adj=.05,cex=1.5,'A.',line=-1.2)

# Probabilities given data:
plot(z,Epred.s,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
lines(z,Ppred.s,lty=3,lwd=2,col='darkred') # data given point-alternative
obs = .4
abline(v=obs,col='grey')
points(obs,EDens(obs,N.s),cex=1.3,pch=19)
points(obs,PDens(obs,N.s),cex=1.3,pch=19)
PDens(obs,N.s)/EDens(obs,N.s)
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
plot(z,Epred.s,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
lines(z,Cpred.s,lty=3,lwd=2,col='darkred') # data given JZS-prior alternative
obs = .4
abline(v=obs,col='grey')
points(obs,EDens(obs,N.s),cex=1.3,pch=19)
points(obs,CDens(obs,N.s),cex=1.3,pch=19)
text(-.4,.2,expression(paste("Model ",M[0])),adj=1)
text(.4,.1,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
mtext(side=3,adj=.05,cex=1.5,'D.',line=-1.2)
EDens(obs,N.s)/CDens(obs,N.s)

# Plotting, row 3
# H0 vs Craig meta with big N
plot(z,Epred.l,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
lines(z,Ppred.l,lty=3,lwd=2,col='darkred') # data given point-alternative
obs = .4
abline(v=obs,col='grey')
points(obs,EDens(obs,N.l),cex=1.3,pch=19)
points(obs,PDens(obs,N.l),cex=1.3,pch=19)
PDens(obs,N.l)/EDens(obs,N.l)
text(-.35,.3,expression(paste("Model ",M[0])),adj=1)
text(.9,.15,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
mtext(side=3,adj=.05,cex=1.5,'E.',line=-1.2)
# H0 vs vague alternative with big N
plot(z,Epred.l,typ='l',xlab=expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),ylab="Density")
lines(z,Cpred.l,lty=3,lwd=2,col='darkred') # data given JZS-prior alternative
obs = .4
abline(v=obs,col='grey')
points(obs,EDens(obs,N.l),cex=1.3,pch=19)
points(obs,CDens(obs,N.l),cex=1.3,pch=19)
text(-.4,.2,expression(paste("Model ",M[0])),adj=1)
text(.4,.1,expression(paste("Model ",M[a])),col="darkgreen",adj=0)
mtext(side=3,adj=.05,cex=1.5,'F.',line=-1.2)
EDens(obs,N.l)/CDens(obs,N.l)

dev.off()

## Figure 2:

# Prepare plot space
pdf('BFFigure2.pdf') # attempting to scale approximately to A4 dimensions #width=12,height=24) 
par(mfrow=c(1,2),cex=1,mar=c(4,4,.5,1),mgp=c(2.2,1,0)) # 

# BF10 vs observed effect size:
plot(x=effect, xlab = expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),
     y=BFList.s, log="y", ylim=c(10, .00001), ylab="Bayes factor (alternative/null)",
     typ='l',
     yaxt='n'
     )
axis(2, at=c(.0001, .001, .01, .1, 1, 10), 
     labels=c(".0001", ".001", ".01", ".1", "1", "10"),
     )
abline(h=1)
lines(x=effect, y=BFList.m, lty=2, col=2)
lines(x=effect, y=BFList.l, lty=3, col=3)
lines(x=effect, y=BFList.xl, lty=4, col=4)
#mtext(side=3,adj=.05,cex=1.5,'A.',line=-1.2)

# BF20 vs observed effect size:
plot(x=effectNarrow, xlab = expression(paste("Observed Effect Size, ",(bar(y)-bar(x))/s)),
     y=BFList2.xl, log="y", #ylim=c(1000, .00001)
     , ylab="Bayes factor (alternative/null)",
     typ='l'
     ,yaxt='n'
)
abline(h=1)
axis(2, at=c(.0001, .001, .01, .1, 1, 10, 100, 1000), 
     labels=c(".0001", ".001", ".01", ".1", "1", "10", "100", "1000"),
)
lines(x=effectNarrow, y=BFList2.s, lty=2, col="darkred")
lines(x=effectNarrow, y=BFList2.m, lty=3, col="darkgreen")
lines(x=effectNarrow, y=BFList2.l, lty=4, col="darkgreen")
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
