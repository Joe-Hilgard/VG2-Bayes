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
effectNarrow = seq(-0, .8, .005)
BFList2.s = BFList2.m = BFList2.l = BFList2.xl = 1:length(effectNarrow)
for (i in 1:length(effectNarrow)) BFList2.s[i] = EDens(effectNarrow[i],N.s)/NDens(effectNarrow[i],N.s)
for (i in 1:length(effectNarrow)) BFList2.m[i] = EDens(effectNarrow[i],N.m)/NDens(effectNarrow[i],N.m)
for (i in 1:length(effectNarrow)) BFList2.l[i] = EDens(effectNarrow[i],N.l)/NDens(effectNarrow[i],N.l)
for (i in 1:length(effectNarrow)) BFList2.xl[i] = EDens(effectNarrow[i],N.xl)/NDens(effectNarrow[i],N.xl)############
