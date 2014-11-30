source("~/R/bayesFactorCalc2.R") # load Christie, Kaye, & Baguley's version of Dienes' BF calculator

#Bf(sd, obtained, uniform = FALSE, lower=0, upper=1, meanoftheory=0, sdtheory=1, tails=2)

# could also use my equivTest function but that requires the exact per-cell N
# List of BF01's from JZS Bayes, scale=.4
esciTest = function(obtained, n, rscale, paired=F) {
  se = sqrt((1-obtained^2)/(n-2))
  t = obtained / se
  if (paired==F) return(equivTest(ceiling(n/2), floor(n/2), t=t, rscale=rscale))
  else return(equivTestPaired(n, t, rscale=rscale))
} 
# Aggressive affect
# Valadez & Ferguson, interaction effect violent (RDR either version) vs nonviolent (FIFA)
esciTest(.22, 100, .4) # MAY BE INVALIDATED BY UNEQUAL CELL SIZES?
# Valdez & Ferguson, interaction effect as originally reported
esciTest(.17, 100, .4)
# Przybylski et al., Study 1
esciTest(.004, 100, .4)
# Przybylski et al., Study 2
esciTest(.08, 100, .4) # or is it -.08?
# Ivory & Kalyaraman, 2007
esciTest(.13, 120, .4)

# Aggressive behavior
# Elson et al. 2013 noise intensity
esciTest(.2, 84, .4)
# Elson et al. 2013 noise duration
esciTest(.11, 84, .4)
# Ferguson et al. 2008
esciTest(.02, 50, .4)
# Ferguson & Rueda, Violent vs Nonviolent game
esciTest(.01, 77, .4)
# Adachi & Willoughby 2011 exp 1
esciTest(0, 42, .4)
# Adachi & Willoughby 2011 exp 2
esciTest(.03, 60, .4)

# Aggressive cognition
# Ivory & Kalyaraman, 2007
esciTest(-.08, 120, .4)


# List of BF02's from Dienes-Christie calculator
BF02 = function(mean, sd, lower, meanoftheory, sdtheory) {
  return(1/Bf(obtained=mean, sd=sd, lower=lower, meanoftheory=meanoftheory, sdtheory=sdtheory)$BayesFactor)
}
  
# All effect sizes should be in r-to-Z units  
# Aggressive affect
# Valdez & Ferguson, interaction effect violent (RDR either version) vs nonviolent (FIFA)
BF02(sd=.101535, .223656, lower=-1, meanoftheory=.298566, sdtheory=.01996)

# Valdez & Ferguson, interaction effect as originally reported
BF02(sd=.101535, .171667, lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Przybylski et al., Study 1
BF02(sd=.102062, 0, lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Przybylski et al., Study 2
BF02(sd=.101015, .080171, lower=-1, meanoftheory=.298566, sdtheory=.01996)
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


## BF01 / BF10 for study 3 of Elson's CRTT complaint
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
