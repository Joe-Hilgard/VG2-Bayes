source("~/R/bayesFactorCalc2.R") # load Christie, Kaye, & Baguley's version of Dienes' BF calculator

#Bf(sd, obtained, uniform = FALSE, lower=0, upper=1, meanoftheory=0, sdtheory=1, tails=2)

# Aggressive affect
# Valdez & Ferguson, interaction effect violent (RDR either version) vs nonviolent (FIFA)
Bf(.101535, .223656, lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Valdez & Ferguson, interaction effect as originally reported
Bf(.101535, .171667, lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Przybylski et al., Study 1
Bf(.102062, 0, lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Przybylski et al., Study 2
Bf(.101015, .080171, lower=-1, meanoftheory=.298566, sdtheory=.01996)
# Ivory & Kalyraman, 2007
Bf(.09245, .13074, lower=-1, meanoftheory=.298566, sdtheory=.01996)

# Aggressive behavior
# Elson et al. 2013 noise intensity
Bf(.111111, .202733, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Elson et al. 2013 noise duration
Bf(.111111, .110447, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Ferguson et al. 2008
Bf(.145865, .020003, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Ferguson & Rueda, Violent vs Nonviolent game
Bf(.116248, .01, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Adachi & Willoughby 2011 exp 1
Bf(.160128, .0, lower=-1, meanoftheory=.213171, sdtheory=.026252)
# Adachi & Willoughby 2011 exp 2
Bf(.132453, .030009, lower=-1, meanoftheory=.213171, sdtheory=.026252)

# Aggressive cognition
# Ivory & Kalyaraman, 2007
Bf(.09245, -.08017, lower=-1, meanoftheory=.223656, sdtheory=.018621)


