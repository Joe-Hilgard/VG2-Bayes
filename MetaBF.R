#Jeff Rouder
#Nov., 2010
#Aug, 2013, added "r" option.

#Bayes factor across a series of experiments

#see bottom for usage


########################################################
#Internal Functions

t.value.bf=function(t,n1,n2=0,sd=1,prior.cauchy=T)
{
nu=ifelse(n2==0 | is.null(n2),n1-1,n1+n2-2)
n=ifelse(n2==0 | is.null(n2),n1,(n1*n2)/(n1+n2))
r2=sd*sd
marg.like.0=(1+t^2/(nu))^(-(nu+1)/2)
marg.like.1=ifelse(prior.cauchy,
 integrate(t.joint,lower=0,upper=Inf,t=t,n=n,nu=nu,r2=r2)$value,
  (1+n*r2)^(-.5)*(1+t^2/((1+n*r2)*(nu)))^(-(nu+1)/2))
return(marg.like.0/marg.like.1)
}


wrapper=function(tval,delta,N) dt(tval,df=N-1,ncp=delta*sqrt(N),log=T)

mv.integ = function(delta,tval,r,N){
temp=outer(tval,delta,wrapper,N=N)
like=apply(temp,2,sum)
exp(like + dcauchy(delta,scale=r,log=TRUE))
}

metaBF=function(tval,N,half=T,r=1){
num = exp(sum(dt(tval,df=N-1,log=T)))
lower=ifelse(half,0,-Inf)
coeff=ifelse(half,2,1)
denom = coeff*integrate(mv.integ,lower,Inf,r,N=N,tval=tval)[[1]]
return(num/denom)
}



dinvgamma=function (x, shape, scale = 1) 
{
    if (shape <= 0 | scale <= 0) {
        stop("Shape or scale parameter negative in dinvgamma().\n")
    }
    alpha <- shape
    beta <- scale
    log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 
        1) * log(x) - (beta/x)
    return(exp(log.density))
}




t.joint=function(g,t,n,nu,r2)
{
t1=-.5*log(1+n*g*r2)
t2=(-(nu+1)/2)*log(1+t^2/((1+n*g*r2)*(nu)))
return(dinvgamma(g,.5,.5)*exp(t1+t2))
}

#END INTERNAL FUNCTIONS
##################################################


##################################################
#Usage
##################################################


#output is null over alternative
#USAGE For TABLE 2 of Rouder and Morey, 2011
N=c(100,100,100,100)
t=c(2.16,1.25,2.80,1.44)
alt.2.tail=metaBF(t,N,half=F)
tab2.bf=1/alt.2.tail  #near value reported by Rouder and Morey (difference is rounding)

#for one-tailed alternatives (null is point plus all negative effects)
alt.1.tail=metaBF(t,N,half=T)
#or
1/alt.1.tail

#for one-tailed alternatives w smaller r
alt.1.tail=metaBF(t,N,half=T,r=.5)
#or
1/alt.1.tail


