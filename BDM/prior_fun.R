betaprior <- function(m,v){
  a <- m^2*(1-m)/v - m; b <- m*(1-m)^2/v - (1-m)
  return(c(a=a,b=b))
}
#KT.prior <- betaprior(KT.mean, KT.var)
#KT.prior
#OD.prior <- betaprior(OD.mean, OD.var)
#OD.prior
#binom.test(x=round(KT.n*KT.mean),n = KT.n, p=0.2)
#binom.test(x=round(OD.n*OD.mean),n = OD.n, p=0.2)


optprior <- function(x){
  return(c(a=x+1,b=2-x))
}

ORRNprior <- function(x,y){
  return(c(a=x+1+x*y, b=(1-x)+1+y*(1-x)))
}

library(nleqslv)
ORRWprior <- function(mu,W, xstart = c(0.5,0.5)){
  prior <- function(x){
    y <- numeric(2)
    y[1] <- x[1]/(x[1]+x[2])-mu
    y[2] <- qbeta(p = 0.975, shape1 = x[1], shape2 = x[2])-qbeta(p = 0.025, shape1 = x[1], shape2 = x[2])-W
    y
  }
  nleqslv(xstart, prior, control=list(trace=1,btol=.01,delta="newton"))$'x'
}