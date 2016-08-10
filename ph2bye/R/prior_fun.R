prior <- function(type=c("MeanVar", "Optimist", "ORRN", "ORRW"), mu, v, N, W, init=c(0.5,0.5)){
  type <- match.arg(type)
  switch (type,
    MeanVar = {c(a = mu^2*(1-mu)/v - mu, b = muu*(1-mu)^2/v - (1-mu))},
    Optimist = {c(a = mu+1,b=2-mu)},
    ORRN = {c(a=mu+1+mu*N, b=(1-mu)+1+N*(1-mu))},
    ORRW = {
      sol <- function(x){
      y <- numeric(2)
      y[1] <- x[1]/(x[1]+x[2])-mu
      y[2] <- qbeta(p = 0.975, shape1 = x[1], shape2 = x[2])-qbeta(p = 0.025, shape1 = x[1], shape2 = x[2])-W
      y
    }
      res <- nleqslv::nleqslv(init, sol, control=list(trace=1,btol=.01,delta="newton"))$'x'
      c(a=res[1],b=res[2])
  })
}
