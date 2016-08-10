## Double threshold design showing futility and efficacy boundary together

DT.desgin <- function(type=c("PostP","PredP"), a, b, nmax, p.vec, theta.vec, theta_t=0.9, optimize=FALSE){
  type <- match.arg(type)
  if (theta.vec[1]>theta.vec[2]) stop("The efficacy threshold should be greater than futility threshold.")
  res <- switch(type,
  PostP = {
    merge(PostP.design(type = "futility", a = a, b = b, nmax = nmax, p0 = p.vec[1], theta = theta.vec[1], optimize=optimize),
          PostP.design(type = "efficacy", a = a, b = b, nmax = nmax, p0 = p.vec[2], theta = theta.vec[2], optimize=optimize),
          by="n")},
  PredP = {
    merge(PredP.design(type = "futility", a = a, b = b, nmax = nmax, p0 = p.vec[1], theta = theta.vec[1], theta_t = theta_t, optimize=optimize),
          PredP.design(type = "efficacy", a = a, b = b, nmax = nmax, p0 = p.vec[2], theta = theta.vec[2], theta_t = theta_t, optimize=optimize),
          by="n")})
  names(res) <- c("n", "futility", "efficacy")
  return(res)
  }

a <- prior(type = "ORRN", mu = TT.mean, N = TT.n)[1]
b <- prior(type = "ORRN", mu = TT.mean, N = TT.n)[2]


a <- prior(type = "ORRN", mu = KT.mean, N = KT.n)[1]
b <- prior(type = "ORRN", mu = KT.mean, N = KT.n)[2]

a <- prior(type = "ORRW", mu = TT.mean, W = TT.W)[1]
b <- prior(type = "ORRW", mu = TT.mean, W = TT.W)[2]
