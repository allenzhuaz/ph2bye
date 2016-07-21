#' The posterior probability criterion function for Phase II single-arm design
#'
#' Thall and Simon's criterion function for determining the trial decision boundaries based on the posterior probability.
#'
#' @usage
#' PostP(x, n, a, b, p0)
#' @param x the number of responses among \eqn{n} patients treated by the experimental drug at a certain stage of the trial.
#' @param n the number of patients treated by the experimental drug at a certain stage of the trial.
#' @param a the hyperparameter (shape1) of the Beta prior for the experimental drug.
#' @param b the hyperparameter (shape2) of the Beta prior for the experimental drug.
#' @param p0 the prespecified reseponse rate.
#' @return
#' \item{prob}{the posterior probability: \eqn{Pr(p > p_0  | X=x)}}
#' @references
#' Berry, S. M., Carlin, B. P., Lee, J. J., & Muller, P. (2010).
#' \emph{Bayesian adaptive methods for clinical trials.}
#' CRC press.
#'
#' Thall, P. F., Simon, R. (1994).
#' Practical Bayesian guidelines for phase IIB clinical trials.
#' \emph{Biometrics} \strong{50}: 337-349.
#'
#' Yin, G. (2013).
#' \emph{Clinical Trial Design: Bayesian and Frequentist Adaptive Methods.}
#' New York: Wiley.
#' @examples
#' PostP(8,15,1,1,0.8)
#' @importFrom stats pbeta
#' @export
PostP <- function(x, n, a, b, p0) {
  pbeta(p0, a+x, b+n-x, lower.tail = F)
}




#' The stopping boundaries based on the posterior probability criterion
#'
#' The design function to calucluate the sequential monitor sample size and boundary based on Thall and Simon's criterion.
#'
#' @usage
#' PostP.design(type, nmax, a, b, p0, delta, theta)
#' @param type type of boundaries: "superiority" or "futility".
#' @param nmax the maximum number of patients treated by the experimental drug.
#' @param a the hyperparameter (shape1) of the Beta prior for the experimental drug.
#' @param b the hyperparameter (shape2) of the Beta prior for the experimental drug.
#' @param p0 the pre-specified reseponse rate.
#' @param delta the minimally acceptable increment of the response rate for the experimental drug compared with the standard drug.
#' @param theta the cutoff probability: typically, \eqn{\theta = [0.95, 0.99]} for superiority, \eqn{\theta = [0.01, 0.05]} for futility.
#' @return
#' \item{boundset}{the boundaries set; \eqn{U_n} or \eqn{L_n}}
#' @references
#' Thall, P. F., Simon, R. (1994).
#' Practical Bayesian guidelines for phase IIB clinical trials.
#' \emph{Biometrics} \strong{50}: 337-349.
#'
#' Yin, G. (2012).
#' \emph{Clinical Trial Design: Bayesian and Frequentist Adaptive Methods.}
#' New York: Wiley.
#' @examples
#' ## Using vague prior Unif(0,1)
#' PostP.design(type = "futility", nmax=100, a=1, b=1, p0=0.15, delta=0.15, theta=0.05)
#' PostP.design(type = "efficacy", nmax=100, a=1, b=1, p0=0.15, delta=0.15, theta=0.9)
#' ## Or using Jeffery prior with Beta(0.5,0.5)
#' PostP.design(type = "futility", nmax=100, a=0.5, b=0.5, p0=0.15, delta=0.15, theta=0.05)
#' PostP.design(type = "efficacy", nmax=100, a=0.5, b=0.5, p0=0.15, delta=0.15, theta=0.9)
#' @export
PostP.design <- function(type = c("efficacy", "futility"), nmax, a, b, p0, delta=0, theta) {
  type <- match.arg(type)
  bound <- rep(NA, nmax)
  for (n in 1:nmax) {
    if (type == "efficacy") {
      for (x in 0:n) {
        if (PostP(x, n, a, b, p0) > theta) {
          bound[n] <- x
          break
        }
      }
    } else {
      for (x in n:0) {
        if (PostP(x, n, a, b, p0+delta) < theta) {
          bound[n] <- x
          break
        }
      }
    }
  }

  boundset <- data.frame(n = 1:nmax, bound = bound)

  return(boundset[!duplicated(boundset[, 2]), ])

}

PostBound <- function(type = c("efficacy", "futility"),a,b,p0, theta){
  type <- match.arg(type)

}


oc.bdry <- function(pu, pa, r1, n1, r, n){
  pet <- err0 <- pbinom(r1,n1,pu)
  ess <- n1 + (n-n1)*(1-pet)
  err1 <- pbinom(r1,n1,pa)
  for(i in (r1+1):r) {
    err0 <- err0 + dbinom(i,n1,pu)*pbinom(r-i,n-n1,pu)
    err1 <- err1 + dbinom(i,n1,pa)*pbinom(r-i,n-n1,pa)
  }
  out <- c(1-err0, 1-err1, pet, ess)
  names(out) <- c("P(reject H0 | p0)","P(reject H0 | p1)","PET(p0)","EN(p0)")
  out
}

