#' The posterior probability criterion function
#'
#' Thall and Simon's criterion function for determining the
#' trial decision cutoffs based on the posterior probability.
#'
#' @usage
#' postprob(y, n, alpha_e, beta_e, alpha_s, beta_s, delta)
#' @param x the number of responses among \eqn{n} patients treated by the experimental drug at a certain stage of the trial.
#' @param n the number of patients treated by the experimental drug at a certain stage of the trial.
#' @param a the hyperparameter (shape1) of the Beta prior for the experimental drug.
#' @param b the hyperparameter (shape2) of the Beta prior for the experimental drug.
#' @param alpha_s the hyperparameter (shape1) of the Beta prior for the standard drug.
#' @param beta_s the hyperparameter (shape2) of the Beta prior for the standard drug.
#' @param delta the minimally acceptable increment of the response rate for the experimental drug compared with the standard drug.
#' @return
#' \item{prob}{the posterior probability: \eqn{Pr(p_E > p_S + \delta | y)}}
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
#' @importFrom stats integrate pbeta dbeta
#' @export
PostProb <- function(x, n, a, b, p0) {
  
  pbeta(p0, a+x, b+n-x, lower.tail = F)

}




#' The stopping boundaries based on Thall and Simon's criterion
#'
#' The stopping boundaries based on Thall and Simon's criterion.
#'
#' @usage
#' stopbound_post(theta, type, nmax, alpha_e, beta_e, alpha_s, beta_s, delta)
#' @param theta the cutoff probability: typically, \eqn{\theta = [0.95, 0.99]} for superiority, \eqn{\theta = [0.01, 0.05]} for futility.
#' @param type type of boundaries: "superiority" or "futility".
#' @param nmax the maximum number of patients treated by the experimental drug.
#' @param alpha_e the hyperparameter (shape1) of the Beta prior for the experimental drug.
#' @param beta_e the hyperparameter (shape2) of the Beta prior for the experimental drug.
#' @param alpha_s the hyperparameter (shape1) of the Beta prior for the standard drug.
#' @param beta_s the hyperparameter (shape2) of the Beta prior for the standard drug.
#' @param delta the minimally acceptable increment of the response rate for the experimental drug compared with the standard drug.
#' Note: if type = "superiority", then delta is set to 0.
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
#' stopbound_post(0.05, "futility", 40, 0.6, 1.4, 15, 35, 0)
#' stopbound_post(0.05, "futility", 30, 0.4, 1.6, 10, 40, 0)
#' stopbound_post(0.95, "superiority", 40, 0.6, 1.4, 15, 35, 0)
#' @export
PostP.design <- function(type = c("efficacy", "futility"), nmax, a, b, p0, delta, theta) {

  type <- match.arg(type)


  bound <- rep(NA, nmax)
  for (n in 1:nmax) {
    if (type == "efficacy") {
      for (x in 0:n) {
        if (PostProb(x, n, a, b, p0) > theta) {
          bound[n] <- x
          break
        }
      }
    } else {
      for (x in n:0) {
        if (PostProb(x, n, a, b, p0+delta) < theta) {
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

