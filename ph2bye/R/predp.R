#' The predictive probability criterion function
#'
#' Lee and Liu's criterion function for determining the
#' trial decision cutoffs based on the predictive probability.
#'
#' @usage
#' predprob(y, n, nmax, alpha_e, beta_e, p_s, theta_t)
#' @param y the number of responses among \eqn{n} patients treated by the experimental drug at a certain stage of the trial.
#' @param n the number of patients treated by the experimental drug at a certain stage of the trial.
#' @param nmax the maximum number of patients treated by the experimental drug.
#' @param alpha_e the hyperparameter (shape1) of the Beta prior for the experimental drug.
#' @param beta_e the hyperparameter (shape2) of the Beta prior for the experimental drug.
#' @param p_s the the response rate for the standard drug.
#' @param theta_t the prespecified target probability; tipically, \eqn{\theta_T = [0.85, 0.95]}.
#' @return
#' \item{prob}{the predictive probability: \eqn{PP = \sum_{x=0}^{n_{max}-n} P(x | y) I(\Pr(p_E > p_S | y, x) \geq \theta_T) }}
#' @references
#' Lee, J. J., Liu, D. D. (2008).
#' A predictive probability design for phase II cancer clinical trials.
#' \emph{Clinical Trials} \strong{5}: 93-106.
#'
#' Yin, G. (2012).
#' \emph{Clinical Trial Design: Bayesian and Frequentist Adaptive Methods.}
#' New York: Wiley.
#' @examples
#' # p. 97, PP = 0.5656
#' predprob(16, 23, 40, 0.6, 0.4, 0.6, 0.9)
#' @export

PredProb <- function(y, n, nmax, alpha_e, beta_e, p_s, theta_t) {
  
  return(predprobCpp(y, n, nmax, alpha_e, beta_e, p_s, theta_t))
  
}


#' The stopping boundaries based on Lee and Liu's criterion
#'
#' The stopping boundaries based on Lee and Liu's criterion.
#'
#' @usage
#' stopbound_pred(theta, type, nmax, alpha_e, beta_e, p_s, theta_t)
#' @param theta the cutoff probability: typically, \eqn{\theta = [0.95, 0.99]} for superiority, \eqn{\theta = [0.01, 0.05]} for futility.
#' @param type type of boundaries: "superiority" or "futility".
#' @param nmax the maximum number of patients treated by the experimental drug.
#' @param alpha_e the hyperparameter (shape1) of the Beta prior for the experimental drug.
#' @param beta_e the hyperparameter (shape2) of the Beta prior for the experimental drug.
#' @param p_s the the response rate for the standard drug.
#' @param theta_t the prespecified target probability; tipically, \eqn{\theta_T = [0.85, 0.95]}.
#' @return
#' \item{boundset}{the boundaries set: \eqn{U_n} or \eqn{L_n}}
#' @references
#' Lee, J. J., Liu, D. D. (2008).
#' A predictive probability design for phase II cancer clinical trials.
#' \emph{Clinical Trials} \strong{5}: 93-106.
#'
#' Yin, G. (2012).
#' \emph{Clinical Trial Design: Bayesian and Frequentist Adaptive Methods.}
#' New York: Wiley.
#' @examples
#' stopbound_pred(0.05, "futility", 40, 0.6, 1.4, 0.3, 0.85)
#' stopbound_pred(0.05, "futility", 30, 0.4, 1.6, 0.2, 0.85)
#' stopbound_pred(0.95, "superiority", 40, 0.6, 1.4, 0.3, 0.85)
#' @export
PredP.design <- function(type = c("efficacy", "futility"), nmax, a, b, p0, theta_t=0.9, delta, theta) {
  
  type <- match.arg(type)
  if (type == "efficacy") {
  return(ph2bayes::stopbound_pred(theta = theta,type = "superiority",nmax = nmax,alpha_e = a,beta_e = b,p_s = p0,theta_t = theta_t)) 
  } else {
  return(ph2bayes::stopbound_pred(theta = theta,type = "futility",nmax = nmax,alpha_e = a,beta_e = b,p_s = p0+delta,theta_t = theta_t)) 
  }  
}
  
