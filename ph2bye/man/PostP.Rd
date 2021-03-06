% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postp.R
\name{PostP}
\alias{PostP}
\title{The posterior probability criterion function for Phase II single-arm design}
\usage{
PostP(x, n, a, b, p0)
}
\arguments{
\item{x}{the number of responses among \eqn{n} patients treated by the experimental drug.}

\item{n}{the number of patients treated by the experimental drug.}

\item{a}{the hyperparameter (shape1) of the Beta prior for the experimental drug.}

\item{b}{the hyperparameter (shape2) of the Beta prior for the experimental drug.}

\item{p0}{the prespecified reseponse rate.}
}
\value{
\item{prob}{the posterior probability: \eqn{Pr(p > p_0  | X=x)}}
}
\description{
Thall and Simon's criterion function for determining the trial decision boundaries based on the posterior probability.
}
\examples{
PostP(8,15,1,1,0.8)
}
\references{
Berry, S. M., Carlin, B. P., Lee, J. J., & Muller, P. (2010).
\emph{Bayesian adaptive methods for clinical trials.}
CRC press.

Thall, P. F., Simon, R. (1994).
Practical Bayesian guidelines for phase IIB clinical trials.
\emph{Biometrics} \strong{50}: 337-349.

Yin, G. (2013).
\emph{Clinical Trial Design: Bayesian and Frequentist Adaptive Methods.}
New York: Wiley.
}

