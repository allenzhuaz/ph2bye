#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double predprobCpp(int y, int n, int nmax, double a, double b, double p0, double theta_t) {

  double prob = 0.0;
  double eps = std::numeric_limits<double>::epsilon();
  double pxy;

  for (int x = 0; x < nmax - n + 1; x++) {
    pxy = (1.0 - R::pbeta(p0, a + y + x, b + nmax - y - x, 1, 0));
    if (pxy > theta_t || std::abs(pxy - theta_t) < eps) {
      prob += exp(
        R::lchoose(nmax - n, x) +
          R::lbeta(a + y + x, b + nmax - y - x) -
          R::lbeta(a + y, b + n - y)
      );
    }
  }

  return prob;

}
