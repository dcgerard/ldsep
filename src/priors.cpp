#include <RcppArmadillo.h>
using namespace Rcpp;

const double TOL = std::sqrt(DOUBLE_EPS);
arma::vec real_to_simplex(const arma::vec y);
arma::mat dreal_to_simplex_dy(const arma::vec y);

//' Prior probability for haplotype frequencies.
//'
//' @param prob The vector of probabilities for haplotypes (ab, Ab, aB, AB).
//' @param alpha The prior counts for \code{prob}.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double lprior(const arma::vec prob, const arma::vec alpha) {
  if (prob.n_elem != 4) {
    Rcpp::stop("lprior: prob must be of length 4");
  }
  if (alpha.n_elem != 4) {
    Rcpp::stop("lprior: prob must be of length 4");
  }
  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("lprior: prob should sum to 1");
  }

  double pl = std::lgamma(arma::sum(alpha)) -
    (double)arma::sum(arma::lgamma(alpha)) +
    (double)arma::sum((alpha - 1.0) % arma::log(prob));

  return pl;
}

//' Derivative of lprior with respect to prob.
//'
//' @param prob The vector of probabilities for haplotypes (ab, Ab, aB, AB).
//' @param alpha The prior counts for \code{prob}.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dlprior_dprob(const arma::vec prob, const arma::vec alpha) {
  if (prob.n_elem != 4) {
    Rcpp::stop("lprior: prob must be of length 4");
  }
  if (alpha.n_elem != 4) {
    Rcpp::stop("lprior: prob must be of length 4");
  }
  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("lprior: prob should sum to 1");
  }

  arma::vec deriv = (alpha - 1.0) / prob;

  return deriv;
}

//' Prior probability on real-line scale.
//'
//' @param par A vector of length 3 containing real numbers that are to
//'     be transformed into the simplex of prob (ab, Ab, aB, AB).
//' @param alpha The prior counts for \code{prob}.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double lprior_par(const arma::vec par,
                  const arma::vec alpha) {
  if (par.n_elem != 3) {
    Rcpp::stop("lprior_par: par needs to be length 3");
  }

  arma::vec prob = real_to_simplex(par);

  double pl = lprior(prob, alpha);

  return pl;
}


//' Derivative of \code{\link{lprior_par}()} with respect to \code{par}.
//'
//' @inheritParams lprior_par
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dlprior_par_dprob(const arma::vec par,
                            const arma::vec alpha) {
  if (par.n_elem != 3) {
    Rcpp::stop("dlprior_par_dprob: par needs to be length 3");
  }

  arma::mat dp_dy = dreal_to_simplex_dy(par);
  arma::vec prob = real_to_simplex(par);
  arma::vec df_dp = dlprior_dprob(prob, alpha);

  arma::vec deriv = {0.0, 0.0, 0.0};
  for (int i = 0; i < 4; i++) {
    deriv[0] += df_dp[i] * dp_dy(i, 0);
    deriv[1] += df_dp[i] * dp_dy(i, 1);
    deriv[2] += df_dp[i] * dp_dy(i, 2);
  }

  return deriv;
}


