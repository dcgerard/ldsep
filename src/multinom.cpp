// vectorized functions for multinomial distribution

#include <RcppArmadillo.h>
using namespace Rcpp;

double log_sum_exp_2(double x, double y);
arma::vec real_to_simplex(const arma::vec y);
const double TOL = std::sqrt(DOUBLE_EPS);
double lprior(const arma::vec prob, const arma::vec alpha);

//' Multinomial pdf
//'
//' @param x A vector of counts
//' @param prob A vector of probabilities
//' @param log_p A logical. Should we return the log (\code{TRUE}) or not
//'     (\code{FALSE})?
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double dmulti_double(const arma::vec x,
                     const arma::vec prob,
                     bool log_p = true) {

  if (x.n_elem != prob.n_elem) {
    Rcpp::stop("dmulti_double: x and prob must have the same length");
  }
  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("dmulti_double: prob must sum to 1");
  }

  double size = arma::sum(x);

  int n = x.n_elem;

  double lval = R::lgammafn(size + 1.0) - arma::sum(arma::lgamma(x + 1.0));

  for (int i = 0; i < n; i++) {
    if (x[i] > 0) {
      lval += x[i] * std::log(prob[i]);
    }
  }

  if (!log_p) {
    lval = std::exp(lval);
  }

  return lval;
}

//' Probability of genotypes given haplotype frequencies for one individual
//'
//' @param gA The number of A alleles.
//' @param gB The number of B alleles.
//' @param K The ploidy of the species.
//' @param prob A numeric vector with the probabilities of haplotypes
//'     (in order) of (ab, Ab, aB, AB).
//' @param log_p A logical. log_p A logical. Should we return the log
//'      (\code{TRUE}) or not (\code{FALSE})?
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double probgeno(const int &gA,
                const int &gB,
                const int K,
                const arma::vec prob,
                bool log_p = true) {

  int minz = std::round(std::max(0, gA + gB - K));
  int maxz = std::round(std::min(gA, gB));

  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("probgeno: prob should sum to 1");
  }

  double lp = -arma::datum::inf;

  arma::vec x(4);
  for (int z = minz; z <= maxz; z++) {
    x[0] = K + z - gA - gB; // number ab
    x[1] = gA - z; // number Ab
    x[2] = gB - z; // number aB
    x[3] = z; // number AB

    lp = log_sum_exp_2(lp, dmulti_double(x, prob, true));
  }

  if (!log_p) {
    lp = std::exp(lp);
  }

  return lp;
}

//' Probability of genotypes given haplotype frequencies for all individuals
//'
//' @param gA Vector of number of A alleles.
//' @param gB Vector of number of B alleles.
//' @inheritParams probgeno
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double proballgeno(const arma::vec &gA,
                   const arma::vec &gB,
                   const int K,
                   const arma::vec prob,
                   bool log_p = true) {

  if (gA.n_elem != gB.n_elem) {
    Rcpp::stop("proballgeno: gA and gB need to be the same length");
  }
  int n = gA.n_elem;

  double lp = 0.0;
  for (int i = 0; i < n; i++) {
    lp += probgeno(gA[i], gB[i], K, prob, log_p = true);
  }

  if (!log_p) {
    lp = std::exp(lp);
  }

  return lp;
}

//' Likelihood function when estimating LD directly from genotypes
//'
//' @param par A vector of length 3 containing real numbers that are to
//'     be transformed into the simplex of prob (ab, Ab, aB, AB).
//' @param alpha The prior sample size used in the penalty.
//' @inheritParams proballgeno
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double llike_geno(const arma::vec par,
                  const arma::vec &gA,
                  const arma::vec &gB,
                  const int K,
                  const arma::vec alpha) {
  if (par.n_elem != 3) {
    Rcpp::stop("llike_geno: par needs to be length 3");
  }

  arma::vec prob = real_to_simplex(par);

  double llike = proballgeno(gA, gB, K, prob, true) +
    lprior(prob, alpha);

  return llike;
}



