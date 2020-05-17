// vectorized functions for multinomial distribution

#include <Rcpp.h>
using namespace Rcpp;

double log_sum_exp_2(double x, double y);
NumericVector real_to_simplex(const NumericVector &y);
const double TOL = std::sqrt(DOUBLE_EPS);

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
double dmulti_double(const IntegerVector &x,
                     const NumericVector &prob,
                     bool log_p = true) {

  if (x.length() != prob.length()) {
    Rcpp::stop("x and prob must have the same length");
  }
  if (std::abs(Rcpp::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("prob must sum to 1");
  }

  double size = Rcpp::sum(x);

  int n = x.length();

  double lval = R::lgammafn(size + 1.0) - Rcpp::sum(Rcpp::lgamma(x + 1.0));

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
                const int &K,
                const NumericVector &prob,
                bool log_p = true) {

  int minz = std::max(0, gA + gB - K);
  int maxz = std::min(gA, gB);

  if (std::abs(Rcpp::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("prob should sum to 1");
  }

  double lp = R_NegInf;

  IntegerVector x(4);
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
double proballgeno(const IntegerVector &gA,
                   const IntegerVector &gB,
                   const int &K,
                   const NumericVector &prob,
                   bool log_p = true) {

  if (gA.length() != gB.length()) {
    Rcpp::stop("gA and gB need to be the same length");
  }
  int n = gA.length();

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
//' @inheritParams proballgeno
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double llike_geno(const NumericVector &par,
                  const IntegerVector &gA,
                  const IntegerVector &gB,
                  const int &K) {
  if (par.length() != 3) {
    Rcpp::stop("par needs to be length 3");
  }

  NumericVector prob = real_to_simplex(par);

  double llike = proballgeno(gA, gB, K, prob, true);

  return llike;
}


//' Negative Likelihood function meant to be used in lbfgs package by pointer.
//'
//' @param xs This is \code{par} from \code{\link{llike_geno}()}.
//' @param env This is an environtment containing \code{gA},
//'     \code{gB}, and \code{K} from \code{\link{llike_geno}()}.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericVector nllike_geno_p(SEXP xs, SEXP env) {
  NumericVector par(xs);
  Environment e = as<Environment>(env);
  IntegerVector gA = as<IntegerVector>(e["gA"]);
  IntegerVector gB = as<IntegerVector>(e["gB"]);
  int K = as<int>(e["K"]);

  if (par.length() != 3) {
    stop("par needs to be length 3");
  }

  NumericVector prob = real_to_simplex(par);

  NumericVector llike = {-1.0 * proballgeno(gA, gB, K, prob, true)};

  return llike;
}


