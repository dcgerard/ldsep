// derivative functions

#include <Rcpp.h>
using namespace Rcpp;

double expit(double x);
double log_sum_exp_2(double x, double y);
double dmulti_double(const IntegerVector &x, const NumericVector &prob, bool log_p);
NumericVector plog_sum_exp(const NumericVector &x, const NumericVector &y);
NumericVector real_to_simplex(const NumericVector &y);

//' Derivative of multinomial pdf
//'
//' Gradient of multinomial pdf (NOT log-pdf) with respect to prob,
//' optionally returned on the log-scale.
//'
//' @inheritParams dmulti_double
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericVector dmulti_dprob(const IntegerVector &x,
                           const NumericVector &prob,
                           bool log_p = true) {
  if (x.length() != prob.length()) {
    Rcpp::stop("x and prob must have the same length");
  }

  double size = Rcpp::sum(x);
  int n = x.length();

  NumericVector deriv(4);
  deriv = deriv + R::lgammafn(size + 1.0) - Rcpp::sum(Rcpp::lgamma(x + 1.0));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) {
        deriv[j] += (x[i] - 1) * std::log(prob[i]) + std::log(x[i]);
      } else {
        deriv[j] += x[i] * std::log(prob[i]);
      }
    }
  }

  if (!log_p) {
    deriv = Rcpp::exp(deriv);
  }

  return deriv;
}

//' Derivative of \code{\link{probgeno}()} with respect to \code{prob}.
//'
//' @inheritParams probgeno
//'
//' @author David Gerard
//'
//' @return A vector of length 4 containing the gradient of
//'     \code{\link{probgeno}()} with respect to \code{prob}.
//'
//' @noRd
// [[Rcpp::export]]
NumericVector dprobgeno_dprob(const int &gA,
                              const int &gB,
                              const int &K,
                              const NumericVector &prob) {
  int minz = std::max(0, gA + gB - K);
  int maxz = std::min(gA, gB);

  NumericVector deriv = NumericVector::create(R_NegInf, R_NegInf, R_NegInf, R_NegInf);
  double logdenom = R_NegInf;
  IntegerVector x(4);
  for (int z = minz; z <= maxz; z++) {
    x[0] = K + z - gA - gB; // number ab
    x[1] = gA - z; // number Ab
    x[2] = gB - z; // number aB
    x[3] = z; // number AB
    logdenom = log_sum_exp_2(logdenom, dmulti_double(x, prob, true));
    deriv = plog_sum_exp(deriv, dmulti_dprob(x, prob, true));
  }
  deriv = Rcpp::exp(deriv - logdenom);
  return deriv;
}

//' Derivative of \code{\link{proballgeno}()} with respect to \code{prob}.
//'
//' @inheritParams proballgeno
//'
//' @author David Gerard
//'
//' @return A vector of length 4 containing the gradient of
//'     \code{\link{proballgeno}()} with respect to \code{prob}.
//'
//' @noRd
// [[Rcpp::export]]
NumericVector dproballgeno_dprob(const IntegerVector &gA,
                                 const IntegerVector &gB,
                                 const int &K,
                                 const NumericVector &prob) {
  if (gA.length() != gB.length()) {
    Rcpp::stop("gA and gB need to be the same length");
  }
  int n = gA.length();

  NumericVector deriv(4);
  for (int i = 0; i < n; i++) {
    deriv = deriv + dprobgeno_dprob(gA[i], gB[i], K, prob);
  }

  return deriv;
}


//' Derivative of \code{\link{real_to_simplex}()} with respect to \code{y}.
//'
//' @param y A numeric matrix.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericMatrix dreal_to_simplex_dy(const NumericVector &y) {
  int K = y.length() + 1;

  NumericMatrix jacob(K, K - 1);

  double recsum;
  double drecsum;
  double zk;
  for (int j = 0; j < K - 1; j++) {
    recsum = 0.0;
    drecsum = 0.0;
    for (int i = 0; i < K; i++) {
      zk = expit(y[i] + std::log(1.0 / ((double)K - ((double)i + 1))));
      if (i < j) {
        jacob(i, j) = 0.0;
      } else if (i == j) {
        jacob(i, j) = (1 - recsum) * zk * (1.0 - zk);
      } else {
        jacob(i, j) = - drecsum * zk;
      }
      recsum += (1.0 - recsum) * zk;
      drecsum += jacob(i, j);
    }
  }

  return jacob;
}


//' Derivative of \code{\link{llike_geno}()} with respect to par.
//'
//'
//' @inheritParams llike_geno
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericVector dllike_geno_dpar(const NumericVector &par,
                               const IntegerVector &gA,
                               const IntegerVector &gB,
                               const int &K) {
  if (par.length() != 3) {
    Rcpp::stop("par needs to be length 3");
  }

  NumericMatrix dp_dy = dreal_to_simplex_dy(par);
  NumericVector prob = real_to_simplex(par);
  NumericVector df_dp = dproballgeno_dprob(gA, gB, K, prob);

  NumericVector deriv(3);
  for (int i = 0; i < 4; i++) {
    deriv[0] += df_dp[i] * dp_dy(i, 0);
    deriv[1] += df_dp[i] * dp_dy(i, 1);
    deriv[2] += df_dp[i] * dp_dy(i, 2);
  }

  return deriv;
}

//' Negative Derivative function meant to be used in lbfgs package by pointer.
//'
//' @param xs This is \code{par} from \code{\link{dllike_geno_dpar}()}.
//' @param env This is an environtment containing \code{gA},
//'     \code{gB}, and \code{K} from \code{\link{dllike_geno_dpar}()}.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericVector dnllike_geno_dpar_p(SEXP xs, SEXP env) {
  NumericVector par(xs);
  Environment e = as<Environment>(env);
  IntegerVector gA = as<IntegerVector>(e["gA"]);
  IntegerVector gB = as<IntegerVector>(e["gB"]);
  int K = as<int>(e["K"]);

  NumericVector deriv = -1.0 * dllike_geno_dpar(par, gA, gB, K);

  return deriv;
}
