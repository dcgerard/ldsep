// derivative functions

#include <RcppArmadillo.h>
using namespace Rcpp;

double expit(double x);
double log_sum_exp_2(double x, double y);
double dmulti_double(const arma::vec x, const arma::vec prob, bool log_p);
arma::vec plog_sum_exp(const arma::vec &x, const arma::vec &y);
arma::vec real_to_simplex(const arma::vec y);
const double TOL = std::sqrt(DOUBLE_EPS);
arma::vec dlprior_dprob(const arma::vec prob, const arma::vec alpha);

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
arma::vec dmulti_dprob(const arma::vec x,
                       const arma::vec prob,
                       bool log_p = true) {
  if (x.n_elem != prob.n_elem) {
    Rcpp::stop("x and prob must have the same length");
  }

  double size = arma::sum(x);
  int n = x.n_elem;

  arma::vec deriv = {0.0, 0.0, 0.0, 0.0};
  deriv = deriv + R::lgammafn(size + 1.0) - arma::sum(arma::lgamma(x + 1.0));

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
    deriv = arma::exp(deriv);
  }

  return deriv;
}

//' Derivative of \code{\link{probgeno}(log = TRUE)} with respect to
//' \code{prob}.
//'
//' @inheritParams probgeno
//' @param log_d A logical. Should we return the log of the derivative or not?
//'
//' @author David Gerard
//'
//' @return A vector of length 4 containing the gradient of
//'     \code{\link{probgeno}()} with respect to \code{prob}.
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dprobgeno_dprob(const int &gA,
                          const int &gB,
                          const int K,
                          const arma::vec prob) {
  int minz = std::max(0, gA + gB - K);
  int maxz = std::min(gA, gB);

  arma::vec deriv(4);
  deriv.fill(-arma::datum::inf);
  double logdenom = -arma::datum::inf;
  arma::vec x(4);
  for (int z = minz; z <= maxz; z++) {
    x[0] = K + z - gA - gB; // number ab
    x[1] = gA - z; // number Ab
    x[2] = gB - z; // number aB
    x[3] = z; // number AB
    logdenom = log_sum_exp_2(logdenom, dmulti_double(x, prob, true));
    deriv = plog_sum_exp(deriv, dmulti_dprob(x, prob, true));
  }
  deriv = arma::exp(deriv - logdenom);
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
arma::vec dproballgeno_dprob(const arma::vec &gA,
                             const arma::vec &gB,
                             const int K,
                             const arma::vec prob) {
  if (gA.n_elem != gB.n_elem) {
    Rcpp::stop("gA and gB need to be the same length");
  }
  int n = gA.n_elem;

  arma::vec deriv = {0.0, 0.0, 0.0, 0.0};
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
arma::mat dreal_to_simplex_dy(const arma::vec y) {
  int K = y.n_elem + 1;

  arma::mat jacob(K, K - 1);

  double recsum;
  double drecsum;
  double zk;
  for (int j = 0; j < K - 1; j++) {
    recsum = 0.0;
    drecsum = 0.0;
    for (int i = 0; i < K; i++) {
      if (i < K - 1) {
        zk = expit(y[i] + std::log(1.0 / ((double)K - ((double)i + 1))));
      } else{
        zk = expit(std::log(1.0 / ((double)K - ((double)i + 1))));
      }
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

//' Derivative of \code{\link{simplex_to_real}()} with respect to \code{x}.
//'
//' @param x A simplex elements.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::mat dsimplex_to_real_dx(const arma::vec x) {
  int K = x.n_elem;

  arma::mat jacob(K - 1, K);

  double recsum = 0.0;
  double zk = 0.0;
  double multval = 0.0;
  for (int i = 0; i < K - 1; i++) {
    zk = x[i] / (1.0 - recsum);
    multval = 1.0 / (zk * (1.0 - zk));
    for (int j = 0; j < K; j++) {
      if (i == j) {
        jacob(i, j) = 1.0 / (1.0 - recsum);
      } else if (i < j) {
        jacob(i, j) = 0.0;
      } else {
        jacob(i, j) = x[i] / std::pow(1.0 - recsum, 2.0);
      }
      jacob(i, j) = jacob(i, j) * multval;
    }
    recsum += x[i];
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
arma::vec dllike_geno_dpar(const arma::vec par,
                           const arma::vec &gA,
                           const arma::vec &gB,
                           const int K,
                           const arma::vec alpha) {
  if (par.n_elem != 3) {
    Rcpp::stop("par needs to be length 3");
  }

  arma::mat dp_dy = dreal_to_simplex_dy(par);
  arma::vec prob = real_to_simplex(par);
  arma::vec df_dp = dproballgeno_dprob(gA, gB, K, prob) +
    dlprior_dprob(prob, alpha);

  arma::vec deriv = {0.0, 0.0, 0.0};
  for (int i = 0; i < 4; i++) {
    deriv[0] += df_dp[i] * dp_dy(i, 0);
    deriv[1] += df_dp[i] * dp_dy(i, 1);
    deriv[2] += df_dp[i] * dp_dy(i, 2);
  }

  return deriv;
}


//' Derivative of prob[[4]] - (prob[[2]] + prob[[4]]) * (prob[[3]] + prob[[4]])
//' with respect to prob.
//'
//' These derivatives are with repsect to prob[1], prob[2], and prob[3].
//' prob[4] is defined in terms of those other three values as
//' prob[4] = 1 - sum(prob[1:3]).
//'
//' @param prob Probability vector in order of (ab, Ab, aB, AB)
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dD_dprob(const arma::vec prob) {

  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("dD_dprob: prob must sum to 1.");
  }
  if (prob.n_elem != 4) {
    Rcpp::stop("dD_dprob: prob must have 4 elements.");
  }

  arma::vec deriv = {0.0, 0.0, 0.0};

  deriv[0] = prob[1] + prob[2] + 2.0 * prob[3] - 1.0;
  deriv[1] = prob[1] + prob[3] - 1.0;
  deriv[2] = prob[2] + prob[3] - 1.0;

  return deriv;
}

//' Derivative of squared correlation with respect to prob.
//'
//' These derivatives are with repsect to prob[1], prob[2], and prob[3].
//' prob[4] is defined in terms of those other three values as
//' prob[4] = 1 - sum(prob[1:3]).
//'
//' @param prob Probability vector in order of (ab, Ab, aB, AB)
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dr2_dprob(const arma::vec prob) {

  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("dD_dprob: prob must sum to 1.");
  }
  if (prob.n_elem != 4) {
    Rcpp::stop("dD_dprob: prob must have 4 elements.");
  }

  arma::vec deriv = {0.0, 0.0, 0.0};

  double pA = prob[1] + prob[3];
  double pB = prob[2] + prob[3];
  double D = prob[3] - pA * pB;
  double denom = pA * (1.0 - pA) * pB * (1.0 - pB);
  double ddenom_dp1 = -(1.0 - pA) * pB * (1.0 - pB) +
    pA * pB * (1.0 - pB) -
    pA * (1.0 - pA) * (1.0 - pB) +
    pA * (1.0 - pA) * pB;
  double ddenom_dp2 = -pA * (1.0 - pA) * (1.0 - pB) +
    pA * (1.0 - pA) * pB;
  double ddenom_dp3 = -(1.0 - pA) * pB * (1.0 - pB) +
    pA * pB * (1.0 - pB);
  arma::vec dDdp = dD_dprob(prob);

  deriv = 2.0 * D * dDdp / denom;
  deriv[0] = deriv[0] - std::pow(D, 2.0) * ddenom_dp1/ std::pow(denom, 2.0);
  deriv[1] = deriv[1] - std::pow(D, 2.0) * ddenom_dp2/ std::pow(denom, 2.0);
  deriv[2] = deriv[2] - std::pow(D, 2.0) * ddenom_dp3/ std::pow(denom, 2.0);
  return deriv;
}


//' Derivative of D'with respect to prob.
//'
//' These derivatives are with repsect to prob[1], prob[2], and prob[3].
//' prob[4] is defined in terms of those other three values as
//' prob[4] = 1 - sum(prob[1:3]).
//'
//' @param prob Probability vector in order of (ab, Ab, aB, AB)
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dDprime_dprob(const arma::vec prob) {
  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("dD_dprob: prob must sum to 1.");
  }
  if (prob.n_elem != 4) {
    Rcpp::stop("dD_dprob: prob must have 4 elements.");
  }

  arma::vec deriv = {0.0, 0.0, 0.0};

  double pA = prob[1] + prob[3];
  double pB = prob[2] + prob[3];
  double D = prob[3] - pA * pB;
  double denom;
  double ddenom_dp1;
  double ddenom_dp2;
  double ddenom_dp3;

  if ((D < 0.0) & (pA * pB > (1.0 - pA) * (1.0 - pB))) {
    denom = (1.0 - pA) * (1.0 - pB);
    ddenom_dp1 = (1.0 - pB) + (1.0 - pA);
    ddenom_dp2 = (1.0 - pA);
    ddenom_dp3 = (1.0 - pB);
  } else if ((D < 0.0) & (pA * pB < (1.0 - pA) * (1.0 - pB))) {
    denom = pA * pB;
    ddenom_dp1 = -pB - pA;
    ddenom_dp2 = -pA;
    ddenom_dp3 = -pB;
  } else if ((D > 0.0) & (pA * (1.0 - pB) > (1.0 - pA) * pB)) {
    denom = (1.0 - pA) * pB;
    ddenom_dp1 = pB - (1.0 - pA);
    ddenom_dp2 = -(1.0 - pA);
    ddenom_dp3 = pB;
  } else {
    denom = pA * (1.0 - pB);
    ddenom_dp1 = -(1.0 - pB) + pA;
    ddenom_dp2 = pA;
    ddenom_dp3 = -(1.0 - pB);
  }

  arma::vec dDdp = dD_dprob(prob);
  deriv = dDdp / denom;
  deriv[0] = deriv[0] - D * ddenom_dp1 / std::pow(denom, 2.0);
  deriv[1] = deriv[1] - D * ddenom_dp2 / std::pow(denom, 2.0);
  deriv[2] = deriv[2] - D * ddenom_dp3 / std::pow(denom, 2.0);

  return deriv;
}








