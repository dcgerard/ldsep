#include <RcppArmadillo.h>
using namespace Rcpp;

const double TOL = std::sqrt(DOUBLE_EPS);

//' Log-sum-exponential trick.
//'
//' @param x A vector to log-sum-exp.
//'
//' @return The log of the sum of the exponential
//'     of the elements in \code{x}.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double log_sum_exp(const arma::vec &x) {
  double max_x = arma::max(x);
  double lse; // the log-sum-exp
  // if all -Inf, need to treat this special to avoid -Inf + Inf.
  if (max_x == -arma::datum::inf) {
    lse = -arma::datum::inf;
  } else {
    lse = max_x + std::log(arma::sum(arma::exp(x - max_x)));
  }
  return lse;
}

//' log-sum-expontential of a matrix.
//'
//' @noRd
// [[Rcpp::export]]
double log_sum_exp_mat(const arma::mat &x) {
  double max_x = x.max();
  double lse; // the log-sum-exp
  // if all -Inf, need to treat this special to avoid -Inf + Inf.
  if (max_x == -arma::datum::inf) {
    lse = -arma::datum::inf;
  } else {
    lse = max_x + std::log(arma::accu(arma::exp(x - max_x)));
  }
  return lse;
}

//' Log-sum-exponential trick using just two doubles.
//'
//' @param x A double.
//' @param y Another double.
//'
//' @return The log of the sum of the exponential of x and y.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double log_sum_exp_2(double x, double y) {
  double z = std::max(x, y);
  double finalval;
  if (z == -arma::datum::inf) {
    finalval = -arma::datum::inf;
  } else {
    finalval = std::log(std::exp(x - z) + std::exp(y - z)) + z;
  }
  return finalval;
}

//' Pair-wise log-sum-exponential
//'
//' Does pair-wise log-sum-exponential on two vectors.
//'
//' @param x A numeric vector.
//' @param y Another numeric vector.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec plog_sum_exp(const arma::vec &x,
                       const arma::vec &y) {

  if (x.n_elem != y.n_elem) {
    Rcpp::stop("x and y must have the same length");
  }

  int n = x.n_elem;

  arma::vec z(n);

  for (int i = 0; i < n; i++) {
    z[i] = log_sum_exp_2(x[i], y[i]);
  }

  return z;
}

//' Parallel log-sum-exp of two matrices
//'
//' @noRd
// [[Rcpp::export]]
arma::mat plog_sum_exp_mat(const arma::mat &x,
                           const arma::mat &y) {

  if (x.n_rows != y.n_rows) {
    Rcpp::stop("plog_sum_exp_mat: x and y must have the same number of rows");
  }
  if (x.n_cols != y.n_cols) {
    Rcpp::stop("plog_sum_exp_mat: x and y must have the same number of columns");
  }

  int nrow = x.n_rows;
  int ncol = x.n_cols;

  arma::mat z(nrow, ncol);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      z(i, j) = log_sum_exp_2(x(i, j), y(i, j));
    }
  }

  return z;
}

//' The logit function.
//'
//' @param x A double between 0 and 1.
//'
//' @return The logit of \code{x}.
//'
//' @author David Gerard
//'
//' @keywords internal
//' @noRd
//'
// [[Rcpp::export]]
double logit(double x) {
  if ((x < TOL) | ((1.0 - x) < TOL)) {
    Rcpp::stop("logit: x must be between 0 and 1.");
  }
  double lv = std::log(x / (1.0 - x));
  return lv;
}

//' The expit (logistic) function.
//'
//' @param x A double.
//'
//' @return The expit (logistic) of \code{x}.
//'
//' @keywords internal
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
double expit(double x) {
  double ev = 1.0 / (1.0 + std::exp(-x));
  return ev;
}


//' Convert real line to simplex using Stan technique
//'
//' @param y A vector of numbers of length K-1
//'
//' @return A vector on the simplex of length K
//'
//' @author David Gerard
//'
//' @seealso \code{\link{simplex_to_real}()} for the inverse function.
//'
//' @noRd
// [[Rcpp::export]]
arma::vec real_to_simplex(const arma::vec y) {
  int K = y.n_elem + 1;
  arma::vec x(K);

  double recsum = 0.0;
  for (int k = 0; k < K - 1; k++) {
    x[k] = (1.0 - recsum) * expit(y[k] + std::log(1.0 / ((double)K - ((double)k + 1))));
    recsum += x[k];
  }

  x[K - 1] = 1.0 - recsum;

  return x;
}

//' Convert simplex to real-line using Stan technique
//'
//' @param x A vector on the simplex of length K.
//'
//' @return A vector of numbers of length K-1
//'
//' @author David Gerard
//'
//' @seealso \code{\link{real_to_simplex}()} for the inverse function.
//'
//' @noRd
// [[Rcpp::export]]
arma::vec simplex_to_real(const arma::vec x) {
  if (std::abs(arma::sum(x) - 1.0) > TOL) {
    Rcpp::stop("x should sum to 1");
  }

  int K = x.n_elem;
  arma::vec y(K - 1);

  double recsum = 0.0;
  for (int k = 0; k < K - 1; k++) {
    y[k] = logit(x[k] / (1.0 - recsum)) - std::log(1 / ((double)K - ((double)k + 1)));
    recsum += x[k];
  }

  return y;
}















