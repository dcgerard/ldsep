#include <RcppArmadillo.h>
using namespace Rcpp;

double log_sum_exp_mat(const arma::mat &x);
double log_sum_exp_2(double x, double y);

//' Density function of the multivariate normal distribution
//'
//' @param x A vector containing the quantile.
//' @param mu A vector containing the mean.
//' @param sigma A positive definite covariance matrix
//' @param log A logical. If \code{TRUE}, log probabilities are returned.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double dmvnorm(arma::vec x, arma::vec mu, arma::mat sigma, bool log = false) {
  if (x.n_elem != mu.n_elem) {
    Rcpp::stop("dmvnorm: x and mu must have same length.");
  }
  if (!sigma.is_sympd()) {
    Rcpp::Rcout << "Sigma:"
                << std::endl
                << sigma
                << std::endl;
    Rcpp::stop("dmvnorm: sigma must be symmetic positive definite.");
  }
  if (x.n_elem != sigma.n_cols) {
    Rcpp::stop("Sigma should have the same number of columns/rows as the length of x.");
  }

  int k = x.n_elem;
  arma::vec resid = x - mu;
  arma::mat cp = 0.5 * resid.t() * arma::inv_sympd(sigma) * resid;
  std::complex<double> ldet = arma::log_det(sigma);

  double lval = -(double)k * 0.5 * std::log(2.0 * arma::datum::pi) -
    0.5 * ldet.real() -
    cp(0,0);

  if (!log) {
    lval = std::exp(lval);
  }

  return lval;
}


//' Returns distribution of proportional bivariate normal.
//'
//' @param mu A vector of length 2 containing the mean.
//' @param sigma A 2-by-2 positive definite covariance matrix
//' @param K The ploidy of the individual.
//' @param log A logical. If \code{TRUE}, log probabilities are returned.
//'
//' @return A matrix. Element (i,j) is the (log) probability of genotype
//'     i-1 at locus 1 and j-1 at locus 2.
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat pbnorm_dist(arma::vec mu, arma::mat sigma, int K, bool log = false) {
  if (mu.n_elem != 2) {
    Rcpp::stop("dpbnorm: mu should have length 2.");
  }
  if ((sigma.n_rows != 2) | (sigma.n_cols != 2)) {
    Rcpp::stop("dpbnorm: sigma should be a 2-by-2 matrix.");
  }
  if (!sigma.is_sympd()) {
    Rcpp::Rcout << "Sigma:"
                << std::endl
                << sigma
                << std::endl;
    Rcpp::stop("pbnorm_dist: sigma must be symmetic positive definite.");
  }

  arma::mat distmat(K + 1, K + 1);
  arma::vec x(2);
  for (int i = 0; i <= K; i++) {
    for (int j = 0; j <= K; j++) {
      x(0) = (double)i;
      x(1) = (double)j;
      distmat(i, j) = dmvnorm(x, mu, sigma, true);
    }
  }

  distmat = distmat - log_sum_exp_mat(distmat);

  if (!log) {
    distmat = arma::exp(distmat);
  }

  return distmat;
}


//' Log-likelihood for joint genotype distribution when using a proportional
//' normal model.
//'
//' @inheritParams pbnorm_dist
//' @param pgA The matrix of genotype log-likelihoods for locus 1.
//'     The rows index the individuals and the columns index the genotypes.
//' @param pgB The matrix of genotype log-likelihoods for locus 2.
//'     The rows index the individuals and the columns index the genotypes.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double llike_pbnorm_genolike(const arma::mat &pgA,
                             const arma::mat &pgB,
                             const arma::vec &mu,
                             const arma::mat &sigma) {
  if (pgA.n_rows != pgB.n_rows) {
    Rcpp::stop("llike_pbnorm_genolike: dimensions of pgA and pgB are different");
  }
  if (pgA.n_cols != pgB.n_cols) {
    Rcpp::stop("llike_pbnorm_genolike: dimensions of pgA and pgB are different");
  }
  int K = pgA.n_cols - 1;
  int n = pgA.n_rows;
  arma::mat distmat = pbnorm_dist(mu, sigma, K, true);

  double lval = 0.0;
  double li;
  for (int ell = 0; ell < n; ell++) {
    li = -arma::datum::inf;
    for (int i = 0; i <= K; i++) {
      for (int j = 0; j <= K; j++) {
        li = log_sum_exp_2(li, pgA(ell, i) + pgB(ell, j) + distmat(i, j));
      }
    }
    lval += li;
  }

  return lval;
}

// [[Rcpp::export]]
double prior_mu(arma::vec mu, int K) {
  if (mu.n_elem != 2) {
    Rcpp::stop("prior_mu: mu not of length 2.");
  }
  double lval = R::dnorm(mu(0), (double)K / 2.0, (double)K, true) +
    R::dnorm(mu(1), (double)K / 2.0, (double)K, true);
  return lval;
}

// [[Rcpp::export]]
double prior_sigma(arma::vec lvec) {
  if (lvec.n_elem != 3) {
    Rcpp::stop("prior_sigma: lvec not of length 3.");
  }
  double cval = 1.5;
  double lval = R::dchisq(lvec(0) / cval, 5.0, true) + std::log(2.0 * lvec(0) / cval) +
    R::dnorm(lvec(1) / cval, 0.0, 1.0, true) +
    R::dchisq(lvec(2) / cval, 4.0, true) + std::log(2.0 * lvec(2) / cval) -
    3.0 * std::log(cval);
  return lval;
}

//' Wrapper for \code{\link{llike_pbnorm_genolike}()} that tkaes
//' a vector of parameters as input for optim().
//'
//' Also includes prior probs
//'
//' @inheritParams llike_pbnorm_genolike
//' @param par A vector of length 5. The first two elements are \code{mu}. The
//'     last three elements are c(l11, l12, l22), the lower three elements of
//'     the cholesky decomposition of sigma.
//'
//' @author David Gerard
//'
//' @noRd
//'
// [[Rcpp::export]]
double obj_pbnorm_genolike(const arma::vec &par,
                           const arma::mat &pgA,
                           const arma::mat &pgB) {
  arma::vec mu(2);
  arma::mat sigma(2, 2);
  double lval;
  int K = pgA.n_cols - 1;

  mu(0) = par(0);
  mu(1) = par(1);
  sigma(0, 0) = par(2);
  sigma(1, 0) = par(3);
  sigma(0, 1) = 0.0;
  sigma(1, 1) = par(4);

  sigma = sigma * sigma.t();

  lval = llike_pbnorm_genolike(pgA, pgB, mu, sigma) +
    prior_mu(mu, K) +
    prior_sigma(par.tail(3));

  return lval;
}


//' Derivative of \code{\link{dmvnorm}()} (not log) with respect to \code{R}
//' where \code{sigma = solve(R) \%*\% solve(t(R))}
//'
//'
//' @inheritParams dmvnorm
//' @param R The inverse of the cholesky square root of sigma.
//'     \code{sigma = solve(R) \%*\% solve(t(R))}
//'
//' @noRd
// arma::mat dmvnorm_dR(arma::vec x, arma::vec mu, arma::mat R) {
//   if (!R.is_trimatl()) {
//     Rcpp::stop("dmvnorm_dR: R needs to be lower triangular.");
//   }
//   if (R.n_cols != R.n_rows) {
//     Rcpp::stop("dmvnorm_dR: R needs to be square.");
//   }
//   if (x.n_elem != mu.n_elem) {
//     Rcpp::stop("dmvnorm_dR: x and mu should have same length.");
//   }
//   if (x.n_elem != R.n_cols) {
//     Rcpp::stop("dmvnorm_dR: R should have same number of columns as length of x.");
//   }
//
//   int k = R.n_cols;
//   arma::mat dR(k, k);
//   arma::vec resid = x - mu;
//   arma::mat ss = resid.t() * R.t() * R * resid;
//
//   2.0 * arma::datum::pi *
//     std::exp(-0.5 * ss(0, 0)) *
//     ;
//
//   return dR;
// }

