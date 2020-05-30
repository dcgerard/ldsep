#include <RcppArmadillo.h>
using namespace Rcpp;

const double TOL = std::sqrt(DOUBLE_EPS);
double dmulti_double(const arma::vec x,
                     const arma::vec prob,
                     bool log_p = true);
double log_sum_exp(const arma::vec &x);
double log_sum_exp_2(double x, double y);
arma::vec plog_sum_exp(const arma::vec &x,
                       const arma::vec &y);



//' Get a matrix of all possible haplotype numbers
//'
//' The numbers are returned in lexicographical order.
//'
//' @param K The ploidy of the species.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::mat get_Amat(int K) {
 int numa = R::choose(K + 4 - 1, K); // K indistinguishable balls into 4 distinguishable urns. "stars and bars"
 arma::mat Amat(4, numa);

 int i = 0;
 for (int a1 = K; a1 >= 0; a1--) {
   for (int a2 = K - a1; a2 >= 0; a2--) {
     for (int a3 = K - a1 - a2; a3 >= 0; a3--) {
       Amat(0, i) = a1;
       Amat(1, i) = a2;
       Amat(2, i) = a3;
       Amat(3, i) = K - a1 - a2 - a3;
       i++;
     }
   }
 }

 return Amat;
}

//' EM algorithm to estimate haplotype frequencies
//'
//' This runs an EM algorithm to obtain the maximum likelihood estimates
//' of the haplotype frequencies for two loci when one has access
//' to genotype likelihoods.
//'
//' @param p A vector of length 4. The intialization for the
//'     haplotype frequencies.
//' @param pgA The matrix of genotype log-likelihoods for locus 1.
//'     The rows index the individuals and the columns index the genotypes.
//' @param pgB The matrix of genotype log-likelihoods for locus 2.
//'     The rows index the individuals and the columns index the genotypes.
//' @param alpha The prior sample size used in the penalty.
//' @param maxit The maximum number of EM iterations.
//' @param tol The convergence tolerance.
//' @param verbose Should we output more (\code{TRUE}) or less
//'     (\code{FALSE}).
//'
//'
//' @author David Gerard
//'
//' @export
//'
// [[Rcpp::export]]
arma::vec genolike_em(arma::vec p,
                       const arma::mat &pgA,
                       const arma::mat &pgB,
                       const arma::vec &alpha,
                       const int maxit = 100,
                       const double tol = 0.001,
                       bool verbose = false) {

  if (pgA.n_rows != pgB.n_rows) {
    Rcpp::stop("genolike_em: dimensions of pgA and pgB are different");
  }
  if (pgA.n_cols != pgB.n_cols) {
    Rcpp::stop("genolike_em: dimensions of pgA and pgB are different");
  }
  if (std::abs(arma::sum(p) - 1.0) > TOL) {
    Rcpp::stop("genolike_em: p should sum to 1");
  }
  if (p.n_elem != 4) {
    Rcpp::stop("genolike_em: p should have exactly 4 elements");
  }
  if (alpha.n_elem != 4) {
    Rcpp::stop("genolike_em: alpha should have exactly 4 elements");
  }

  int K = pgA.n_cols - 1; // ploidy
  int n = pgA.n_rows; // number of individuals

  arma::mat Amat = get_Amat(K); // possible "a" values
  int na = Amat.n_cols; // number of possible "a" values

  arma::mat Wmat(na, n); // log of EM weights
  arma::vec multivec(na); // multinomial log-probabilities

  int gA; // current A genotype
  int gB; // current B genotype

  arma::vec pold(4);

  double err = tol + 1.0;
  int iternum = 0;
  while ((err > tol) & (iternum < maxit)) {
    pold = p;

    // get multinomial probabilities
    for (int a = 0; a < na; a++) {
      multivec(a) = dmulti_double(Amat.col(a), p, true);
    }

    // E-step
    for(int i = 0; i < n; i++) {
      for (int a = 0; a < na; a++) {
        gA = Amat(1, a) + Amat(3, a);
        gB = Amat(2, a) + Amat(3, a);
        Wmat(a, i) = pgA(i, (int)gA) +
          pgB(i, (int)gB) +
          multivec(a);
      }
      Wmat.col(i) = Wmat.col(i) - log_sum_exp(Wmat.col(i));
    }
    // M-step
    p.fill(-arma::datum::inf);
    for(int i = 0; i < n; i++) {
      for (int a = 0; a < na; a++) {
        if (Amat(0, a) > 0.5) {
          p(0) = log_sum_exp_2(p(0), Wmat(a, i) + std::log(Amat(0, a)));
        }
        if (Amat(1, a) > 0.5) {
          p(1) = log_sum_exp_2(p(1), Wmat(a, i) + std::log(Amat(1, a)));
        }
        if (Amat(2, a) > 0.5) {
          p(2) = log_sum_exp_2(p(2), Wmat(a, i) + std::log(Amat(2, a)));
        }
        if (Amat(3, a) > 0.5) {
          p(3) = log_sum_exp_2(p(3), Wmat(a, i) + std::log(Amat(3, a)));
        }
      }
    }
    if ((alpha(0) > 1.0) & (alpha(1) > 1.0) & (alpha(2) > 1.0) & (alpha(3) > 1.0)) {
      p = plog_sum_exp(p, arma::log(alpha - 1.0));
      p = arma::exp(p - log_sum_exp(p));
    } else {
      p = arma::exp(p) + alpha - 1.0;
      p = p / arma::sum(p);
    }

    // check convergence
    iternum++;
    err = arma::max(arma::abs(p - pold));

    // print
    if (verbose) {
      Rcpp::Rcout << "Iteration: "
                  << iternum
                  << std::endl
                  << "Criteria: "
                  << err
                  << std::endl
                  << "p: "
                  << p.t()
                  << std::endl
                  << std::endl;
    }
  }

  return p;
}
