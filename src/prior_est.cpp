#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double log_sum_exp_2(double x, double y);
double log_sum_exp(const arma::vec &x);

//' Objective function for \code{\link{em_li}()}
//'
//' @inheritParams em_li
//' @param pivec The current prior probability vector.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double llike_li(const arma::mat &A, const arma::vec &pivec) {
  double ll = 0.0;
  int K = A.n_cols - 1;
  int n = A.n_rows;

  for (int i = 0; i < n; i++) {
    double ival = 0.0;
    for (int j = 0; j <= K; j++) {
      ival += A(i, j) * pivec(j);
    }
    ll += std::log(ival);
  }
  return ll;
}


//' Same as \code{\link{llike_li}()} but using genotype log likelihoods
//'
//' @param B The log-likelihood matrix. Rows are individuals columns are
//'     genotypes.
//' @param lpivec The log prior vector.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double llike_li_log(const arma::mat &B, const arma::vec &lpivec) {

  double ll = 0.0;
  int K = B.n_cols - 1;
  int n = B.n_rows;

  if (B.n_cols != lpivec.n_elem) {
    Rcpp::stop("Number of columns in B should equal length of lpivec");
  }

  for (int i = 0; i < n; i++) {
    double ival = R_NegInf;
    for (int j = 0; j <= K; j++) {
      ival = log_sum_exp_2(ival, B(i, j) + lpivec(j));
    }
    ll += ival;
  }
  return ll;
}

//' EM algorithm from Li (2011)
//'
//' EM algorithm to estimate prior genotype probabilities from genotype
//' likelihoods.
//'
//' @param A The genotype likelihood (not log-likelihood) matrix. The rows
//'     index the individuals and the columns index the genotypes. So
//'     \code{A[i, k]} is the genotype likelihood for genotype \code{k} and
//'     individual \code{i}.
//' @param itermax The maximum number of iterations.
//' @param eps The stopping criteria.
//'
//' @return A vector of prior probabilities for each genotype.
//'
//' @references
//' \itemize{
//'   \item{Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. \emph{Bioinformatics}, 27(21), 2987-2993. \doi{10.1093/bioinformatics/btr509}}
//' }
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec em_li(const arma::mat &A, int itermax = 100, double eps = 1e-5) {
  int K = A.n_cols - 1;
  int n = A.n_rows;
  double valinit = 1 / ((double)K + 1);

  arma::vec pivec(K + 1, arma::fill::value(valinit));

  double llold;
  double llnew = llike_li(A, pivec);
  double err = R_PosInf;
  arma::vec w(K + 1);

  int i = 0;
  while((i < itermax) & (err > eps)) {
    llold = llnew;
    w.fill(0.0);

    for (int j = 0; j < n; j++) {
      arma::vec api(K + 1);
      for (int k = 0; k <= K; k++) {
        api(k) = A(j, k) * pivec(k);
      }
      api = api / sum(api);
      for (int k = 0; k <= K; k++) {
        w(k) += api(k);
      }
    }

    pivec = w / arma::sum(w);
    llnew = llike_li(A, pivec);
    if (llnew < llold) {
      Rcpp::stop("log-likelihood not increasing");
    }

    err = std::abs(llnew - llold);
    i++;
  }
  return pivec;
}

//' Same as \code{\link{em_li}()} but using genotype log-likelihoods
//'
//' @param B Matrix of genotype log-likelihoods. The rows index the individuals
//'     and the columsn index the genotypes.
//' @inheritParams em_li
//'
//' @return A vector of log prior probabilities for each genotype.
//'
//' @author David Gerard
//'
//' @references
//' \itemize{
//'   \item{Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. \emph{Bioinformatics}, 27(21), 2987-2993. \doi{10.1093/bioinformatics/btr509}}
//' }
//'
//' @noRd
// [[Rcpp::export]]
arma::vec em_li_log(const arma::mat &B, int itermax = 100, double eps = 1e-5) {
  int K = B.n_cols - 1;
  int n = B.n_rows;
  double valinit = -std::log((double)K + 1);

  arma::vec lpivec(K + 1, arma::fill::value(valinit));

  double llold;
  double llnew = llike_li_log(B, lpivec);
  double err = R_PosInf;
  arma::vec lw(K + 1);

  int i = 0;
  while((i < itermax) & (err > eps)) {
    llold = llnew;
    lw.fill(-arma::datum::inf);

    for (int j = 0; j < n; j++) {
      arma::vec api(K + 1);
      for (int k = 0; k <= K; k++) {
        api(k) = B(j, k) + lpivec(k);
      }
      api = api - log_sum_exp(api);
      for (int k = 0; k <= K; k++) {
        lw(k) = log_sum_exp_2(lw(k), api(k));
      }
    }

    lpivec = lw - log_sum_exp(lw);
    llnew = llike_li_log(B, lpivec);
    if (llnew < llold) {
      Rcpp::stop("log-likelihood not increasing");
    }

    err = std::abs(llnew - llold);
    i++;
  }

  return lpivec;
}

