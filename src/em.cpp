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


arma::vec em_fix(const arma::vec &p,
                 arma::vec &multivec,
                 const arma::mat &Amat,
                 const arma::mat &pgA,
                 const arma::mat &pgB,
                 arma::mat &Wmat,
                 const arma::vec &alpha) {
  int n = pgA.n_rows; // number of individuals
  int na = Amat.n_cols; // number of possible "a" values
  int gA; // current A genotype
  int gB; // current B genotype

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
  arma::vec pnew(4);
  pnew.fill(-arma::datum::inf);
  for(int i = 0; i < n; i++) {
    for (int a = 0; a < na; a++) {
      if (Amat(0, a) > 0.5) {
        pnew(0) = log_sum_exp_2(pnew(0), Wmat(a, i) + std::log(Amat(0, a)));
      }
      if (Amat(1, a) > 0.5) {
        pnew(1) = log_sum_exp_2(pnew(1), Wmat(a, i) + std::log(Amat(1, a)));
      }
      if (Amat(2, a) > 0.5) {
        pnew(2) = log_sum_exp_2(pnew(2), Wmat(a, i) + std::log(Amat(2, a)));
      }
      if (Amat(3, a) > 0.5) {
        pnew(3) = log_sum_exp_2(pnew(3), Wmat(a, i) + std::log(Amat(3, a)));
      }
    }
  }
  if ((alpha(0) > 1.0) & (alpha(1) > 1.0) & (alpha(2) > 1.0) & (alpha(3) > 1.0)) {
    pnew = plog_sum_exp(pnew, arma::log(alpha - 1.0));
    pnew = arma::exp(pnew - log_sum_exp(pnew));
  } else {
    pnew = arma::exp(pnew) + alpha - 1.0;
    pnew = pnew / arma::sum(pnew);
  }

  return pnew;
}


// Project onto simplex via algorithm from Chen and Ye (2011)
// "Projection Onto A Simplex"
// [[Rcpp::export]]
arma::vec simplex_proj(arma::vec y) {
  double that;
  arma::vec y_ord = arma::sort(y);
  for (int i = 3; i >= 0; i--) {
    if (i > 0) {
      that = (arma::sum(y_ord.tail(4 - i)) - 1.0) / (4.0 - (double)i);
      if (that >= y_ord(i - 1)) {
        break;
      }
    } else {
      that = (arma::sum(y) - 1.0) / 4.0;
    }
  }
  for (int i = 0; i < 4; i++) {
    y(i) = y(i) - that;
    if (y(i) < 0.0) {
      y(i) = 0.0;
    }
  }
  return y;
}

//' EM algorithm to estimate haplotype frequencies
//'
//' This runs an EM algorithm to obtain the maximum likelihood estimates
//' of the haplotype frequencies for two loci when one has access
//' to genotype likelihoods.
//'
//' EM Squaring is performed via the algorithm of Ravi and Roland (2008)
//' with each iteration's projection onto the simplex using the algorithm
//' of Chen and Ye (2011). Though squaring is turned off be default because
//' it usually doesn't help.
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
//'     (\code{FALSE})?
//' @param square Should we implement squared acceleratred EM (\code{TRUE})
//'     or not (\code{FALSE})?
//'
//'
//' @references
//' \itemize{
//'   \item{Varadhan, Ravi, and Christophe Roland. "Simple and globally convergent methods for accelerating the convergence of any EM algorithm." Scandinavian Journal of Statistics 35.2 (2008): 335-353.}
//'   \item{Chen, Yunmei, and Xiaojing Ye. "Projection onto a simplex." arXiv preprint arXiv:1101.6081 (2011).}
//' }
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
                       const int maxit = 500,
                       const double tol = 0.0001,
                       bool verbose = false,
                       bool square = false) {

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

  arma::vec pold(4);

  // squarem parameters ---
  arma::vec p1(4);
  arma::vec p2(4);
  arma::vec pprime(4);
  arma::vec r(4);
  arma::vec v(4);
  double steplength;
  // end squarem parameters ---

  double err = tol + 1.0;
  int iternum = 0;
  while ((err > tol) & (iternum < maxit)) {
    pold = p;

    if (!square) {
      p = em_fix(p, multivec, Amat, pgA, pgB, Wmat, alpha);
    } else {
      // calculate new point via SQUAREM
      p1 = em_fix(p, multivec, Amat, pgA, pgB, Wmat, alpha);
      p2 = em_fix(p1, multivec, Amat, pgA, pgB, Wmat, alpha);
      r = p1 - p;
      v = p2 + p;
      // steplength = -1.0 * arma::norm(r) / arma::norm(v);
      steplength = arma::dot(r, r) / arma::dot(r, v);
      // steplength = arma::dot(r, v) / arma::dot(v, v);

      pprime = p - 2 * steplength * r + std::pow(steplength, 2.0) * v;
      pprime = simplex_proj(pprime);
      // pprime = pprime / arma::sum(pprime);
      p = em_fix(pprime, multivec, Amat, pgA, pgB, Wmat, alpha);
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
