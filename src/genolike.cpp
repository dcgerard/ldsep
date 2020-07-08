// all functions for dealing with genotype likelihoods

#include <RcppArmadillo.h>
using namespace Rcpp;

const double TOL = std::sqrt(DOUBLE_EPS);
double log_sum_exp_2(double x, double y);
double probgeno(const int &gA,
                const int &gB,
                const int K,
                const arma::vec prob,
                bool log_p);
arma::vec dprobgeno_dprob(const int &gA,
                          const int &gB,
                          const int K,
                          const arma::vec prob);
arma::vec real_to_simplex(const arma::vec y);
arma::mat dreal_to_simplex_dy(const arma::vec y);
arma::vec dmulti_dprob(const arma::vec x,
                       const arma::vec prob,
                       bool log_p);
arma::vec plog_sum_exp(const arma::vec &x,
                       const arma::vec &y);
arma::mat get_prob_array(int K, arma::vec prob, bool log_p);
double lprior(const arma::vec prob, const arma::vec alpha);
arma::vec dlprior_dprob(const arma::vec prob, const arma::vec alpha);

//' Probability of genotype likelihoods given haplotype frequencies for a
//' single individual.
//'
//' The ploidy of the species is assumed to be one less the length of
//' \code{pgA} and \code{pgB} (which must be the same length).
//'
//' @param pgA The genotype log-likelihoods at locus 1. pgA[i] is the
//'     log probability of the data given the genotype at locus 1 is i.
//' @param pgB The genotype log-likelihoods at locus 2. pgA[i] is the
//'     log probability of the data given the genotype at locus 2 is i.
//' @param prob The vector of probabilities for haplotypes (ab, Ab, aB, AB).
//' @param log_p A logical. Should we return the log probability or not?
//'
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double probgenolike(const arma::vec &pgA,
                    const arma::vec &pgB,
                    const arma::vec prob,
                    bool log_p = true) {
  if (pgA.n_elem != pgB.n_elem) {
    Rcpp::stop("probgenolike: pgA and pgB should have the same number of elements");
  }
  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("probgenolike: prob should sum to 1");
  }

  int K = pgA.n_elem - 1;
  double lp = -arma::datum::inf;

  for (int i = 0; i <= K; i++) {
    for (int j = 0; j <= K; j++) {
      lp = log_sum_exp_2(lp, pgA[i] + pgB[j] + probgeno(i, j, K, prob, true));
    }
  }

  if (!log_p) {
    lp = std::exp(lp);
  }

  return lp;
}


//' Probability of genotype likelihoods given haplotype frequencies for all
//' individuals.
//'
//' @param pgA The matrix of genotype log-likelihoods for locus 1.
//'     The rows index the individuals and the columns index the genotypes.
//' @param pgA The matrix of genotype log-likelihoods for locus 2.
//'     The rows index the individuals and the columns index the genotypes.
//' @param prob The vector of probabilities for haplotypes (ab, Ab, aB, AB).
//' @param log_p A logical. Should we return the log probability or not?
//'
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double proballgenolike_old(const arma::mat &pgA,
                           const arma::mat &pgB,
                           const arma::vec prob,
                           bool log_p = true) {
  if (pgA.n_rows != pgB.n_rows) {
    Rcpp::stop("proballgenolike_old: dimensions of pgA and pgB are different");
  }
  if (pgA.n_cols != pgB.n_cols) {
    Rcpp::stop("proballgenolike_old: dimensions of pgA and pgB are different");
  }
  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("proballgenolike_old: prob should sum to 1");
  }
  if (prob.n_elem != 4) {
    Rcpp::stop("proballgenolike_old: prob should have exactly 4 elements");
  }


  int n = pgA.n_rows;
  double lp = 0.0;
  for (int i = 0; i < n; i++) {
    lp += probgenolike(pgA.row(i).t(), pgB.row(i).t(), prob, log_p = true);
  }

  if (!log_p) {
    lp = std::exp(lp);
  }

  return lp;
}

//' Probability of genotype likelihoods given haplotype frequencies for all
//' individuals.
//'
//' @param pgA The matrix of genotype log-likelihoods for locus 1.
//'     The rows index the individuals and the columns index the genotypes.
//' @param pgA The matrix of genotype log-likelihoods for locus 2.
//'     The rows index the individuals and the columns index the genotypes.
//' @param prob The vector of probabilities for haplotypes (ab, Ab, aB, AB).
//' @param log_p A logical. Should we return the log probability or not?
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double proballgenolike(const arma::mat &pgA,
                       const arma::mat &pgB,
                       const arma::vec prob,
                       bool log_p = true) {
  if (pgA.n_rows != pgB.n_rows) {
    Rcpp::stop("proballgenolike_new: dimensions of pgA and pgB are different");
  }
  if (pgA.n_cols != pgB.n_cols) {
    Rcpp::stop("proballgenolike_new: dimensions of pgA and pgB are different");
  }
  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("proballgenolike_new: prob should sum to 1");
  }
  if (prob.n_elem != 4) {
    Rcpp::stop("proballgenolike_new: prob should have exactly 4 elements");
  }

  int n = pgA.n_rows;
  int K = pgA.n_cols - 1;

  arma::mat parray = get_prob_array(K, prob, true);

  double logdenom_ind;
  double logdenom = 0.0;

  for (int i = 0; i < n; i++) {
    logdenom_ind = -arma::datum::inf;
    for (int gA = 0; gA <= K; gA++) {
      for (int gB = 0; gB <= K; gB++) {
        logdenom_ind = log_sum_exp_2(logdenom_ind,
                                     pgA(i, gA) + pgB(i, gB) + parray(gA, gB));
      }
    }
    logdenom += logdenom_ind;
  }

  if (!log_p) {
    logdenom = std::exp(logdenom);
  }

  return logdenom;
}

//' Likelihood function when estimating LD from genotype log-likelihoods
//'
//' @param par A vector of length 3 containing real numbers that are to
//'     be transformed into the simplex of prob (ab, Ab, aB, AB).
//' @param alpha The prior sample size used in the penalty.
//' @inheritParams proballgenolike
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double llike_genolike(const arma::vec par,
                      const arma::mat &pgA,
                      const arma::mat &pgB,
                      const arma::vec alpha) {
  if (par.n_elem != 3) {
    Rcpp::stop("llike_genolike: par needs to be length 3");
  }

  arma::vec prob = real_to_simplex(par);

  double llike = proballgenolike(pgA, pgB, prob, true) +
    lprior(prob, alpha);

  return llike;
}

//' Obtain a matrix of derivatives of p(geno) (NOT log(p(geno)))
//' with respect to prob for all genotypes.
//'
//'
//' @param K the ploidy
//' @param prob Haplotype frequencies in order (ab, Ab, aB, AB).
//'
//' @return Element (i,j,k) is the derivative of \code{\link{probgeno}()}
//'     when gA = i, gB = j with respect to prob[k]
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::cube get_dprobgeno_dprob_array(int K, arma::vec prob) {
  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("get_dprobgeno_dprob_array: prob should sum to 1");
  }
  if (prob.n_elem != 4) {
    Rcpp::stop("get_dprobgeno_dprob_array: prob should have exactly 4 elements");
  }

  int minz;
  int maxz;

  arma::cube derivarray(K + 1, K + 1, 4);
  derivarray.fill(-arma::datum::inf);
  arma::vec x(4);
  for (int gA = 0; gA <= K; gA++) {
    for (int gB = 0; gB <= K; gB++) {
      minz = std::max(0, gA + gB - K);
      maxz = std::min(gA, gB);
      for (int z = minz; z <= maxz; z++) {
        x[0] = K + z - gA - gB; // number ab
        x[1] = gA - z; // number Ab
        x[2] = gB - z; // number aB
        x[3] = z; // number AB
        derivarray.tube(gA, gB) =
          plog_sum_exp(derivarray.tube(gA, gB), dmulti_dprob(x, prob, true));
      }
    }
  }

  return derivarray;
}

//' Obtain the distribution of genotypes given haplotype frequencies under HWE
//'
//' This function will calculate the (log) probabilities for all genotype
//' combinations at two loci given just the haplotype frequencies. This
//' is under the assumptions of HWE.
//'
//' @param K The ploidy of the species.
//' @param prob Haplotype frequencies in the order of (ab, Ab, aB, AB).
//' @param log_p A logical. Should we return the log-probabilities (\code{TRUE})
//'     or the probabilities (\code{FALSE}). Defaults to \code{TRUE}.
//'
//' @return Element (i, j) is the (log) probability of genotype i-1 at locus 1
//'     and genotype j-1 at locus 2.
//'
//' @author David Gerard
//'
//' @examples
//' get_prob_array(K = 6, prob = c(0.1, 0.2, 0.3, 0.4), log_p = FALSE)
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat get_prob_array(int K, arma::vec prob, bool log_p = true) {

  arma::mat probmat(K + 1, K + 1);

  for (int gA = 0; gA <= K; gA++) {
    for (int gB = 0; gB <= K; gB++) {
      probmat(gA, gB) = probgeno(gA, gB, K, prob, true);
    }
  }

  if (!log_p) {
    probmat = arma::exp(probmat);
  }

  return probmat;
}



//' Gradient of \code{\link{probgenolike}()} with respect to \code{prob}.
//'
//' @inheritParams probgenolike
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dprobgenolike_dprob(const arma::vec &pgA,
                              const arma::vec &pgB,
                              const arma::vec prob) {
  int K = pgA.n_elem - 1;
  if (pgA.n_elem != pgB.n_elem) {
    Rcpp::stop("dprobgenolike_dprob: pgA and pgB should have the same number of elements.");
  }

  arma::cube darray = get_dprobgeno_dprob_array(K, prob);
  arma::mat parray = get_prob_array(K, prob, true);

  arma::vec deriv(4);
  deriv.fill(-arma::datum::inf);
  double logdenom = -arma::datum::inf;
  for (int i = 0; i <= K; i++) {
    for (int j = 0; j <= K; j++) {
      deriv = plog_sum_exp(deriv, pgA[i] + pgB[j] + darray.tube(i, j));
      logdenom = log_sum_exp_2(logdenom, pgA[i] + pgB[j] + parray(i, j));
    }
  }

  deriv = arma::exp(deriv - logdenom);

  return deriv;
}


//' Gradient of \code{\link{proballgenolike}(,log = TRUE)} with respect to
//' \code{prob}
//'
//' @inheritParams proballgenolike
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dproballgenolike_dprob(const arma::mat &pgA,
                                 const arma::mat &pgB,
                                 const arma::vec prob) {
  if (pgA.n_rows != pgB.n_rows) {
    Rcpp::stop("dproballgenolike_dprob: dimensions of pgA and pgB are different");
  }
  if (pgA.n_cols != pgB.n_cols) {
    Rcpp::stop("dproballgenolike_dprob: dimensions of pgA and pgB are different");
  }
  if (std::abs(arma::sum(prob) - 1.0) > TOL) {
    Rcpp::stop("dproballgenolike_dprob: prob should sum to 1");
  }
  if (prob.n_elem != 4) {
    Rcpp::stop("dproballgenolike_dprob: prob should have exactly 4 elements");
  }

  int n = pgA.n_rows;
  int K = pgA.n_cols - 1;

  arma::cube darray = get_dprobgeno_dprob_array(K, prob);
  arma::mat parray = get_prob_array(K, prob, true);

  arma::vec deriv = {0.0, 0.0, 0.0, 0.0};

  arma::vec deriv_indi(4);
  double logdenom;

  for (int i = 0; i < n; i++) {
    logdenom = -arma::datum::inf;
    deriv_indi.fill(-arma::datum::inf);
    for (int gA = 0; gA <= K; gA++) {
      for (int gB = 0; gB <= K; gB++) {
        deriv_indi = plog_sum_exp(deriv_indi, pgA(i, gA) + pgB(i, gB) + darray.tube(gA, gB));
        logdenom = log_sum_exp_2(logdenom, pgA(i, gA) + pgB(i, gB) + parray(gA, gB));
      }
    }
    deriv += arma::exp(deriv_indi - logdenom);
  }

  return deriv;
}


//' Derivative of \code{\link{llike_genolike}()} with respect to \code{par}.
//'
//' @inheritParams llike_genolike
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::vec dllike_genolike_dpar(const arma::vec par,
                               const arma::mat &pgA,
                               const arma::mat &pgB,
                               const arma::vec alpha) {
  if (par.n_elem != 3) {
    Rcpp::stop("dllike_genolike_dpar: par needs to be length 3");
  }

  arma::mat dp_dy = dreal_to_simplex_dy(par);
  arma::vec prob = real_to_simplex(par);
  arma::vec df_dp = dproballgenolike_dprob(pgA, pgB, prob) +
    dlprior_dprob(prob, alpha);

  arma::vec deriv = {0.0, 0.0, 0.0};
  for (int i = 0; i < 4; i++) {
    deriv[0] += df_dp[i] * dp_dy(i, 0);
    deriv[1] += df_dp[i] * dp_dy(i, 1);
    deriv[2] += df_dp[i] * dp_dy(i, 2);
  }

  return deriv;
}



