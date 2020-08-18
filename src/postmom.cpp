#include <RcppArmadillo.h>
using namespace Rcpp;

void ds_from_gp(arma::mat &ds, const arma::cube &gp);
void pv_from_gp(arma::vec &pv,
                const arma::cube &gp,
                const arma::mat &ds,
                const int &i);

//' Fast bias-correction for LD
//'
//' This one does not assume any information about prior moments.
//'
//' @inheritParams ldfast
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
List ldfast_post(const arma::cube &gp,
                      char type,
                      Nullable<Rcpp::NumericVector> priorvar_ = R_NilValue) {
  int nsnp = gp.n_rows;
  int nind = gp.n_cols;
  int ploidy = gp.n_slices - 1;
  double pd = (double)ploidy;

  arma::mat ds(nind, nsnp); // posterior means. Notice reversing of snps/inds for armadillo
  arma::vec pv(nind); // posterior variances
  arma::vec rr(nsnp); // reliability ratio

  arma::vec mean_pv(nsnp); // mean of posterior variances
  arma::vec var_ds(nsnp); // variance of posterior means

  ds_from_gp(ds, gp); // populates ds with posterior means

  if (priorvar_.isNotNull()) {
    Rcpp::NumericVector priorvar(priorvar_);
    for (int i = 0; i < nsnp; i++) {
      pv_from_gp(pv, gp, ds, i);
      mean_pv(i) = arma::mean(pv);
      rr(i) = std::sqrt(1.0 / (1.0 - mean_pv(i) / priorvar(i)));
    }
  } else{
    for (int i = 0; i < nsnp; i++) {
      pv_from_gp(pv, gp, ds, i);
      mean_pv(i) = arma::mean(pv);
      var_ds(i) = arma::var(ds.col(i));
      rr(i) = std::sqrt((mean_pv(i) + var_ds(i)) / var_ds(i));
    }
  }

  arma::mat cmat(nsnp, nsnp);
  cmat = arma::diagmat(rr) * arma::cor(ds) * arma::diagmat(rr);

  // Correct for too much expansion
  for (int i = 0; i < nsnp; i++) {
    for (int j = 0; j < nsnp; j++) {
      if (cmat(i, j) < -1.0) {
        cmat(i, j) = -1.0;
      } else if (cmat(i, j) > 1.0) {
        cmat(i, j) = 1.0;
      }
    }
  }

  if (type == 'D') {
    arma::vec sdvec = arma::sqrt(mean_pv + var_ds);
    cmat = (arma::diagmat(sdvec) * cmat * arma::diagmat(sdvec)) / pd;
  }

  List L = List::create(Named("cor") = cmat , _["rr"] = rr);

  return L;
}

//' Calculate posterior mean from posterior probs
//'
//' @param ds The matrix of posterior means to populate of dimension
//'     ind by snp. This is non-standard dimension order and is for
//'     armadillo to calculate covariance efficiently.
//' @inheritParams ldfast_post
//'
//' @noRd
//'
//' @author David Gerard
// [[Rcpp::export]]
void ds_from_gp(arma::mat &ds, const arma::cube &gp) {
  int nsnp = gp.n_rows;
  int nind = gp.n_cols;
  int ploidy = gp.n_slices - 1;

  ds.fill(0.0);
  for (int i = 0; i < nsnp; i++) {
    for (int j = 0; j < nind; j++) {
      for (int k = 0; k <= ploidy; k++) {
        ds(j, i) += (double)k * gp(i, j, k);
      }
    }
  }
}


//' Calculate posterior variance from posterior probs and posterior mean
//'
//' @param pv The posterior variance vector to be filled
//' @param ds The matrix of posterior means. TAKE NOTE:
//'     Columns index SNPs and rows index individuals.
//' @param i The current snp to look at.
//' @inheritParams ldfast_post
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
void pv_from_gp(arma::vec &pv,
                const arma::cube &gp,
                const arma::mat &ds,
                const int &i) {
  int nind = gp.n_cols;
  int ploidy = gp.n_slices - 1;
  pv.fill(0.0);
  for (int j = 0; j < nind; j++) {
    for (int k = 0; k <= ploidy; k++) {
      pv(j) += std::pow((double)k, 2.0) * gp(i, j, k);
    }
    pv(j) -= std::pow(ds(j, i), 2.0);
  }
}
