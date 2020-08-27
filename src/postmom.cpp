#include <RcppArmadillo.h>
using namespace Rcpp;

// pm matrix to be filled with posterior means
// gp posterior probability array. See ldfast_calc
// [[Rcpp::export]]
void fill_pm(NumericMatrix &pm, const arma::cube &gp) {
  int nsnp = gp.n_rows;
  int nind = gp.n_cols;
  int ploidy = gp.n_slices - 1;

  pm.fill(0.0);
  for (int i = 0; i < nsnp; i++) {
    for (int j = 0; j < nind; j++) {
      for (int k = 0; k <= ploidy; k++) {
        pm(i, j) += (double)k * gp(i, j, k);
      }
    }
  }
}

// pv matrix to be filled with posterior variances.
// pm matrix of posterior means.
// gp posterior probability array. See ldfast_calc.
// [[Rcpp::export]]
void fill_pv(NumericMatrix &pv, NumericMatrix &pm, const arma::cube &gp) {
  int nsnp = gp.n_rows;
  int nind = gp.n_cols;
  int ploidy = gp.n_slices - 1;

  pv.fill(0.0);
  for (int i = 0; i < nsnp; i++) {
    for (int j = 0; j < nind; j++) {
      for (int k = 0; k <= ploidy; k++) {
        pv(i, j) += std::pow((double)k - pm(i, j), 2.0) * gp(i, j, k);
      }
    }
  }
}

// M the mean moment vector
// grad the gradient vector to be filled
// pd the ploidy
// [[Rcpp::export]]
void grad_delta_m(const arma::vec &M, arma::vec &grad, double pd) {
  grad(0) = -1.0 * (std::pow(M(0), 4.0) * M(2) - 2.0 * M(0) * M(4) * M(5) + std::pow(M(0), 2.0) * M(2) * (M(5) - 2.0 * M(1)) + M(1) * M(2) * (M(1) + M(5))) * (std::pow(M(2), 2.0) - M(3) - M(6)) / (pd * std::pow(std::pow(M(0), 2.0) - M(1), 2.0) * (std::pow(M(2), 2.0) - M(3)));
  grad(1) = (M(0) * M(2) - M(4)) * M(5) * (std::pow(M(2), 2.0) - M(3) - M(6)) / (pd * std::pow(std::pow(M(0), 2.0) - M(1), 2.0) * (std::pow(M(2), 2.0) - M(3)));
  grad(2) = -1.0 * (std::pow(M(0), 2.0) - M(1) - M(5)) * (-2.0 * M(2) * M(4) * M(6) + M(0) * (std::pow(M(2), 4.0) + std::pow(M(2), 2.0) * (-2.0 * M(3) + M(6)) + M(3) * (M(3) + M(6)))) / (pd * (std::pow(M(0), 2.0) - M(1)) * std::pow(std::pow(M(2), 2.0) - M(3), 2.0));
  grad(3) = (M(0) * M(2) - M(4)) * (std::pow(M(0), 2.0) - M(1) - M(5)) * M(6) / (pd * (std::pow(M(0), 2.0) - M(1)) * std::pow(std::pow(M(2), 2.0) - M(3), 2.0));
  grad(4) = (-1.0 * std::pow(M(0), 2.0) + M(1) + M(5)) * (-1.0 * std::pow(M(2), 2.0) + M(3) + M(6)) / (pd * (std::pow(M(0), 2.0) - M(1)) * (std::pow(M(2), 2.0) - M(3)));
  grad(5) = (-1.0 * M(0) * M(2) + M(4)) * (-1.0 * std::pow(M(2), 2.0) + M(3) + M(6)) / (pd * (std::pow(M(0), 2.0) - M(1)) * (std::pow(M(2), 2.0) - M(3)));
  grad(6) = (-1.0 * M(0) * M(2) + M(4)) * (-1.0 * std::pow(M(0), 2.0) + M(1) + M(5)) / (pd * (std::pow(M(0), 2.0) - M(1)) * (std::pow(M(2), 2.0) - M(3)));
}

// M the mean moment vector
// grad the gradient vector to be filled
// pd the ploidy
// [[Rcpp::export]]
void grad_deltaprime_m(const arma::vec &M, arma::vec &grad, double pd) {
  grad_delta_m(M, grad, pd);
  double c1;
  double c3;

  double Delta = ((M(5) + M(1) - std::pow(M(0), 2.0)) / (M(1) - std::pow(M(0), 2.0))) * ((M(6) + M(3) - std::pow(M(2), 2.0)) / (M(3) - std::pow(M(2), 2.0))) * ((M(4) - M(0) * M(2)) / pd);
  double Delta_m;
  if (M(4) < M(0) * M(2)) {
    Delta_m = std::min(M(0) * M(2), (pd - M(0)) * (pd - M(2))) / std::pow(pd, 2.0);
  } else {
    Delta_m = std::min(M(0) * (pd - M(2)), (pd - M(0)) * M(2)) / std::pow(pd, 2.0);
  }

  if ((M(4) < M(0) * M(2)) && (M(0) * M(2) < (pd - M(0)) * (pd - M(2)))) {
    c1 = M(2);
    c3 = M(0);
  } else if ((M(4) < M(0) * M(2)) && (M(0) * M(2) > (pd - M(0)) * (pd - M(2)))) {
    c1 = -1.0 * (pd - M(2));
    c3 = -1.0 * (pd - M(0));
  } else if ((M(4) > M(0) * M(2)) && (M(0) * (pd - M(2)) > (pd - M(0)) * M(2))) {
    c1 = -1.0 * M(2);
    c3 = pd - M(0);
  } else {
    c1 = pd - M(2);
    c3 = -1.0 * M(0);
  }
  c1 = c1 / std::pow(pd, 2.0);
  c3 = c3 / std::pow(pd, 2.0);

  grad = grad/ Delta_m;

  grad(0) = grad(0) - Delta * c1 / std::pow(Delta_m, 2.0);
  grad(2) = grad(2) - Delta * c3 / std::pow(Delta_m, 2.0);
}

// M the mean moment vector
// grad the gradient vector to be filled
// [[Rcpp::export]]
void grad_rho_m(const arma::vec &M, arma::vec &grad) {
  grad(0) = (std::pow(M(0), 3.0) * M(4) + std::pow(M(0), 2.0) * M(2) * (M(5) - M(1)) + M(1) * M(2) * (M(1) + M(5)) - M(0) * M(4) * (M(1) + 2.0 * M(5))) * std::sqrt(M(3) + M(6) - std::pow(M(2), 2.0)) / (std::pow(std::pow(M(0), 2.0) - M(1), 2.0) * (std::pow(M(2), 2.0) - M(3)) * std::sqrt(M(1) + M(5) - std::pow(M(0), 2.0)));
  grad(1) = (M(0) * M(2) - M(4)) * (std::pow(M(0), 2.0) - M(1) - 2.0 * M(5)) * std::sqrt(M(3) + M(6) - std::pow(M(2), 2.0)) / (2.0 * std::pow(std::pow(M(0), 2.0) - M(1), 2.0) * (std::pow(M(2), 2.0) - M(3)) * std::sqrt(M(1) + M(5) - std::pow(M(0), 2.0)));
  grad(2) = -1.0 * std::sqrt(M(1) + M(5) - std::pow(M(0), 2.0)) * (M(0) * std::pow(M(2), 2.0) * (M(3) - M(6)) - M(0) * M(3) * (M(3) + M(6)) + M(2) * M(4) * (-1.0 * std::pow(M(2), 2.0) + M(3) + 2.0 * M(6))) / ((std::pow(M(0), 2.0) - M(1)) * std::pow(std::pow(M(2), 2.0) - M(3), 2.0) * std::sqrt(M(3) + M(6) - std::pow(M(2), 2.0)));
  grad(3) = (M(0) * M(2) - M(4)) * std::sqrt(M(1) + M(5) - std::pow(M(0), 2.0)) * (std::pow(M(2), 2.0) - M(3) - 2.0 * M(6)) / (2.0 * (std::pow(M(0), 2.0) - M(1)) * std::pow(std::pow(M(2), 2.0) - M(3), 2.0) * std::sqrt(M(3) + M(6) - std::pow(M(2), 2.0)));
  grad(4) = std::sqrt(M(1) + M(5) - std::pow(M(0), 2.0)) * std::sqrt(M(3) + M(6) - std::pow(M(2), 2.0)) / ((std::pow(M(0), 2.0) - M(1)) * (std::pow(M(2), 2.0) - M(3)));
  grad(5) = (M(4) - M(0) * M(2)) * std::sqrt(M(3) + M(6) - std::pow(M(2), 2.0)) / (2.0 * (std::pow(M(0), 2.0) - M(1)) * (std::pow(M(2), 2.0) - M(3)) * std::sqrt(M(1) + M(5) - std::pow(M(0), 2.0)));
  grad(6) = (M(4) - M(0) * M(2)) * std::sqrt(M(1) + M(5) - std::pow(M(0), 2.0)) / (2.0 * (std::pow(M(0), 2.0) - M(1)) * (std::pow(M(2), 2.0) - M(3)) * std::sqrt(M(3) + M(6) - std::pow(M(2), 2.0)));
}


//' Fast bias-correction for LD
//'
//' This one does not assume any information about prior moments.
//' This modifies two empty matrices: the correlation and the standard error
//' matrices.
//'
//' @param cormat The matrix that will hold the correlations.
//' @param semat The matrix that will hold the standard errors.
//' @param gp A three-way array with dimensions SNPs by individuals by dosage.
//'     That is, \code{gp[i, j, k]} is the posterior probability of
//'     dosage \code{k-1} for individual \code{j} at SNP \code{i}.
//' @param type a = D, b = r, c = D'
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
void ldfast_calc(NumericMatrix &cormat,
                 NumericMatrix &semat,
                 const arma::cube &gp,
                 char type) {
  int nsnp = gp.n_rows; // number of SNPs
  int nind = gp.n_cols; //  number of individuals
  int ploidy = gp.n_slices - 1; // ploidy of species
  double pd = (double)ploidy; // double version of ploidy

  if (cormat.nrow() != cormat.ncol()) {
    Rcpp::stop("cormat should have same number of rows as columns");
  }
  if (cormat.nrow() != nsnp) {
    Rcpp::stop("cormat should have dimensions of nsnp");
  }
  if (semat.nrow() != semat.ncol()) {
    Rcpp::stop("semat should have same number of rows as columns");
  }
  if (semat.nrow() != nsnp) {
    Rcpp::stop("semat should have dimensions of nsnp");
  }

  // posterior means and variances
  NumericMatrix pm_mat(nsnp, nind); // posterior mean matrix
  NumericMatrix pv_mat(nsnp, nind); // posterior variance matrix
  fill_pm(pm_mat, gp);
  fill_pv(pv_mat, pm_mat, gp);

  // M_i = (X_{iA}, X_{iA}^2, X_{iB}, X_{iB}^2, X_{iA}X_{iB}, Y_{iA}, Y_{iB})
  arma::vec Mi(7); // moments for each individual between two loci
  arma::mat Omega(7, 7); // Sample covariance between moments
  arma::vec Mbar(7); // Sample mean of moments.
  arma::vec grad(7); // gradient for transformation from M to LD measure.
  int n; // sample size for each snp combo;

  double uxa; // mean of posterior means at locus 1.
  double uxb; // mean of posterior means at locus 2.
  double vxa; // variance of posterior means at locus 1.
  double vxb; // variance of posterior means at locus 2.
  double cx; // covariance of posterior means at loci 1 and 2.
  double uya; // mean of posterior variances at locus 1.
  double uyb; // mean of posterior variances at locus 2.

  double Deltam; // bound of delta conditional on means

  // Fill in correlations
  for (int i = 0; i < nsnp; i++) {
    for (int j = i; j < nsnp; j++) {
      n = 0;
      Mbar.zeros();
      Omega.zeros();
      for (int ell = 0; ell < nind; ell++) {
        if (!NumericMatrix::is_na(pm_mat(i, ell)) &&
            !NumericMatrix::is_na(pm_mat(j, ell))) {
          n++;
          Mi(0) = pm_mat(i, ell);
          Mi(1) = std::pow(pm_mat(i, ell), 2.0);
          Mi(2) = pm_mat(j, ell);
          Mi(3) = std::pow(pm_mat(j, ell), 2.0);
          Mi(4) = pm_mat(i, ell) * pm_mat(j, ell);
          Mi(5) = pv_mat(i, ell);
          Mi(6) = pv_mat(j, ell);
          Mbar = Mbar + Mi;
          Omega = Omega + (Mi * Mi.t());
        }
      }
      Mbar = Mbar / (double)n;
      Omega = Omega / ((double)n - 1.0) -
        ((double)n) / ((double)n - 1.0) * (Mbar * Mbar.t());

      // Rcpp::Rcout << Omega << std::endl << std::endl;

      uxa = Mbar(0);
      uxb = Mbar(2);
      vxa = (Mbar(1) - std::pow(Mbar(0), 2.0)) * ((double)n) / ((double)n - 1.0);
      vxb = (Mbar(3) - std::pow(Mbar(2), 2.0)) * ((double)n) / ((double)n - 1.0);
      cx =  (Mbar(4) - Mbar(0) * Mbar(2)) * ((double)n) / ((double)n - 1.0);
      uya = Mbar(5);
      uyb = Mbar(6);

      if (type == 'a') {
        // D
        if (i == j) {
          cormat(i, j) = (uya + vxa) / pd;
        } else {
          cormat(i, j) = ((uya + vxa) / vxa) * ((uyb + vxb) / vxb) * (cx / pd);

          if (cormat(i, j) > std::sqrt((uya + vxa) * (uyb + vxb)) / pd) {
            cormat(i, j) = std::sqrt((uya + vxa) * (uyb + vxb)) / pd;
          } else if (cormat(i, j) < -1.0 * std::sqrt((uya + vxa) * (uyb + vxb)) / pd) {
            cormat(i, j) = -1.0 * std::sqrt((uya + vxa) * (uyb + vxb)) / pd;
          }

          // populate grad and calculate SE
          grad_delta_m(Mbar, grad, pd);
          semat(i, j) = std::sqrt((grad.t() * Omega * grad).eval()(0, 0) / (double)n);
        }
      } else if (type == 'b') {
        // r
        if (i == j) {
          cormat(i, j) = 1.0;
        } else {
          cormat(i, j) = std::sqrt((uya + vxa) / vxa) * std::sqrt((uyb + vxb) / vxb) * (cx / std::sqrt(vxa * vxb));
          if (cormat(i, j) > 1.0) {
            cormat(i, j) = 1.0;
          } else if (cormat(i, j) < -1.0) {
            cormat(i, j) = -1.0;
          }

          // populate grad and calculate SE
          grad_rho_m(Mbar, grad);
          semat(i, j) = std::sqrt((grad.t() * Omega * grad).eval()(0, 0) / (double)n);
        }
      } else if (type == 'c') {
        // D'
        if (i == j) {
          cormat(i, j) = pd;
        } else {
          if (cx < 0.0) {
            Deltam = std::min(uxa * uxb, (pd - uxa) * (pd - uxb)) / std::pow(pd, 2.0);
          } else {
            Deltam = std::min(uxa * (pd - uxb), (pd - uxa) * uxb) / std::pow(pd, 2.0);
          }

          cormat(i, j) = ((uya + vxa) / vxa) * ((uyb + vxb) / vxb) * (cx / pd);
          cormat(i, j) = cormat(i, j) / Deltam;

          if (cormat(i, j) > pd) {
            cormat(i, j) = pd;
          } else if (cormat(i, j) < -1.0 * pd) {
            cormat(i, j) = -1.0 * pd;
          }

          // populate grad and calculate SE
          grad_deltaprime_m(Mbar, grad, pd);
          semat(i, j) = std::sqrt((grad.t() * Omega * grad).eval()(0, 0) / (double)n);
        }
      } else {
        Rcpp::stop("type should be a, b, or c");
      }
    }
  }
}

// OLD CODE FOR REMOVAL LATER -------------------------------------------------

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


//' Calculate reliability ratio using just posterior moments
//'
//' @param gp The three-way array of posterior probabilities of dimension
//'     SNPs by individuals by dosage.
//' @param ds The matrix of posterior means of dimension individuals by SNPs
//'
//' @return A matrix of two columns. The first column contains the
//'     estimated reliability ratios for the correlations, the second
//'     column contains the estimated prior standard deviations.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::mat post_rr(const arma::cube &gp, const arma::mat &ds) {
  int nsnp = gp.n_rows;
  int nind = gp.n_cols;

  double mean_pv; // mean of posterior variances
  double var_ds; // variance of posterior means
  arma::vec pv(nind); // posterior variances
  arma::mat rr(nsnp, 2); // reliability ratio
  for (int i = 0; i < nsnp; i++) {
    pv_from_gp(pv, gp, ds, i);
    mean_pv = arma::mean(pv);
    var_ds = arma::var(ds.col(i));
    rr(i, 0) = std::sqrt((mean_pv + var_ds) / var_ds);
    rr(i, 1) = std::sqrt(mean_pv + var_ds);
  }
  return rr;
}


//' Calculate reliability ratio using just posterior moments
//'
//' @param gp The three-way array of posterior probabilities of dimension
//'     SNPs by individuals by dosage.
//' @param ds The matrix of posterior means of dimension individuals by SNPs
//' @param priorvar The vector of prior variances
//'
//' @return A matrix of two columns. The first column contains the
//'     estimated reliability ratios for the correlations, the second
//'     column contains the estimated prior standard deviations.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
arma::mat prior_rr(const arma::cube &gp, const arma::mat &ds, const arma::vec &priorvar) {
  int nsnp = gp.n_rows;
  int nind = gp.n_cols;

  double mean_pv; // mean of posterior variances
  arma::vec pv(nind); // posterior variances
  arma::mat rr(nsnp, 2); // reliability ratio
  for (int i = 0; i < nsnp; i++) {
    pv_from_gp(pv, gp, ds, i);
    mean_pv = arma::mean(pv);
    rr(i, 0) = std::sqrt(1.0 / (1.0 - mean_pv / priorvar(i)));
    rr(i, 1) = std::sqrt(priorvar(i));
  }
  return rr;
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
