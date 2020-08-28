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
  int n; // sample size for pairwise complete observations;
  int na; // sample size for SNP A
  int nb; // sample size for SNP B

  // I need this because I didn't have locus 1 (or 2) contiguous in original code.
  arma::Col<int> aind = {0, 1, 5}; // indices of A locus in Mi
  arma::Col<int> bind = {2, 3, 6}; // indices of B locus in Mi
  int crossind = 4; // index of cross product

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
      na = 0;
      nb = 0;
      Mbar.zeros();
      Omega.zeros();

      for (int ell = 0; ell < nind; ell++) {
        if (!NumericMatrix::is_na(pm_mat(i, ell))) {
          // update locus A moments
          na++;
          Mi(0) = pm_mat(i, ell);
          Mi(1) = std::pow(pm_mat(i, ell), 2.0);
          Mi(5) = pv_mat(i, ell);

          Mbar(0) = Mbar(0) * ((double)na - 1.0) / (double)na + Mi(0) / (double)na;
          Mbar(1) = Mbar(1) * ((double)na - 1.0) / (double)na + Mi(1) / (double)na;
          Mbar(5) = Mbar(5) * ((double)na - 1.0) / (double)na + Mi(5) / (double)na;

          for (int ai = 0; ai < 3; ai++) {
            for (int aj = 0; aj < 3; aj++) {
              Omega(aind(ai), aind(aj)) = Omega(aind(ai), aind(aj)) * ((double)na - 1.0) / (double)na +
                Mi(aind(ai)) * Mi(aind(aj)) / (double)na;
            }
          }

        }

        if (!NumericMatrix::is_na(pm_mat(j, ell))) {
          // update locus B moments
          nb++;
          Mi(2) = pm_mat(j, ell);
          Mi(3) = std::pow(pm_mat(j, ell), 2.0);
          Mi(6) = pv_mat(j, ell);

          Mbar(2) = Mbar(2) * ((double)nb - 1.0) / (double)nb + Mi(2) / (double)nb;
          Mbar(3) = Mbar(3) * ((double)nb - 1.0) / (double)nb + Mi(3) / (double)nb;
          Mbar(6) = Mbar(6) * ((double)nb - 1.0) / (double)nb + Mi(6) / (double)nb;

          for (int bi = 0; bi < 3; bi++) {
            for (int bj = 0; bj < 3; bj++) {
              Omega(bind(bi), bind(bj)) = Omega(bind(bi), bind(bj)) * ((double)nb - 1.0) / (double)nb +
                Mi(bind(bi)) * Mi(bind(bj)) / (double)nb;
            }
          }

        }

        if (!NumericMatrix::is_na(pm_mat(i, ell)) &&
            !NumericMatrix::is_na(pm_mat(j, ell))) {
          // update cross product moments
          n++;
          Mi(4) = pm_mat(i, ell) * pm_mat(j, ell);
          Mbar(4) = Mbar(4) * ((double)n - 1.0) / (double)n + Mi(4) / (double)n;

          Omega(crossind, crossind) = Omega(crossind, crossind) * ((double)n - 1.0) / (double)n +
            std::pow(Mi(4), 2.0) / (double)n;

          for (int ai = 0; ai < 3; ai++) {
            Omega(aind(ai), crossind) = Omega(aind(ai), crossind) * ((double)n - 1.0) / (double)n +
              Mi(aind(ai)) * Mi(crossind) / (double)n;
            Omega(crossind, aind(ai)) = Omega(aind(ai), crossind);
          }

          for (int bi = 0; bi < 3; bi++) {
            Omega(bind(bi), crossind) = Omega(bind(bi), crossind) * ((double)n - 1.0) / (double)n +
              Mi(bind(bi)) * Mi(crossind) / (double)n;
            Omega(crossind, bind(bi)) = Omega(bind(bi), crossind);
          }

          for (int ai = 0; ai < 3; ai++) {
            for (int bi = 0; bi < 3; bi++) {
              Omega(aind(ai), bind(bi)) = Omega(aind(ai), bind(bi)) * ((double)n - 1.0) / (double)n +
                Mi(aind(ai)) * Mi(bind(bi)) / (double)n;
              Omega(bind(bi), aind(ai)) = Omega(aind(ai), bind(bi));
            }
          }
        }
      }
      Omega = Omega - Mbar * Mbar.t();

      // Update cov denominator to make unbiased
      Omega(crossind, crossind) = Omega(crossind, crossind) * (double)n / ((double)n - 1.0);
      for (int ai = 0; ai < 3; ai++) {
        for (int aj = 0; aj < 3; aj++) {
          Omega(aind(ai), aind(aj)) = Omega(aind(ai), aind(aj)) * (double)na / ((double)na - 1.0);
        }
      }
      for (int bi = 0; bi < 3; bi++) {
        for (int bj = 0; bj < 3; bj++) {
          Omega(bind(bi), bind(bj)) = Omega(bind(bi), bind(bj)) * (double)nb / ((double)nb - 1.0);
        }
      }
      for (int ai = 0; ai < 3; ai++) {
        for (int bi = 0; bi < 3; bi++) {
          Omega(aind(ai), bind(bi)) = Omega(aind(ai), bind(bi)) * (double)n / ((double)n - 1.0);
          Omega(bind(bi), aind(ai)) = Omega(aind(ai), bind(bi));
        }
      }
      for (int ab = 0; ab < 3; ab++) {
        Omega(aind(ab), crossind) = Omega(aind(ab), crossind) * (double)n / ((double)n - 1.0);
        Omega(crossind, aind(ab)) = Omega(aind(ab), crossind);

        Omega(bind(ab), crossind) = Omega(bind(ab), crossind) * (double)n / ((double)n - 1.0);
        Omega(crossind, bind(ab)) = Omega(bind(ab), crossind);
      }

      // Calculate central posterior moments
      uxa = Mbar(0);
      uxb = Mbar(2);
      vxa = (Mbar(1) - std::pow(Mbar(0), 2.0)) * ((double)na) / ((double)na - 1.0);
      vxb = (Mbar(3) - std::pow(Mbar(2), 2.0)) * ((double)nb) / ((double)nb - 1.0);
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
            semat(i, j) = NA_REAL;
          } else if (cormat(i, j) < -1.0 * std::sqrt((uya + vxa) * (uyb + vxb)) / pd) {
            cormat(i, j) = -1.0 * std::sqrt((uya + vxa) * (uyb + vxb)) / pd;
            semat(i, j) = NA_REAL;
          } else {
            // populate grad and calculate SE
            grad_delta_m(Mbar, grad, pd);
            semat(i, j) = std::sqrt((grad.t() * Omega * grad).eval()(0, 0) / (double)n);
          }
        }
      } else if (type == 'b') {
        // r
        if (i == j) {
          cormat(i, j) = 1.0;
        } else {
          cormat(i, j) = std::sqrt((uya + vxa) / vxa) * std::sqrt((uyb + vxb) / vxb) * (cx / std::sqrt(vxa * vxb));

          if (cormat(i, j) > 1.0) {
            cormat(i, j) = 1.0;
            semat(i, j) = NA_REAL;
          } else if (cormat(i, j) < -1.0) {
            cormat(i, j) = -1.0;
            semat(i, j) = NA_REAL;
          } else {
            // populate grad and calculate SE
            grad_rho_m(Mbar, grad);
            semat(i, j) = std::sqrt((grad.t() * Omega * grad).eval()(0, 0) / (double)n);
          }
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
            semat(i, j) = NA_REAL;
          } else if (cormat(i, j) < -1.0 * pd) {
            cormat(i, j) = -1.0 * pd;
            semat(i, j) = NA_REAL;
          } else {
            // populate grad and calculate SE
            grad_deltaprime_m(Mbar, grad, pd);
            semat(i, j) = std::sqrt((grad.t() * Omega * grad).eval()(0, 0) / (double)n);
          }

        }
      } else {
        Rcpp::stop("type should be a, b, or c");
      }
    }
  }
}
