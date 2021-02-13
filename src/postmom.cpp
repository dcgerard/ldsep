#include <RcppArmadillo.h>
using namespace Rcpp;

// pm matrix to be filled with posterior means
// gp posterior probability array. See ldfast()
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
// gp posterior probability array. See ldfast()
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

//' Calculate just the standard errors from genotype posterior array.
//'
//' Only pairwise complete observations are used to calculate standard errors.
//'
//' @param gp A three-way array with dimensions SNPs by individuals by dosage.
//'     That is, \code{gp[i, j, k]} is the posterior probability of
//'     dosage \code{k-1} for individual \code{j} at SNP \code{i}.
//' @param pm_mat The matrix of posterior mean genotypes for each individual.
//'     Rows index SNPs and columns index individuals.
//' @param pv_mat The matrix of posterior variances for each individual.
//'     Rows index SNPs and columns index individuals.
//' @param type a = D, b = r, c = D'
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericMatrix secalc(const arma::cube &gp,
                     const NumericMatrix &pm_mat,
                     const NumericMatrix &pv_mat,
                     char type) {
  if ((type != 'a') && (type != 'b') && (type != 'c')) {
    Rcpp::stop("type should be a, b, or c");
  }

  int nsnp = gp.n_rows; // number of SNPs
  int nind = gp.n_cols; //  number of individuals
  int ploidy = gp.n_slices - 1; // ploidy of species
  double pd = (double)ploidy; // double version of ploidy

  arma::vec Mi(7); // moments for each individual between two loci
  arma::mat Omega(7, 7); // Sample covariance between moments using just pairwise complete observations
  arma::vec Mbar(7); // Sample mean using just pairwise complete obs.
  arma::vec grad(7); // gradient for transformation from M to LD measure.
  int n; // sample size for pairwise complete observations;
  double one_over_n; // One over n
  double nm1_over_n; // (n-1)/n

  NumericMatrix semat(nsnp, nsnp);
  std::fill(semat.begin(), semat.end(), NA_REAL);

  for (int i = 0; i < nsnp; i++) {
    for (int j = i; j < nsnp; j++) {
      n = 0;
      Mbar.zeros();
      Omega.zeros();
      for (int ell = 0; ell < nind; ell++) {
        if (!NumericMatrix::is_na(pm_mat(i, ell)) &&
            !NumericMatrix::is_na(pm_mat(j, ell))) {
            n++;
          one_over_n = 1.0 / (double)n;
          nm1_over_n = ((double)n - 1.0) / (double)n;

          Mi(0) = pm_mat(i, ell);
          Mi(1) = std::pow(pm_mat(i, ell), 2.0);
          Mi(2) = pm_mat(j, ell);
          Mi(3) = std::pow(pm_mat(j, ell), 2.0);
          Mi(4) = pm_mat(i, ell) * pm_mat(j, ell);
          Mi(5) = pv_mat(i, ell);
          Mi(6) = pv_mat(j, ell);

          Mbar = Mbar * nm1_over_n + Mi * one_over_n;
          Omega = Omega * nm1_over_n + (Mi * Mi.t()) * one_over_n;
        }
      }
      Omega = (Omega - (Mbar * Mbar.t())) * (double)n / ((double)n - 1.0);

      if ((type == 'a') && (i != j)) {
        // D
        grad_delta_m(Mbar, grad, pd);
        semat(i, j) = std::sqrt((grad.t() * Omega * grad).eval()(0, 0) / (double)n);
      } else if ((type == 'b') && (i != j)) {
        // r
        grad_rho_m(Mbar, grad);
        semat(i, j) = std::sqrt((grad.t() * Omega * grad).eval()(0, 0) / (double)n);
      } else if ((type == 'c') && (i != j)) {
        // D'
        grad_deltaprime_m(Mbar, grad, pd);
        semat(i, j) = std::sqrt((grad.t() * Omega * grad).eval()(0, 0) / (double)n);
      }
      semat(j, i) = semat(i, j);
    }
  }

  return semat;
}
