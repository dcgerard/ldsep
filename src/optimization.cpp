#include <RcppArmadillo.h>
#include <roptim.h>
using namespace Rcpp;

double llike_geno(const arma::vec par,
                  const arma::vec &gA,
                  const arma::vec &gB,
                  const int K);

arma::vec dllike_geno_dpar(const arma::vec par,
                           const arma::vec &gA,
                           const arma::vec &gB,
                           const int K);

double llike_genolike(const arma::vec par,
                      const arma::mat &pgA,
                      const arma::mat &pgB);

arma::vec dllike_genolike_dpar(const arma::vec par,
                               const arma::mat &pgA,
                               const arma::mat &pgB);


struct genoCor : public roptim::Functor {
  public:
    arma::vec gA;
    arma::vec gB;
    int K;
    double operator()(const arma::vec &par) override {
      return llike_geno(par, gA, gB, K);
    }; // objective function
    void Gradient(const arma::vec &par, arma::vec &grad) override {
      grad = dllike_geno_dpar(par, gA, gB, K);
    }; // gradient
};

struct genoLikeCor : public roptim::Functor {
public:
  arma::mat pgA;
  arma::mat pgB;
  double operator()(const arma::vec &par) override {
    return llike_genolike(par, pgA, pgB);
  }; // objective function
  void Gradient(const arma::vec &par, arma::vec &grad) override {
    grad = dllike_genolike_dpar(par, pgA, pgB);
  }; // gradient
};

//' Find LD estimates using genotype likelihoods
//'
//'
//' @param par The parameters on the real-scale.
//' @param pgA The genotype log-liielihoods at locus 1.
//' @param pgB The genotype log-likelihoods at locus 2.
//' @param reltol The stopping criterion for the gradient ascent.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
List optimize_genolikecor(arma::vec &par,
                          const arma::mat &pgA,
                          const arma::mat &pgB) {
  genoLikeCor gc;
  gc.pgA = pgA;
  gc.pgB = pgB;
  roptim::Roptim<genoLikeCor> opt("BFGS");
  opt.control.fnscale = -1.0;
  opt.set_hessian(true);
  opt.minimize(gc, par);


  List L = List::create(Named("par") = opt.par() ,
                        _["value"] = opt.value(),
                        _["convergence"] = opt.convergence(),
                        _["message"] = opt.message(),
                        _["hessian"] = -1.0 * opt.hessian(),
                        _["grcount"] = opt.grcount(),
                        _["fncount"] = opt.fncount());

  return L;
}

//' Find LD estimates using just the genotypes.
//'
//'
//' @param par The parameters on the real-scale.
//' @param gA The genotypes at locus 1.
//' @param gB The genotypes at locus 2.
//' @param K The ploidy for the species. Assumed to be the same for
//'     all individuals.
//' @param reltol The stopping criterion for the gradient ascent.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
List optimize_genocor(arma::vec &par,
                      const arma::vec &gA,
                      const arma::vec &gB,
                      const int &K,
                      double reltol = 10.0e-08) {
  genoCor gc;
  gc.gA = gA;
  gc.gB = gB;
  gc.K = K;
  roptim::Roptim<genoCor> opt("BFGS");
  opt.control.fnscale = -1.0;
  opt.control.reltol = reltol;
  opt.set_hessian(true);
  opt.minimize(gc, par);

  List L = List::create(Named("par") = opt.par() ,
                        _["value"] = opt.value(),
                        _["convergence"] = opt.convergence(),
                        _["message"] = opt.message(),
                        _["hessian"] = -1.0 * opt.hessian());

  return L;
}

