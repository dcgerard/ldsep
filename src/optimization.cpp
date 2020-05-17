#include <RcppArmadillo.h>
#include <roptim.h>
using namespace Rcpp;

double llike_geno(const NumericVector &par,
                  const IntegerVector &gA,
                  const IntegerVector &gB,
                  const int &K);

NumericVector dllike_geno_dpar(const NumericVector &par,
                               const IntegerVector &gA,
                               const IntegerVector &gB,
                               const int &K);


struct genoCor : public roptim::Functor {
  public:
    arma::vec gA;
    arma::vec gB;
    int K;
    double operator()(const arma::vec &par) override {
      return (double)K;
    }; // objective function
    void Gradient(const arma::vec &par, arma::vec &grad) override; // gradient
    };
