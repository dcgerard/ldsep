context("Environment API")


test_that("llike_geno is unchanged", {
  gA <- rep(0:4, each = 5)
  gB <- rep(0:4, 5)
  K <- 4
  par <- c(-1, 2, 3)

  l1 <- llike_geno(par = par, gA = gA, gB = gB, K = K)
  d1 <- dllike_geno_dpar(par = par, gA = gA, gB = gB, K = K)

  env <- new.env()
  env[["gA"]] <- gA
  env[["gB"]] <- gB
  env[["K"]] <- K
  l2 <- -1 * nllike_geno_p(xs = par, env = env)
  d2 <- -1 * dnllike_geno_dpar_p(xs = par, env = env)

  expect_equal(l1, l2)
  expect_equal(d1, d2)

})


# likelihood.include <-'Rcpp::NumericVector lhood(SEXP xs, SEXP env){arma::vec par = Rcpp::as<arma::vec>(xs);Rcpp::Environment e = Rcpp::as<Rcpp::Environment>(env);arma::mat X = Rcpp::as<arma::mat>(e["X"]);arma::vec y = Rcpp::as<arma::vec>(e["y"]);double prec = Rcpp::as<double>(e["prec"]);arma::mat Xbeta = X * par;double sum1 = sum(y % Xbeta - log(1 + exp(Xbeta)));arma::mat sum2 = sum(pow(par, 2 * prec));arma::vec out = -(sum1 - 0.5 * sum2);Rcpp::NumericVector ret = Rcpp::as<Rcpp::NumericVector>(wrap(out));return ret;}'
#
# likelihood.body <-'typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);return(XPtr<funcPtr>(new funcPtr(&lhood)));'
#
# likelihood.CPP <- cxxfunction(signature(), body=likelihood.body,inc=likelihood.include, plugin="RcppArmadillo")
#






