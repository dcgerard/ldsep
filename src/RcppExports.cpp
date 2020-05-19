// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dmulti_dprob
arma::vec dmulti_dprob(const arma::vec x, const arma::vec prob, bool log_p);
RcppExport SEXP _ldsep_dmulti_dprob(SEXP xSEXP, SEXP probSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    Rcpp::traits::input_parameter< bool >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(dmulti_dprob(x, prob, log_p));
    return rcpp_result_gen;
END_RCPP
}
// dprobgeno_dprob
arma::vec dprobgeno_dprob(const int& gA, const int& gB, const int K, const arma::vec prob);
RcppExport SEXP _ldsep_dprobgeno_dprob(SEXP gASEXP, SEXP gBSEXP, SEXP KSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type gA(gASEXP);
    Rcpp::traits::input_parameter< const int& >::type gB(gBSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(dprobgeno_dprob(gA, gB, K, prob));
    return rcpp_result_gen;
END_RCPP
}
// dproballgeno_dprob
arma::vec dproballgeno_dprob(const arma::vec& gA, const arma::vec& gB, const int K, const arma::vec prob);
RcppExport SEXP _ldsep_dproballgeno_dprob(SEXP gASEXP, SEXP gBSEXP, SEXP KSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type gA(gASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gB(gBSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(dproballgeno_dprob(gA, gB, K, prob));
    return rcpp_result_gen;
END_RCPP
}
// dreal_to_simplex_dy
arma::mat dreal_to_simplex_dy(const arma::vec y);
RcppExport SEXP _ldsep_dreal_to_simplex_dy(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(dreal_to_simplex_dy(y));
    return rcpp_result_gen;
END_RCPP
}
// dsimplex_to_real_dx
arma::mat dsimplex_to_real_dx(const arma::vec x);
RcppExport SEXP _ldsep_dsimplex_to_real_dx(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(dsimplex_to_real_dx(x));
    return rcpp_result_gen;
END_RCPP
}
// dllike_geno_dpar
arma::vec dllike_geno_dpar(const arma::vec par, const arma::vec& gA, const arma::vec& gB, const int K);
RcppExport SEXP _ldsep_dllike_geno_dpar(SEXP parSEXP, SEXP gASEXP, SEXP gBSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gA(gASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gB(gBSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(dllike_geno_dpar(par, gA, gB, K));
    return rcpp_result_gen;
END_RCPP
}
// dD_dprob
arma::vec dD_dprob(const arma::vec prob);
RcppExport SEXP _ldsep_dD_dprob(SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(dD_dprob(prob));
    return rcpp_result_gen;
END_RCPP
}
// dr2_dprob
arma::vec dr2_dprob(const arma::vec prob);
RcppExport SEXP _ldsep_dr2_dprob(SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(dr2_dprob(prob));
    return rcpp_result_gen;
END_RCPP
}
// dDprime_dprob
arma::vec dDprime_dprob(const arma::vec prob);
RcppExport SEXP _ldsep_dDprime_dprob(SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(dDprime_dprob(prob));
    return rcpp_result_gen;
END_RCPP
}
// probgenolike
double probgenolike(const arma::vec& pgA, const arma::vec& pgB, const arma::vec prob, bool log_p);
RcppExport SEXP _ldsep_probgenolike(SEXP pgASEXP, SEXP pgBSEXP, SEXP probSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pgA(pgASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pgB(pgBSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    Rcpp::traits::input_parameter< bool >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(probgenolike(pgA, pgB, prob, log_p));
    return rcpp_result_gen;
END_RCPP
}
// proballgenolike
double proballgenolike(const arma::mat& pgA, const arma::mat& pgB, const arma::vec prob, bool log_p);
RcppExport SEXP _ldsep_proballgenolike(SEXP pgASEXP, SEXP pgBSEXP, SEXP probSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type pgA(pgASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pgB(pgBSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    Rcpp::traits::input_parameter< bool >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(proballgenolike(pgA, pgB, prob, log_p));
    return rcpp_result_gen;
END_RCPP
}
// llike_genolike
double llike_genolike(const arma::vec par, const arma::mat& pgA, const arma::mat& pgB);
RcppExport SEXP _ldsep_llike_genolike(SEXP parSEXP, SEXP pgASEXP, SEXP pgBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pgA(pgASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pgB(pgBSEXP);
    rcpp_result_gen = Rcpp::wrap(llike_genolike(par, pgA, pgB));
    return rcpp_result_gen;
END_RCPP
}
// get_dprobgeno_dprob_array
arma::cube get_dprobgeno_dprob_array(int K, arma::vec prob);
RcppExport SEXP _ldsep_get_dprobgeno_dprob_array(SEXP KSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(get_dprobgeno_dprob_array(K, prob));
    return rcpp_result_gen;
END_RCPP
}
// get_prob_array
arma::mat get_prob_array(int K, arma::vec prob);
RcppExport SEXP _ldsep_get_prob_array(SEXP KSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(get_prob_array(K, prob));
    return rcpp_result_gen;
END_RCPP
}
// dprobgenolike_dprob
arma::vec dprobgenolike_dprob(const arma::vec& pgA, const arma::vec& pgB, const arma::vec prob);
RcppExport SEXP _ldsep_dprobgenolike_dprob(SEXP pgASEXP, SEXP pgBSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pgA(pgASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pgB(pgBSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(dprobgenolike_dprob(pgA, pgB, prob));
    return rcpp_result_gen;
END_RCPP
}
// dproballgenolike_dprob
arma::vec dproballgenolike_dprob(const arma::mat& pgA, const arma::mat& pgB, const arma::vec prob);
RcppExport SEXP _ldsep_dproballgenolike_dprob(SEXP pgASEXP, SEXP pgBSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type pgA(pgASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pgB(pgBSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(dproballgenolike_dprob(pgA, pgB, prob));
    return rcpp_result_gen;
END_RCPP
}
// dllike_genolike_dpar
arma::vec dllike_genolike_dpar(const arma::vec par, const arma::mat& pgA, const arma::mat& pgB);
RcppExport SEXP _ldsep_dllike_genolike_dpar(SEXP parSEXP, SEXP pgASEXP, SEXP pgBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pgA(pgASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pgB(pgBSEXP);
    rcpp_result_gen = Rcpp::wrap(dllike_genolike_dpar(par, pgA, pgB));
    return rcpp_result_gen;
END_RCPP
}
// dmulti_double
double dmulti_double(const arma::vec x, const arma::vec prob, bool log_p);
RcppExport SEXP _ldsep_dmulti_double(SEXP xSEXP, SEXP probSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    Rcpp::traits::input_parameter< bool >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(dmulti_double(x, prob, log_p));
    return rcpp_result_gen;
END_RCPP
}
// probgeno
double probgeno(const int& gA, const int& gB, const int K, const arma::vec prob, bool log_p);
RcppExport SEXP _ldsep_probgeno(SEXP gASEXP, SEXP gBSEXP, SEXP KSEXP, SEXP probSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type gA(gASEXP);
    Rcpp::traits::input_parameter< const int& >::type gB(gBSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    Rcpp::traits::input_parameter< bool >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(probgeno(gA, gB, K, prob, log_p));
    return rcpp_result_gen;
END_RCPP
}
// proballgeno
double proballgeno(const arma::vec& gA, const arma::vec& gB, const int K, const arma::vec prob, bool log_p);
RcppExport SEXP _ldsep_proballgeno(SEXP gASEXP, SEXP gBSEXP, SEXP KSEXP, SEXP probSEXP, SEXP log_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type gA(gASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gB(gBSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type prob(probSEXP);
    Rcpp::traits::input_parameter< bool >::type log_p(log_pSEXP);
    rcpp_result_gen = Rcpp::wrap(proballgeno(gA, gB, K, prob, log_p));
    return rcpp_result_gen;
END_RCPP
}
// llike_geno
double llike_geno(const arma::vec par, const arma::vec& gA, const arma::vec& gB, const int K);
RcppExport SEXP _ldsep_llike_geno(SEXP parSEXP, SEXP gASEXP, SEXP gBSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gA(gASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gB(gBSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(llike_geno(par, gA, gB, K));
    return rcpp_result_gen;
END_RCPP
}
// optimize_genocor
List optimize_genocor(arma::vec& par, const arma::vec& gA, const arma::vec& gB, const int& K, double reltol);
RcppExport SEXP _ldsep_optimize_genocor(SEXP parSEXP, SEXP gASEXP, SEXP gBSEXP, SEXP KSEXP, SEXP reltolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gA(gASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gB(gBSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type reltol(reltolSEXP);
    rcpp_result_gen = Rcpp::wrap(optimize_genocor(par, gA, gB, K, reltol));
    return rcpp_result_gen;
END_RCPP
}
// log_sum_exp
double log_sum_exp(const arma::vec x);
RcppExport SEXP _ldsep_log_sum_exp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_exp(x));
    return rcpp_result_gen;
END_RCPP
}
// log_sum_exp_2
double log_sum_exp_2(double x, double y);
RcppExport SEXP _ldsep_log_sum_exp_2(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_exp_2(x, y));
    return rcpp_result_gen;
END_RCPP
}
// plog_sum_exp
arma::vec plog_sum_exp(const arma::vec x, const arma::vec y);
RcppExport SEXP _ldsep_plog_sum_exp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(plog_sum_exp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// logit
double logit(double x);
RcppExport SEXP _ldsep_logit(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logit(x));
    return rcpp_result_gen;
END_RCPP
}
// expit
double expit(double x);
RcppExport SEXP _ldsep_expit(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(expit(x));
    return rcpp_result_gen;
END_RCPP
}
// real_to_simplex
arma::vec real_to_simplex(const arma::vec y);
RcppExport SEXP _ldsep_real_to_simplex(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(real_to_simplex(y));
    return rcpp_result_gen;
END_RCPP
}
// simplex_to_real
arma::vec simplex_to_real(const arma::vec x);
RcppExport SEXP _ldsep_simplex_to_real(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(simplex_to_real(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ldsep_dmulti_dprob", (DL_FUNC) &_ldsep_dmulti_dprob, 3},
    {"_ldsep_dprobgeno_dprob", (DL_FUNC) &_ldsep_dprobgeno_dprob, 4},
    {"_ldsep_dproballgeno_dprob", (DL_FUNC) &_ldsep_dproballgeno_dprob, 4},
    {"_ldsep_dreal_to_simplex_dy", (DL_FUNC) &_ldsep_dreal_to_simplex_dy, 1},
    {"_ldsep_dsimplex_to_real_dx", (DL_FUNC) &_ldsep_dsimplex_to_real_dx, 1},
    {"_ldsep_dllike_geno_dpar", (DL_FUNC) &_ldsep_dllike_geno_dpar, 4},
    {"_ldsep_dD_dprob", (DL_FUNC) &_ldsep_dD_dprob, 1},
    {"_ldsep_dr2_dprob", (DL_FUNC) &_ldsep_dr2_dprob, 1},
    {"_ldsep_dDprime_dprob", (DL_FUNC) &_ldsep_dDprime_dprob, 1},
    {"_ldsep_probgenolike", (DL_FUNC) &_ldsep_probgenolike, 4},
    {"_ldsep_proballgenolike", (DL_FUNC) &_ldsep_proballgenolike, 4},
    {"_ldsep_llike_genolike", (DL_FUNC) &_ldsep_llike_genolike, 3},
    {"_ldsep_get_dprobgeno_dprob_array", (DL_FUNC) &_ldsep_get_dprobgeno_dprob_array, 2},
    {"_ldsep_get_prob_array", (DL_FUNC) &_ldsep_get_prob_array, 2},
    {"_ldsep_dprobgenolike_dprob", (DL_FUNC) &_ldsep_dprobgenolike_dprob, 3},
    {"_ldsep_dproballgenolike_dprob", (DL_FUNC) &_ldsep_dproballgenolike_dprob, 3},
    {"_ldsep_dllike_genolike_dpar", (DL_FUNC) &_ldsep_dllike_genolike_dpar, 3},
    {"_ldsep_dmulti_double", (DL_FUNC) &_ldsep_dmulti_double, 3},
    {"_ldsep_probgeno", (DL_FUNC) &_ldsep_probgeno, 5},
    {"_ldsep_proballgeno", (DL_FUNC) &_ldsep_proballgeno, 5},
    {"_ldsep_llike_geno", (DL_FUNC) &_ldsep_llike_geno, 4},
    {"_ldsep_optimize_genocor", (DL_FUNC) &_ldsep_optimize_genocor, 5},
    {"_ldsep_log_sum_exp", (DL_FUNC) &_ldsep_log_sum_exp, 1},
    {"_ldsep_log_sum_exp_2", (DL_FUNC) &_ldsep_log_sum_exp_2, 2},
    {"_ldsep_plog_sum_exp", (DL_FUNC) &_ldsep_plog_sum_exp, 2},
    {"_ldsep_logit", (DL_FUNC) &_ldsep_logit, 1},
    {"_ldsep_expit", (DL_FUNC) &_ldsep_expit, 1},
    {"_ldsep_real_to_simplex", (DL_FUNC) &_ldsep_real_to_simplex, 1},
    {"_ldsep_simplex_to_real", (DL_FUNC) &_ldsep_simplex_to_real, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_ldsep(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
