// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fit_mingle
List fit_mingle(const int N_iterations, const NumericVector mu, arma::vec phi, arma::vec alpha, arma::vec beta, arma::vec p1, arma::vec sigsq, arma::vec c, const NumericVector w, const double alpha_p1, const double beta_p1, const double alpha_phi, const double beta_phi, const double alpha_alpha, const double beta_alpha, const double mu_beta, const double sigsq_beta, const double a_sigsq, const double b_sigsq, const double a_c, const double b_c, const double beta_p1_tune, const double beta_phi_tune, const double sigsq_alpha_tune);
RcppExport SEXP mingle_fit_mingle(SEXP N_iterationsSEXP, SEXP muSEXP, SEXP phiSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP p1SEXP, SEXP sigsqSEXP, SEXP cSEXP, SEXP wSEXP, SEXP alpha_p1SEXP, SEXP beta_p1SEXP, SEXP alpha_phiSEXP, SEXP beta_phiSEXP, SEXP alpha_alphaSEXP, SEXP beta_alphaSEXP, SEXP mu_betaSEXP, SEXP sigsq_betaSEXP, SEXP a_sigsqSEXP, SEXP b_sigsqSEXP, SEXP a_cSEXP, SEXP b_cSEXP, SEXP beta_p1_tuneSEXP, SEXP beta_phi_tuneSEXP, SEXP sigsq_alpha_tuneSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const int >::type N_iterations(N_iterationsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigsq(sigsqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha_p1(alpha_p1SEXP);
    Rcpp::traits::input_parameter< const double >::type beta_p1(beta_p1SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha_phi(alpha_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type beta_phi(beta_phiSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha_alpha(alpha_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta_alpha(beta_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_beta(mu_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type sigsq_beta(sigsq_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type a_sigsq(a_sigsqSEXP);
    Rcpp::traits::input_parameter< const double >::type b_sigsq(b_sigsqSEXP);
    Rcpp::traits::input_parameter< const double >::type a_c(a_cSEXP);
    Rcpp::traits::input_parameter< const double >::type b_c(b_cSEXP);
    Rcpp::traits::input_parameter< const double >::type beta_p1_tune(beta_p1_tuneSEXP);
    Rcpp::traits::input_parameter< const double >::type beta_phi_tune(beta_phi_tuneSEXP);
    Rcpp::traits::input_parameter< const double >::type sigsq_alpha_tune(sigsq_alpha_tuneSEXP);
    __result = Rcpp::wrap(fit_mingle(N_iterations, mu, phi, alpha, beta, p1, sigsq, c, w, alpha_p1, beta_p1, alpha_phi, beta_phi, alpha_alpha, beta_alpha, mu_beta, sigsq_beta, a_sigsq, b_sigsq, a_c, b_c, beta_p1_tune, beta_phi_tune, sigsq_alpha_tune));
    return __result;
END_RCPP
}
