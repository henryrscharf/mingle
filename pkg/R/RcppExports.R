# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' fit the model using mcmc
#'
#' @param N_iterations if fitting all parameters, this takes approximately 1
#' minutes per 100 iterations on my laptop.
#' @param mu this is an array of dimension \code{c(n.indiv, T, 2, K)} where
#' \code{K} is some number of draws from the distribution of \eqn{s|\mu}.
#' @param phi vector of length \code{N_iterations}
#' @param alpha vector of length \code{N_iterations}
#' @param beta vector of length \code{N_iterations}
#' @param p1 vector of length \code{N_iterations}
#' @param sigsq vector of length \code{N_iterations}
#' @param c vector of length \code{N_iterations}
#' @param w array of dimension \code{(n.indiv, n.indiv, T, N_iterations)}
#' @param alpha_p1 hyperparameters
#' @param beta_p1 hyperparameters
#' @param alpha_phi hyperparameters
#' @param beta_phi hyperparameters
#' @param alpha_alpha hyperparameters
#' @param beta_alpha hyperparameters
#' @param mu_beta hyperparameters
#' @param sigsq_beta hyperparameters
#' @param a_sigsq hyperparameters
#' @param b_sigsq hyperparameters
#' @param a_c hyperparameters
#' @param b_c hyperparameters
#' @param beta_p1_tune tuning parameter for adjusted beta proposal
#' @param beta_phi_tune tuning parameter for beta proposal
#' @param sigsq_alpha_tune tuning parameter for normal proposal
#' @return list of chains for each parameter
#' @useDynLib mingle
#' @importFrom Rcpp sourceCpp
#' @export
fit_mingle <- function(N_iterations, mu, phi, alpha, beta, p1, sigsq, c, w, alpha_p1, beta_p1, alpha_phi, beta_phi, alpha_alpha, beta_alpha, mu_beta, sigsq_beta, a_sigsq, b_sigsq, a_c, b_c, beta_p1_tune, beta_phi_tune, sigsq_alpha_tune) {
    .Call('mingle_fit_mingle', PACKAGE = 'mingle', N_iterations, mu, phi, alpha, beta, p1, sigsq, c, w, alpha_p1, beta_p1, alpha_phi, beta_phi, alpha_alpha, beta_alpha, mu_beta, sigsq_beta, a_sigsq, b_sigsq, a_c, b_c, beta_p1_tune, beta_phi_tune, sigsq_alpha_tune)
}

