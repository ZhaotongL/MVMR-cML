# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' MVMRcML method with BIC selection
#'
#' This is the main function of MVMR-cML.
#'
#' @param b_exp A m*L matrix of SNP effects on the exposure variable.
#' @param b_out A m*1 matrix of SNP effects on the outcome variable.
#' @param se_bx A m*L matrix of standard errors of \code{b_exp}.
#' @param Sig_inv_l A list of the inverse of m covariance matrices, generated by \code{invcov_mvmr()}.
#' @param n The smallest sample size of the L+1 GWAS dataset.
#' @param K_vec Sets of candidate K's, the constraint parameter representing number of invalid IVs.
#' @param random_start Number of random start points, default is 1.
#' @param maxit Maximum number of iterations for each optimization, default is 100.
#' @param thres Threshold for convergence criterion.
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @return A list
#' \describe{
#' \item{BIC_theta}{Estimated causal effect from MVMR-cML-BIC}
#' \item{BIC_invalid}{Invalid IVs selected by MVMR-cML-BIC}
#' \item{K_vec}{K_vec to be searched from}
#' \item{Converge}{0 if converge; 1 otherwise}
#' \item{BIC_vec}{BIC vector corresponding to K_vec}
#' }
#'
#' @export
MVmr_cML <- function(b_exp, b_out, se_bx, Sig_inv_l, n, K_vec = as.numeric( c()), random_start = 1L, maxit = 100L, thres = 1e-4) {
    .Call('_MVMRcML_MVmr_cML', PACKAGE = 'MVMRcML', b_exp, b_out, se_bx, Sig_inv_l, n, K_vec, random_start, maxit, thres)
}

#' MVMRcML method with Data Perturbation
#'
#' This is the main function of MVMR-cML.
#'
#' @param b_exp A m*L matrix of SNP effects on the exposure variable.
#' @param b_out A m*1 matrix of SNP effects on the outcome variable.
#' @param se_bx A m*L matrix of standard errors of \code{b_exp}.
#' @param Sig_inv_l A list of the inverse of m covariance matrices, generated by \code{invcov_mvmr()}.
#' @param n The smallest sample size of the L+1 GWAS dataset.
#' @param K_vec Sets of candidate K's, the constraint parameter representing number of invalid IVs.
#' @param random_start Number of random start points, default is 1.
#' @param num_pert Number of perturbation, default is 100.
#' @param maxit Maximum number of iterations for each optimization, default is 100.
#' @param thres Threshold for convergence criterion.
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @return A list
#' \describe{
#' \item{BIC_theta}{Estimated causal effect from MVMR-cML-BIC}
#' \item{BIC_invalid}{Invalid IVs selected by MVMR-cML-BIC}
#' \item{BIC_DP_theta}{Estimated causal effect from MVMR-cML-DP }
#' \item{BIC_DP_se}{Estimate standard error for \code{BIC_DP_theta}}
#' \item{eff_DP_B}{Data perturbation with successful convergence}
#' }
#'
#' @export
MVmr_cML_DP <- function(b_exp, b_out, se_bx, Sig_inv_l, n, K_vec = as.numeric( c()), random_start = 1L, num_pert = 100L, maxit = 100L, thres = 1e-4) {
    .Call('_MVMRcML_MVmr_cML_DP', PACKAGE = 'MVMRcML', b_exp, b_out, se_bx, Sig_inv_l, n, K_vec, random_start, num_pert, maxit, thres)
}

