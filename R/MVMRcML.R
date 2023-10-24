#' Profile likelihood of valid IVs
#'
#' calculate the profile likelihood of valid IVs up to a constant
#'
#' @keywords internal
#'
pl <- function(x,b_exp_v,b_out_v,Sig_inv_v){
  k = ncol(b_exp_v)
  m_valid = length(b_out_v)
  pll = 0
  for(i in 1:m_valid){
    W = Sig_inv_v[[i]]
    b_exp_i = b_exp_v[i,]
    b_out_i = b_out_v[i]
    beta = c(b_exp_i,b_out_i)
    B = W[-(k+1),] %*% beta + c(W[(k+1),] %*% beta) * x
    A = W[-(k+1),-(k+1)] + W[-(k+1),k+1] %*% t(x) + x %*% t(W[k+1,-(k+1)]) + W[k+1,k+1] * x %*% t(x)
    bhat_xi = solve(A) %*% B
    b_i = c(bhat_xi,t(bhat_xi) %*% x)
    pll = pll  -1/2 * t(b_i) %*% W %*% b_i + t(beta) %*% W %*% b_i
  }
  return(-pll)
}

#' Generate the list of inverse of covariance matrices used in \code{MVMR-cML-DP}
#'
#' @param se_bx A m*L matrix of standard errors of SNP-exposure association
#' @param se_by A vector of standard errors of SNP-outcome association
#' @param rho_mat A (L+1)\times(L+1) correlation matrix among the L exposures and the outcome, note that the upper-left L\times L submatrix corresponding to the L exposures.
#' @return A list of inverse of covariance matrices with respect to each genetic variant, retaining the ordering in \code{se_bx}
#'
#' @export
#'
invcov_mvmr <- function(se_bx, se_by, rho_mat){
  Sig_l = Sig_inv_l = list()
  m = nrow(se_bx)
  for(i in 1:m){
    Sig = rho_mat*crossprod(t(c(se_bx[i,],se_by[i])))
    Sig_l[[i]] = Sig
    Sig_inv_l[[i]] = solve(Sig)
  }
  return(Sig_inv_l)
}

#' Standard error estimate for MVMR-cML-BIC
#'
#' This is based on the profile likelihood of the set of valid IVs, which is not robust to uncertainty in model selection.
#'
#' @param b_exp A matrix of SNP effects on the exposure variable.
#' @param b_out A vector of SNP effects on the outcome variable.
#' @param Sig_inv_l A list of inverse of covariance matrix.
#' @param theta A vector of final estimates of causal effect of each exposure by MVMR-cML-BIC obtained from \code{MVmr_cML_DP}.
#' @param zero_ind A vector of the index of valid IVs.
#' @param r_vec A vector of estimated horizontal pleiotropic effects.
#'
#'
#' @import numDeriv
#' @return A vector
#'
#' @export
#'
MVcML_SdTheta <- function(b_exp,b_out,Sig_inv_l,theta,zero_ind,r_vec=NULL){
  if(!is.null(r_vec)){
      zero_ind = which(r_vec==0)
    }
  H = hessian(pl,x=theta,b_exp_v = b_exp[zero_ind,,drop=FALSE], b_out_v=b_out[zero_ind],Sig_inv_v =Sig_inv_l[zero_ind])
  se_theta = sqrt(diag(solve(H)))
  # ores = optim(theta,pl,method = 'BFGS',b_exp_v=b_exp[zero_ind,],b_out_v=b_out[zero_ind],Sig_inv_v= Sig_inv_l[zero_ind],
  #              hessian = TRUE)
  # sqrt(diag(solve(ores$hessian)))
  return(se_theta)
}


