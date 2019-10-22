##################################
# Post MCMC analysis
# Functions to print parameter estimates
# Naihui Zhou (ashley.n.zhou@gmail.com)
# Last updated 20191021
##################################

#' Percentage of jumps in MCMC iterations
#'
#' @param ret returned list from \code{pmrf()}
#' @param Ti number of total iterations
#' @param N datasize
#' @return A list with jump percentages for \code{w}, \code{alpha}, \code{eta} and \code{tau2}
#' @examples
#' adj_trans = transform_small(irregular_lattice)
#' n = dim(irregular_lattice)[1]
#' y = simulate_y(n = n, adj_mat = adj_trans, alpha = 2, eta = 2, tau2 = 3, m = 2, M = 2000)
#' bounds_e = get_eta_param_space_small(adj_trans)
#' ret = pmrf(total_iter = 1000, N = n , y, adj_trans, c(1, 0.5, NA, 2), bounds_e)
#' get_jump_frequency(ret, 1000, n)
#' @author Naihui Zhou (ashley.n.zhou@gmail.com)
#'
#' @export
get_jump_frequency<-function(ret, Ti, N){
  #ret is a list of returned values from c
  #Ti is total number of iterations
  #N is data size

  jumpcount = ret[[5]]
  rate_w = jumpcount[1]/(N*Ti)
  rate_alpha = jumpcount[2]/Ti
  rate_eta = jumpcount[3]/Ti
  rate_tau2 = jumpcount[4]/Ti

  return(list(w = rate_w, alpha = rate_alpha, eta = rate_eta, tau2 = rate_tau2))
}

#' Delete burn in iterations in MCMC
#'
#' @param ret returned list from \code{pmrf()}
#' @param burn_in number of iterations to throw away
#' @return A list with MCMC iterations after burn in for \code{w}, \code{alpha}, \code{eta} and \code{tau2}
#' @examples
#' adj_trans = transform_small(irregular_lattice)
#' n = dim(irregular_lattice)[1]
#' y = simulate_y(n = n, adj_mat = adj_trans, alpha = 2, eta = 2, tau2 = 3, m = 2, M = 2000)
#' bounds_e = get_eta_param_space_small(adj_trans)
#' ret = pmrf(total_iter = 1000, N = n , y, adj_trans, c(1, 0.5, NA, 2), bounds_e)
#' newret = delete_burn_in(ret, 200)
#' @author Naihui Zhou (ashley.n.zhou@gmail.com)
#'
#'
#' @export
delete_burn_in <-function(ret, burn_in){
  total_iter = length(ret$alpha)
  new_w = ret$w[,(burn_in:total_iter)]
  new_alpha = ret$alpha[burn_in:total_iter]
  new_eta = ret$eta[burn_in:total_iter]
  new_tau2 = ret$tau2[burn_in:total_iter]
  return(list(w=new_w, alpha=new_alpha, eta=new_eta, tau2=new_tau2))
}

#' Print parameter estimates to console
#'
#' Parameter estimates include mean and 95\% credible interval.
#'
#' It is STRONGLY recommended that you output the estimates after throwing away burn-in iterations.
#'
#' @param ret returned list from \code{pmrf()}
#' @param burn_in number of iterations to throw away. Default 200.
#' @param digit number of digits after the decimal point to display. Default 3.
#' @examples
#' adj_trans = transform_small(irregular_lattice)
#' n = dim(irregular_lattice)[1]
#' y = simulate_y(n = n, adj_mat = adj_trans, alpha = 2, eta = 2, tau2 = 3, m = 2, M = 2000)
#' bounds_e = get_eta_param_space_small(adj_trans)
#' ret = pmrf(total_iter = 1000, N = n , y, adj_trans, c(1, 0.5, NA, 2), bounds_e)
#' print_param_estimates(ret, 300, 4)
#' print_param_estimates(ret, 0, 2)
#' @author Naihui Zhou (ashley.n.zhou@gmail.com)
#'
#' @export
print_param_estimates<-function(ret, burn_in = 200, digit = 3){
  newret = delete_burn_in(ret, burn_in)
  alpha_percentile = quantile(newret$alpha, probs=c(0.025, 0.975))
  cat(paste(round(mean(newret$alpha),digit), " (", round(alpha_percentile[1],digit), ", ", round(alpha_percentile[2],digit), ") ", sep=""))
  cat("\n")
  eta_percentile = quantile(newret$eta, probs=c(0.025, 0.975))
  cat(paste(round(mean(newret$eta),digit), " (", round(eta_percentile[1],digit), ", ", round(eta_percentile[2],digit), ") ", sep=""))
  cat("\n")
  tau2_percentile = quantile(newret$tau2, probs=c(0.025, 0.975))
  cat(paste(round(mean(newret$tau2),digit), " (", round(tau2_percentile[1],digit), ", ", round(tau2_percentile[2],digit), ") ", sep=""))
  cat("\n")
}
