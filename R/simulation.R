#########################################
# simulation function
# Naihui Zhou (ashley.n.zhou@gmail.com)
# Last updated 20191021
########################################

#' Simulate data with given parameters
#'
#' @param n Number of locations.
#' @param adj_mat \code{n*n} adjacency matrix.
#' @param alpha True parameter value for \code{alpha}.
#' @param eta True parameter value for \code{eta}
#' @param tau2 True parameter value for \code{tau2}
#' @param m Number of observations per location.
#' @param M Number of Gibbs iterations
#' @return A \code{n*m} matrix of data y
#' @examples
#' adj_trans = transform_small(irregular_lattice)
#' y = simulate_y(dim(irregular_lattice)[1], adj_trans, alpha = 2, eta = 2, tau2 = 3, m = 2, M = 2000)
#' @author Naihui Zhou (ashley.n.zhou@gmail.com)
#'
#' @export
simulate_y<-function(n, adj_mat, alpha, eta, tau2, m, M){
  w = rep(0,n)
  for (t in 1:M){
    if (t%%100==1) cat('Simulation iteration', t, '\n')
    for (i in 1:n){
      mu = alpha+eta*adj_mat[i,]%*%(w-alpha)
      #updated 20190402: if w too big rpois(1, exp(w)) will produce NA
      w[i] = rnorm(1,mu,sqrt(tau2))
      while (w[i]>21) {
        w[i] = rnorm(1,mu,sqrt(tau2))
      }
    }
  }
  lambda = exp(w)
  y = matrix(NA, nrow=n, ncol=m)
  for (i in 1:n){
    y[i,] = rpois(m,lambda[i])
  }
  return(y)
}
