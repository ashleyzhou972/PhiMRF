############################################
# R wrapper function calling for C
# To perform the main pmrf MCMC inference
# updated 20191016
############################################

#' Poisson mixed Markov Random Field model inference
#'
#' @description A function to perform a double Metropolis-Hastings algorithm for
#' the parameter estimates of spatial count data observed on a network under a Poisson hierarchical Markov Random Field model.
#' The parameters are \code{w}, \code{alpha}, \code{eta} and \code{tau2}. \code{alpha} is connected with the conditional mean,
#' \code{eta} is connected with the level of spatial dependency and \code{tau2} is connected with conditional variance.
#' See References for detailed model specification.
#'
#' @param total_iter Total number of MCMC iterations.
#' @param N Number of locations that data is observed. Multiple observations can be made at each location.
#' @param y Observed data at each location. \code{y} should be a matrix with \code{N} rows, and number of columns equal to the number of observations at each location.
#' @param adj_mat Transformed adjacency matrix for the locations. See Details on how to transform the adjacency matrix.
#' Must be in the 'dgCMatrix' class. Use \code{as(yourmatrix, 'dgCMatrix')} to convert.
#' @param vars Variances for the random walk proposal for the MCMC for each parameter,
#' \code{w}, \code{alpha}, \code{eta} and \code{tau2}.
#' @param bounds_e Lower and upper bound for \code{eta}. This should be decided by the eigenvalues of the adjacency matrix. See details. No default supplied.
#' @param bounds_a Lower and upper bound for \code{alpha}. Default \code{c(-10,10)}.
#' @param bounds_t Lower and upper bound for \code{tau} (Not \code{tau2}, see Details for the prior distributions.)
#' Default \code{c(0,10)}.
#' @param inis Initial values for \code{alpha}, \code{eta} and \code{tau2}. Default \code{c(0.1,0.0,0.1)}.
#' @param wInis Initial values for \code{w}. Must be a vector of length \code{N}. Default \code{rnorm(N, inis[1], sqrt(inis[3]))}.
#'
#' @return A list with the posterior results.  \code{w} is a \code{N} by \code{total_iters} matrix,
#' \code{alpha}, \code{eta} and \code{tau2} are vectors with length \code{N} respectively.
#' Last element in the list is a vector of length four that counts the number of jumps made for each parameter.
#'
#' @examples
#' pmrf(5000, nrow(y), y, adj_mat, c(0.5,0.5,NA,1), c(-2,2))
#' pmrf(2000, nrow(y), y, adj_mat, c(0.5,0.5,NA,1), c(-2,2), inis = c(0,0,0), wInis = rnorm(nrow(y), 0, 1))
#'
#' @details The adjacency matrix \code{adj_mat} needs to be transformed by dividing each nonzero element (default is \code{1}) by the sum of the row sum and column sum.
#' See  and  for  ways to transform.
#'
#' The parameter space of \code{eta} is defined as the reciprocal of the maximum and minimum of the eigenvalues of the (transformed) adjacency matrix.
#' (i.e. \code{c(1/min(eigen(adj_mat)$values), 1/max(eigen(adj_mat)$values))})
#' See and for ways to compute.
#'
#' The \code{vars} argument supplies the variance used in the random walk proposal.
#' In general, the larger the variance, the smaller the jump frequency.
#'
#' For \code{eta}, an independence proposal was used. So the third element of \code{vars} is never used.
#'
#' For \code{alpha} and \code{eta}, uniform priors were with the specified bounds were used for each.
#' For \code{tau2}, a uniform prior is specified for \code{tau},
#' and the prior distribution for \code{tau2} is derived using change-of-variables.
#'
#' You can supply the inital values using the last iteration of a previous MCMC run to continue to process.
#'
#' @author Naihui Zhou (ashley.n.zhou@gmail.com)
#' @references
#' Zhou, N., Friedberg, I. & Kaiser, M.S. (2019).  ""
#'
#' Liang, F. (2010). A double Metropolis–Hastings sampler for spatial models with intractable normalizing constants. \emph{Journal of Statistical Computation and Simulation}, 80(9), 1007-1022.
#'
#' @useDynLib PhiMRF double_metropolis_cont
#' @export
pmrf<-function(total_iter, N, y, adj_mat, vars, bounds_e, bounds_a=c(-10,10), bounds_t=c(0,10), inis=c(0.1,0.0,0.1), wInis=rnorm(N, inis[1], sqrt(inis[3]))){
	if (!is.numeric(total_iter) || !is.numeric(y) || !is.numeric(N) || !is.numeric(vars) || !is.numeric(bounds_a) || !is.numeric(bounds_e) || !is.numeric(bounds_t) || !is.numeric(inis) || !is.numeric(wInis)){
		stop("Input data not numeric\n")
	}
	else if (class(adj_mat)!='dgCMatrix') {
		stop("Please input neighborhood matrix in 'dgCMatrix' class.\n Use as(yourmatrix, 'dgCMatrix') in the Matrix package\n")
	}
  att = attributes(adj_mat)
	if (att$Dim[1]!=att$Dim[2]){
		stop("Input matrix not square\n")
	}
  dim = att$Dim[1]
  val = att$x
  if (sum(val) == length(val)){##if all non-zero values are 1 then matrix probably not transformed
    warning("Input matrix might not be transformed\n", immediate. = T)
  }
  if (is.null(dim(y))){
    warning("Input y should be a matrix\n")
    y = as.matrix(y, ncol = 1)
  }
  if (N!=dim(y)[1]) {
    stop("Size of input y does not agree with input N\n")
  }
  if (length(inis)!=3){
    stop("Please input 3 initial values for the three parameters, alpha, eta and tau2\n")
  }
  if (length(wInis)!=N){
    stop("Input initial values for w has incorrect size\n")
  }
  row_ind = att$i
  col_ptr = att$p
	ret = .Call(double_metropolis_cont, T_in = as.integer(total_iter),
		    N_in = as.integer(N), y_in = as.double(y) ,
		    dim_in = as.integer(dim),
		    val_in = as.double(val), row_ind_in = as.integer(row_ind),
		    col_ptr_in = as.integer(col_ptr), vars_in =as.double(vars),
		    bounds_alpha = as.double(bounds_a),
		    bounds_eta = as.double(bounds_e),
		    bounds_tau2 = as.double(bounds_t),
		    initials = as.double(inis), wInitials = as.double(wInis))
	return(list(w = ret[[1]], alpha = ret[[2]], eta = ret[[3]],
		    tau2 = ret[[4]], jump_count = ret[[5]]))
}




