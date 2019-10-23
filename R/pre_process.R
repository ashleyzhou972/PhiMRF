##################################
# Preprocess to the MCMC analysis
# Functions to transform adjacency matrices and
# calculate max and min eigen values
# Naihui Zhou (ashley.n.zhou@gmail.com)
# Last Updated 20191018
##################################

#' Transform adjacency matrix by dividing sum of row sum and column sum
#'
#' @details  This is for small matrices that R can handle with its memory.
#' If your matrix is big (n >10000), consider the \code{preprocess_big} function.
#' @param adj_mat A symmetric adjacency matrix.
#' Edges are represented as 1 and non edges 0.
#' @return A transformed matrix in the \code{dgCMatrix} class
#' @examples
#' mat = matrix(c(1,1,0,1,1,0,0,0,1), nrow = 3)
#' transform_small(mat)
#' @author Naihui Zhou (ashley.n.zhou@gmail.com)
#'
#' @export
transform_small<-function(adj_mat){
  if (!Matrix::isSymmetric(adj_mat)) {
    warning("Input adjacency matrix not symmetric\n", immediate. = T)
  }
	dim = dim(adj_mat)[1]
	rowSum_mat = matrix(rep(Matrix::rowSums(adj_mat),dim), byrow=F, nrow=dim)
	colSum_mat = matrix(rep(Matrix::colSums(adj_mat),dim), byrow=T, nrow=dim)
	sum_mat = rowSum_mat + colSum_mat
	sum_mat[sum_mat==0]<--Inf
	new_neighbor = as(adj_mat/sum_mat, "dgCMatrix")
	return(new_neighbor)
}

#' Compute the parameter space of eta from the adjacency matrix
#'
#' @description The bounds of eta is decided by the reciprocal of the maximum and minimum of eigenvalues of the adjacency matrix.
#' @details  This is for small matrices that R can handle with its memory.
#' If your matrix is big (n >10000), consider the \code{preprocess_big} function.
#' The input matrix should be transformed using \code{transform_small()} or \code{preprocess_big}
#' to be used in \code{pmrf()} directly as the \code{bounds_e} argument.
#' @param adj_mat A (transformed) symmetric adjacency matrix.
#' @return A vector of the lower and upper bounds for eta
#' @examples
#' mat = matrix(c(1,1,0,1,1,0,0,0,1), nrow = 3)
#' get_eta_param_space_small(mat)
#' @author Naihui Zhou (ashley.n.zhou@gmail.com)
#'
#' @export
get_eta_param_space_small<-function(adj_mat){
  if (!Matrix::isSymmetric(adj_mat)) {
    warning("Input adjacency matrix not symmetric\n", immediate. = T)
  }
  evalues = eigen(adj_mat, symmetric = T, only.values=T)$values
  min = min(evalues)
  max = max(evalues)
  if (min<0 & max>0){
    eta_l = 1/min
    eta_u = 1/max
  } else {
    eta_l = 1/max
    eta_u = 1/min
  }
  cat("Bounds for eta is (", eta_l, ",", eta_u, ")\n")
  return(c(eta_l, eta_u))
}

#' Matrix transformation and eigenvalue calculation for large matrices
#'
#' @description This function has the combined function of
#'  \code{transform_small()} and \code{get_eta_param_space_small}
#'
#' @details  This function is designed for large matrices that R cannot handle.
#'
#' @param adj_mat An \strong{untransformed} symmetric adjacency matrix.
#' @param savepath The path to save the transformed matrix and bounds_eta.
#' If \code{NULL} or missing, then not saving.
#' @return A list of \code{b_eta}: parameter space of eta and \code{new_mat}: transformed matrix.
#' @author Naihui Zhou (ashley.n.zhou@gmail.com)
#'
#' @useDynLib PhiMRF preprocess
#' @export
preprocess_big<-function(adj_mat, savepath){
  if (class(adj_mat)!='dgCMatrix') {
    stop("Please input neighborhood matrix in 'dgCMatrix' class.\n Use as(yourmatrix, 'dgCMatrix') in the Matrix package to convert\n")
  }
  if (!Matrix::isSymmetric(adj_mat)){
    stop("Please input symmetric matrix\n")
  }
  att = attributes(adj_mat)
  dim = att$Dim[1]
  print(is.integer(dim))
  val = att$x
  row_ind = att$i
  col_ptr = att$p
  ret = .Call("preprocess", dim_in = as.integer(dim),
              val_in = as.double(val), row_ind_in = as.integer(row_ind),
              col_ptr_in = as.integer(col_ptr))
  b_eta = ret[[2]]
  new_val = ret[[1]]
  new_mat = Matrix::sparseMatrix(i = row_ind, p = col_ptr, x = new_val, dims = att$Dim, index1 = F, giveCsparse = T)
  new_mat = as(new_mat, 'dgCMatrix')
  if (!(missing(savepath) || is.null(savepath))){
    saveRDS(b_eta, paste0(savepath, "/bounds_eta.rds"))
    saveRDS(new_mat, paste0(savepath, "/adj_trans.rds"))
  }
  return(list(b_eta = b_eta, new_mat = new_mat))
}


