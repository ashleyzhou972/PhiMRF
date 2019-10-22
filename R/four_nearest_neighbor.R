################################################
# Generate adjacency matrix for regular lattice with four nearest neighbors
# Naihui Zhou (ashley.n.zhou@gmail.com)
# Last Updated 20191021
##################################################


#' Generate adjacency matrix for a four-nearest-neighbor regular lattice
#'
#' @description This function generates the adjacency matrix for a regular rectangular lattice.
#' Each location on the lattice has four neighbors on orthogonal directions
#' The neighbors wraps around on edges of the lattice.
#'
#' @param v Number of locations in the rectangular field; Must be an even number
#' @param k Number of rows in the rectangular field with a total of \code{v} locations.
#'
#' @details for example if \code{v = 30},\code{k = 5}, then the indices of the locations are as follows.
#'
#' 1  6 11 16 21 26
#'
#' 2  7 12 17 22 27
#'
#' 3  8 13 18 23 28
#'
#' 4  9 14 19 24 29
#'
#' 5 10 15 20 25 30
#'
#' @return  A v*v adjacency matrix
#' @author Naihui Zhou (ashley.n.zhou@gmail.com)
#'
#' @export
generate_four_nearest_neighbor_matrix<-function(v,k){
  mat = matrix(0, nrow = v, ncol = v)
  #four vertex cases (1, k, 1+v-k and v)
  #four edge cases ([1,k], [v-k+1,v],[i%%k==1], [i%%k==1])
  #rest are internal cases
  pointy = c(1,k,1+v-k,v)
  edgeL = c(1:k)
  edgeR = c((v-k+1):v)
  rlength = v/k
  edgeT = rep(NA, rlength)
  edgeB = rep(NA, rlength)
  for (j in 1:rlength){
    edgeT[j] = 1+(j-1)*k
    edgeB[j] = j*k
  }
  edges = c(edgeL, edgeR, edgeT, edgeB)
  for (i in 1:v){
    if (i %in% edges){
      if (i %in% edgeL){
        mat[i,i+k] = 1
        mat[i,i-k+v] = 1
        if (!(i %in% pointy)){
          mat[i,i+1] = 1
          mat[i,i-1] = 1
        }
      }
      if (i %in% edgeR){
        mat[i,i+k-v] = 1
        mat[i,i-k] = 1
        if (!(i %in% pointy)){
          mat[i,i+1] = 1
          mat[i,i-1] = 1
        }
      }
      if (i %in% edgeT){
        mat[i,i+1] = 1
        mat[i,i-1+k] = 1
        if (!(i %in% pointy)){
          mat[i,i+k] = 1
          mat[i,i-k] = 1
        }
      }
      if (i %in% edgeB){
        mat[i,i+1-k] = 1
        mat[i,i-1] = 1
        if (!(i %in% pointy)){
          mat[i,i+k] = 1
          mat[i,i-k] = 1
        }
      }
    }
    else {
      mat[i,i+1] = 1
      mat[i,i-1] = 1
      mat[i,i+k] = 1
      mat[i,i-k] = 1
    }
  }
  return(mat)
}
