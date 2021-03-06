% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_analysis.R
\name{delete_burn_in}
\alias{delete_burn_in}
\title{Delete burn in iterations in MCMC}
\usage{
delete_burn_in(ret, burn_in)
}
\arguments{
\item{ret}{returned list from \code{pmrf()}}

\item{burn_in}{number of iterations to throw away}
}
\value{
A list with MCMC iterations after burn in for \code{w}, \code{alpha}, \code{eta} and \code{tau2}
}
\description{
Delete burn in iterations in MCMC
}
\examples{
adj_trans = transform_small(irregular_lattice)
n = dim(irregular_lattice)[1]
y = simulate_y(n = n, adj_mat = adj_trans, alpha = 2, eta = 2, tau2 = 3, m = 2, M = 2000)
bounds_e = get_eta_param_space_small(adj_trans)
ret = pmrf(total_iter = 1000, N = n , y, adj_trans, c(1, 0.5, NA, 2), bounds_e)
newret = delete_burn_in(ret, 200)
}
\author{
Naihui Zhou (ashley.n.zhou@gmail.com)
}
