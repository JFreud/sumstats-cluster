# helper functions

#' @param x N by R table of N variant values and R markers
#' @param res result of gibbs sampler
#' @param niter number of iterations for the gibbs sampling
get_cluster_assignments <- function(x, res, niter=100) {
  assignments <- data.frame(
    ID=rownames(variants_binary),
    assignment=res$z[niter,])
  assignments <- cbind(assignments, t(res$p.z.given.x[,,niter]))
  return(assignments)
}
