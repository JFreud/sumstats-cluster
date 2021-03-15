# helper functions

#' @param data summary stats dataframe
#' @param x N by R table of N variant values and R markers
#' @param res result of gibbs sampler
#' @param niter number of iterations for the gibbs sampling
get_cluster_assignments <- function(data, x, res, niter=100) {
  p_zmax <- apply(res$p.z.given.x[,,niter], 2, which.max)
  # assign to cluster with highest probability
  assignments <- data.frame(
    ID=rownames(x),
    assignment=p_zmax)
  assignments <- cbind(assignments, t(res$p.z.given.x[,,niter]))
  assignments$POS <- sapply(assignments$ID,function(x) {df[df$ID==x,]$POS[1]})
  assignments$POS <- unlist(assignments$POS)
  return(assignments)
}
