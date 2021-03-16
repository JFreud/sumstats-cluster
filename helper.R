# helper functions

#' @param data summary stats dataframe
#' @param x N by R table of N variant values and R markers
#' @param res result of gibbs sampler
#' @param niter number of iterations for the gibbs sampling
get_cluster_assignments <- function(data, x, res, K = 6, niter=100) {
  z_assignments <- apply(res$z[-(1:50),], 2, get_mode)
  # assign to cluster with highest probability
  assignments <- data.frame(
    ID=rownames(x),
    assignment=z_assignments)
  # get samples from posterior for z and proportions of cluster allocations
  z_prob <- apply(res$z[-(1:50),], 2, function(x) {table(factor(x, levels=1:K))})/(niter-50)
  assignments <- cbind(assignments, t(z_prob))
  assignments$POS <- sapply(assignments$ID,function(x) {df[df$ID==x,]$POS[1]})
  assignments$POS <- unlist(assignments$POS)
  return(assignments)
}


get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


mat_mean <- function(mat) {
  X <- list(mat)
  Y <- do.call(cbind, X)
  Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
  return(apply(Y, c(1, 2), mean, na.rm = TRUE))
}
