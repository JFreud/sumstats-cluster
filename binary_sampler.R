# referenced https://stephens999.github.io/fiveMinuteStats/gibbs_structure_simple.html
# and https://stephens999.github.io/fiveMinuteStats/gibbs2.html
source('helper.R')

#' @param x an R vector of data for on variant
#' @param P a K by R matrix of effect frequencies
#' @return the log-likelihood for each of the K clusters
log_pr_x_given_P <- function(x,P){
  tP <- t(P) #transpose P so tP is R by K
  return(colSums(x*log(tP)+(1-x)*log(1-tP)))
}

normalize <- function(x){return(x/sum(x))}

#' @param z an n vector of cluster allocations (1...k)
#' @param k the number of clusters
sample_pi = function(z,k){
  counts = colSums(outer(z,1:k,FUN="=="))
  pi_vec = gtools::rdirichlet(1,counts+1)
  return(pi_vec)
}

#' @param x an n by R matrix of data
#' @param P a K by R matrix of effect frequencies
#' @param pi_vec vector of cluster proportions
#' @return an n vector of cluster memberships
sample_z <- function(x,P,pi_vec){
  K <- nrow(P)
  loglik_matrix <- apply(x, 1, log_pr_x_given_P, P=P)
  lik_matrix <- exp(loglik_matrix)
  p.z.given.x <- sweep(lik_matrix, MARGIN=1, pi_vec, `*`)
  p.z.given.x <- apply(p.z.given.x,2,normalize) # normalize lik * prior
  z <- rep(0, nrow(x))
  for(i in 1:length(z)){
    z[i] <- sample(1:K, size=1,prob=p.z.given.x[,i],replace=TRUE)
  }
  return(list(z=z, p.z.given.x=p.z.given.x))
}


#' @param x an n by R matrix of data
#' @param z an n vector of cluster allocations
#' @param k number of clusters
#' @return a K by R matrix of effect frequencies
sample_P <- function(x, z, k){
  R <- ncol(x)
  P <- matrix(nrow=k,ncol=R)
  for(i in 1:k){
    sample_size <- sum(z==i)
    if(sample_size==0){
      number_of_ones<-rep(0,R)
    } else {
      number_of_ones <- colSums(x[z==i,])
    }
    # note this is making the simplifying assumption that biomarkers are independent which they are clearly not
    P[i,] <- rbeta(R,1+number_of_ones,1+sample_size-number_of_ones)
  }
  return(P)
}

#' @param R number of biomarkers
gibbs <- function(x, K, R, niter = 100){
  # initalize z
  z <- sample(1:K,nrow(x),replace=TRUE)
  res = list(z = matrix(nrow=niter, ncol=nrow(x)),
             P = array(rep(NaN, K*R*niter), c(K, R, niter)),
             pi = matrix(nrow=niter,ncol=K),
             p.z.given.x = array(rep(NaN, nrow(x)*K*niter), c(K, nrow(x), niter)))
  res$z[1,] <- z

  for(i in 2:niter){
    P <- sample_P(x, z, K)
    pi_vec <- sample_pi(z, K)
    zsample <- sample_z(x,P,pi_vec)
    z <- zsample$z
    p.z.given.x <- zsample$p.z.given.x
    res$z[i,] <- z
    res$P[,,i] <- P
    res$p.z.given.x[,,i] <- p.z.given.x
    res$pi[i,] <- pi_vec
  }
  return(res)
}


posterior_K_raw <- function(variants_binary, assignments, P_mean, prior) {
  # refer to appendix of structure paper, not sure this is right
  log_p_xi_given_zP <- rep(0, nrow(assignments))
  for (i in 1:nrow(assignments)) {
    z <- assignments[i,]$assignment
    for (j in 1:35) {
      log_p_xi_given_zP[i] <- log_p_xi_given_zP[i] + dbinom(variants_binary[i,j], 1, prob=P_mean[z,j], log=T)
    }
  }
  mu_hat <- mean(-2 * log_p_xi_given_zP)
  sigma2_hat <- mean((-2 * log_p_xi_given_zP - mu_hat)^2)
  p_x_given_k <- exp(-mu_hat/2 - sigma2_hat/8)
  p_k_given_x <- p_x_given_k * prior
  return(p_k_given_x)
}

K_inference <- function(variants_binary, k_max) {
  k_max <- 10
  K_posteriors <- matrix(,nrow=10,ncol=2)
  for (k in 1:k_max) {
    K_posteriors[k,1] <- k
    res <- gibbs(variants_binary, K=6, R=ncol(variants_binary), niter=100)
    assignments <- get_cluster_assignments(x=variants_binary, res=res, niter=100)
    P_mean <- mat_mean(res$P[,,-(1:50)])
    K_posteriors[k,2] <- posterior_K_raw(variants_binary, assignments, P_mean, prior=1.0/k_max)
  }
  K_posteriors[,2] <- sapply(K_posteriors[,2], function(x) {x / sum(K_posteriors[,2])})
  return(K_posteriors)
}
