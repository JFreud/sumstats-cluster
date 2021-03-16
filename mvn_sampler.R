# referenced https://stephens999.github.io/fiveMinuteStats/gibbs_structure_simple.html
# and https://stephens999.github.io/fiveMinuteStats/gibbs2.html


library(mvtnorm)
library(LaplacesDemon)
library(gtools)

#' @param x an R vector of data for on variant
#' @param mu a K * R matrix of MVN mean parameters
#' @param Sigma a length K * R^2 matrix of MVN covariance parameters
#' @return the log-likelihood for each of the K clusters

log_pr_x_given_P <- function(x, mu, Sigma){
    K <- nrow(mu)
    R <- ncol(mu)
    ret <- vector(length = K)
    for(j in (1:K)){
        Sigma.mat <- matrix(Sigma[j,], nrow = R, ncol = R)
        ret[j] <- log(dmvnorm(x, mu[j,], Sigma.mat))
    }
    return(ret)
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
#' @param mu a K x R matrix of MVN mean params
#' @param Sigma a K x R^2 matrix of MVN cov params
#' @param pi_vec vector of cluster proportions
#' @return an n vector of cluster memberships

sample_z <- function(x, mu, Sigma, pi_vec){
    n <- nrow(x)
    K <- length(pi_vec)
    R <- ncol(x)
    
    loglik_matrix <- apply(x, 1, log_pr_x_given_P, mu=mu, Sigma=Sigma)
    lik_matrix <- exp(loglik_matrix)
    p.z.given.x <- sweep(lik_matrix, MARGIN=1, pi_vec, `*`)
    #p.z.given.x <- apply(p.z.given.x,2,normalize) # normalize lik * prior
    
    z <- rep(0, n)
    for(i in 1:n){
        #print(p.z.given.x[,i])
        z[i] <- sample(1:K, size=1,prob=p.z.given.x[,i],replace=TRUE)
    }
    return(list(z=z, p.z.given.x=p.z.given.x))
}

#' @param x an n by R matrix of data
#' @param z an n vector of cluster allocations
#' @param K number of clusters
#' @return tuple with K x R matrix of means and K x R^2 matrix of covs
sample_params <- function(x, z, K){
    R <- ncol(x)
    mu <- matrix(nrow = K, ncol = R)
    Sigma <- matrix(nrow = K, ncol = R * R)
    
    for(i in 1:K){
        obs <- x[z == i, ]
        sample_size <- sum(z == i)
        
        if(sample_size == 0){
            xbar <- matrix(rep(0, R), nrow = 1)
            S <- matrix(0, nrow = R, ncol = R)
        }
        else if(sample_size == 1 ){
            xbar <- matrix(x[i, ], nrow = 1)
            S <- matrix(0, nrow = R, ncol = R)
        }
        else{
            xbar <- matrix(apply(obs, 2, FUN = mean), nrow =1 )
            diff <- t(apply(obs, 1, FUN = function(x) (x - xbar)))
            S <- t(diff) %*% diff
        }
        
        mu.n <- (1 + sample_size * xbar)/(R + sample_size)
        lambda.n <- (1 + sample_size)
        nu.n <- (R + sample_size)
        xbar <- matrix(xbar, nrow = 1)
        Psi.n <- diag(R) + S + sample_size/(1 + sample_size) * t(xbar) %*% xbar
        param <- rnorminvwishart(mu0 = mu.n,
                                 lambda = lambda.n,
                                 S = Psi.n, 
                                 nu = nu.n)
        mu[i,] <- param$mu
        Sigma[i,] <- as.vector(param$Sigma)
    }
    return(list(mu = mu, Sigma = Sigma))
}

#' @param R number of biomarkers
#' @param K number of clusters
#' @param x the data
gibbs <- function(x, K, niter = 100){
    # initalize z
    R <- ncol(x)
    z <- sample(1:K,nrow(x),replace=TRUE)
    res = list(z = matrix(nrow=niter, ncol=nrow(x)),
               mu = array(rep(NaN, K*R*niter), c(K, R, niter)),
               Sigma = array(rep(NaN, K*R^2*niter), c(K, R^2, niter)),
               pi = matrix(nrow=niter,ncol=K))
    res$z[1,] <- z
    
    for(i in 2:niter){
        if (i %% 10 == 0) {print(i)}
        params <- sample_params(x, z, K)
        mu <- params$mu
        Sigma <- params$Sigma
        
        pi_vec <- sample_pi(z, K)
        zsample <- sample_z(x, mu, Sigma, pi_vec)
        z <- zsample$z
        #print(round(zsample$p.z.given.x, 3))
        
        res$z[i,] <- z
        res$mu[,,i] <- mu
        res$Sigma[,,i] <- Sigma
        res$pi[i,] <- pi_vec
        
        #print(pi_vec) 
    }
    return(res)
}
