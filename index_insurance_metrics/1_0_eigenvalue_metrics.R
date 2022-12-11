#' Figenvalue of the correlation matrix
#' 
#' 
#' @param df_w dataset in a wide format: T rows and N columns, for each field
#' @param to_df return data as data.frame?
#' @param cross Whether to use the cross trick?
mtr_r2_optim <- function(df_w, to_df=FALSE, cross=FALSE){
  
  
  if(cross){
    X_scale_col <- scale(df_w)
    S_cross <- mtr_cor_cross(X_scale_col, is_scaled=TRUE)
    out <- mtr_eigenvec_cross(X_scale_col, S_cross, return_lambda = TRUE)
    vec <- drop(out$vec)
    lambda_1 <- out$lambda
  } else {
    eig_dec <- RSpectra::eigs_sym(cor(df_w), k=1)
    vec <- eig_dec$vectors[,1]
    lambda_1 <- eig_dec$values[1]
  }
  res <-  lambda_1* vec^2
  
  ##  to df?
  if(to_df){
    nam <- colnames(R)
    if(is.null(nam)) nam <- paste0("id_", 1:ncol(X))
    res <- tibble(plot_id=nam, r_2=res)
  }
  
  res
}


## When we want only first eigens, use XX'/N (T x T) instead of X'X/N (N x N)!
## XX' in N x p notation, i.e. X'X in Ishi Yata et al 2014
mtr_cov_cross <- function(X, is_scaled=FALSE){
  n <- nrow(X)
  if(!is_scaled) X <- scale(X, scale=FALSE)
  tcrossprod(X)/(n-1)
}

## same for cor!
mtr_cor_cross <- function(X, is_scaled=FALSE){
  n <- nrow(X)
  if(!is_scaled) X <- scale(X, scale=TRUE)
  tcrossprod(X)/(n-1)
}

## first eigenvector cross, from Wang and Fan 2017, (2.1) p. 1348
mtr_eigenvec_cross <- function(X, S_cross = mtr_cov_cross(X), return_lambda=FALSE){
  N <- nrow(X)
  
  eigen_prob <- RSpectra::eigs_sym(S_cross, k=1)
  vec_1 <- eigen_prob$vectors
  lambda_1 <- eigen_prob$values
  res <- t(X) %*%vec_1 /sqrt((N-1)*lambda_1)
  if(return_lambda) {
    res <- list(vec=res, lambda=lambda_1)
  }
  res
}

#' simulate data
#' 
#' @param n sample size, corresponds to T in the T/N notation
#' @param lambda eigenvalues of thye cov matrix
mtr_mvnorm_sim <- function(n=100, lambda=NULL, mu=rep(0, length(lambdas)),
                           is_diagonal=FALSE,
                           seed = NULL) {
  p <- length(lambda)
  if(!is.null(seed)) set.seed(seed)
  
  
  ## from https://stats.stackexchange.com/questions/215497/how-to-create-an-arbitrary-covariance-matrix
  if(is_diagonal){
    Sig_half <- diag(sqrt(lambda))
    
  } else {
    Q <- clusterGeneration:::genOrthogonal(p)
    Sig_half <- Q*sqrt(lambda)
  }
  
  X <- matrix(rnorm(p * n), n)
  t(mu + Sig_half %*% t(X)) ## from MASS::mvrnorm
  
}

if(FALSE){
  p <- 200
  lambdas <- c(p, rep(1, p-1))
  lambdas_rel <- 100* lambdas/sum(lambdas)
  X_sim <- mtr_mvnorm_sim(n=30, lambdas)
  dim(X_sim)
  
  ## Get all R2 to the optimal correlation-index
  r2_opt <- mtr_r2_optim(X_sim)
  
  ## verify that their mean is the first eigenvalue, computed manually
  fully_manual <- function(df_w) {
    eigs_cor <- eigen(cor(df_w))$value 
    eigs_cor[1]/sum(eigs_cor)
  }
  fully_manual_cross <- function(df_w) {
    eigs_cor <- eigen(mtr_cor_cross(df_w))$value 
    eigs_cor[1]/sum(eigs_cor)
  }
  all.equal(mean(r2_opt),
            fully_manual(X_sim))
  
  ## Check speed
  lambdas_big <- c(p, rep(1, 1000-1))
  X_sim_big <- mtr_mvnorm_sim(n=30, lambdas_big)
  microbenchmark::microbenchmark(fo_cross = mean(mtr_r2_optim(X_sim_big, cross = TRUE)),
                                 fo_no_cross = mean(mtr_r2_optim(X_sim_big, cross = FALSE)),
                                 manual = fully_manual(X_sim_big),
                                 manual_cross = fully_manual_cross(X_sim_big),
                                 times = 20,check = "equal")
}