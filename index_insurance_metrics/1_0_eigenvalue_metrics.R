#' Equation 5
#' 
#' 
#' @examples 
#' X_sim <- idx_mvnorm_sim(n=30, lambda = c(200, rep(1, 199)), seed = 123)
#' 
#' # R2 where index is the mean
#' idx_r2_general(X_sim)
#' 
#' ## compare with manual R2 with `lm()`:
#' all.equal(idx_r2_general(X_sim),
#'           sapply(summary(lm(X_sim~rowMeans(X_sim))), \(x) x$r.squared),
#'           check.attributes=FALSE)
#'           
#'  ## R2 with optimal index
#'  w_opt <- eigen(cor(X_sim))$vectors[,1] * sqrt(1/diag(cov(X_sim))) # Theorem 1, 1.1
#'  idx_r2_general(X_sim, w = w_opt)
#'  
#'  ## compare with eigenvalue
#'  all.equal(mean(idx_r2_general(X_sim, w = w_opt)), 
#'            eigen(cor(X_sim))$values |> (\(x) x[1]/sum(x))())
#'            mean(idx_r2_optimal(X_sim))
idx_r2_general <- function(df_w, w=matrix(1, ncol=1, nrow=ncol(df_w))){
  S <- cov(df_w)
  if(is.matrix(w) && ncol(w)>1) stop("w should be a vector of length ncol(df_w)")
  
  Sw <- S%*% w
  sum_S <- drop(t(w)%*%S%*% w)
  S_pred <- tcrossprod(Sw)
  diag(S_pred)/(diag(S)*sum_S)
}

#' Eigenvalue of the correlation matrix
#' 
#' 
#' @param df_w dataset in a wide format: T rows and N columns, for each field
#' @param to_df return data as data.frame?
#' @param cross Whether to use the cross trick?
idx_r2_optimal <- function(df_w, to_df=FALSE, cross=FALSE){
  

  if(cross){
    X_scale_col <- scale(df_w, scale = TRUE)
    # idx_fo <-  switch(type, "cov" = idx_cov_cross, "cor" = idx_cor_cross)
    S_cross <- idx_cor_cross(X_scale_col, is_scaled=TRUE)
    out <- idx_eigenvec_cross(X_scale_col, S_cross, return_lambda = TRUE)
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

#' Total R2: this is from Jong and Kotz, equ (5), remove n-1 and divide by total trace
#' @param S Covariance matrix
#' @param w weight vector
idx_total_r2_general <- function(df_w, w=matrix(1, nrow = nrow(S))){
  
  S <- cov(df_w)
  # Paper notation: sum(diag(S %*%t(w) %*%solve(w %*%S %*%t(w)) %*% w %*% S ))
  
  ## Change notation (transpose w)
  # Slowest: SSE <- sum(diag(S %*%w %*%solve(t(w) %*%S %*%w) %*% t(w) %*% S )) # this is lambda_1!
  # Slower: SSE <- solve(t(w) %*%S %*%w) %*% t(w) %*% S %*% S %*%w  # this is lambda_1!
  Sw  <- S %*%w
  SSE <- drop(solve(t(w) %*%Sw) %*% crossprod(Sw))
  SST <- sum(diag(S))
  SSE/SST
}

idx_total_r2_optimal <- function(df_w, cross=TRUE){
  
  S <- if(cross) idx_cov_cross(df_w) else cov(df_w)
  eig1 <- RSpectra::eigs_sym(S, k=1, opts = list(retvec=FALSE))
  eig1$values/sum(diag(S))
  
}

################################
#'## Utility functions
################################

## When we want only first eigens, use XX'/N (T x T) instead of X'X/N (N x N)!
## XX' in N x p notation, i.e. X'X in Ishi Yata et al 2014
idx_cov_cross <- function(X, is_scaled=FALSE){
  n <- nrow(X)
  if(!is_scaled) X <- scale(X, scale=FALSE)
  tcrossprod(X)/(n-1)
}

## same for cor!
idx_cor_cross <- function(X, is_scaled=FALSE){
  n <- nrow(X)
  if(!is_scaled) X <- scale(X, scale=TRUE)
  tcrossprod(X)/(n-1)
}

## first eigenvector cross, from Wang and Fan 2017, (2.1) p. 1348
idx_eigenvec_cross <- function(X, S_cross = idx_cov_cross(X), return_lambda=FALSE){
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

if(FALSE){
  p <- 200
  set.seed(1234)
  X_sim <- idx_mvnorm_sim(n=30, lambda = c(p, rep(1, p-1)))
  eig_S_fo <- idx_eigenvec_cross(X_sim)
  eig_R_fo <- idx_eigenvec_cross(scale(X_sim), S_cross = idx_cor_cross(X_sim))

  ## Not always TRUE due to normalization!  
  all.equal(eig_S_fo, eigen(cov(X_sim))$vectors[,1, drop=FALSE])
  all.equal(eig_R_fo, eigen(cor(X_sim))$vectors[,1, drop=FALSE])
}

#' simulate data
#' 
#' @param n sample size, corresponds to T in the T/N notation
#' @param lambda eigenvalues of thye cov matrix
idx_mvnorm_sim <- function(n=100, lambda=NULL, mu=rep(0, length(lambda)),
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
  100* lambdas[1]/sum(lambdas)
  X_sim <- idx_mvnorm_sim(n=30, lambdas)
  dim(X_sim)
  
  ## Get all R2 to the optimal correlation-index
  r2_R_opt <- idx_r2_optimal(X_sim)
  r2_S_opt <- idx_r2_optimal(X_sim, type = "cov")
  
  
  fully_manual <- function(df_w, type = c("cor", "cov"), cross=TRUE) {
    type <- match.arg(type)
    if(cross) {
      covr_fo <- switch(type, "cov" = idx_cov_cross, "cor" = idx_cor_cross)
    } else {
      covr_fo <- switch(type, "cov" = cov, "cor" = cor)
    }
    eigs_covr <- eigen(covr_fo(df_w))$value 
    eigs_covr[1]/sum(eigs_covr)
  }
  
  ## verify that their mean is the first eigenvalue, computed manually
  all.equal(mean(r2_R_opt),
            fully_manual(X_sim))
  v <- apply(X_sim, 2, var)
  w <- v/sum(v)
  idx_total_r2_general(df_w = X_sim)
  all.equal(idx_total_r2_general_optim(X_sim),
            fully_manual(X_sim, type = "cov", cross=FALSE))
  mean(idx_r2_optimal(X_sim, type = "cov", cross=FALSE))
  
  
  ## Check speed
  lambdas_big <- c(p, rep(1, 1000-1))
  X_sim_big <- idx_mvnorm_sim(n=30, lambdas_big)
  microbenchmark::microbenchmark(fo_cross = mean(idx_r2_optimal(X_sim_big, cross = TRUE)),
                                 fo_no_cross = mean(idx_r2_optimal(X_sim_big, cross = FALSE)),
                                 manual_cross = fully_manual(X_sim_big),
                                 manual_no_cross = fully_manual(X_sim_big, cross = FALSE),
                                 times = 20,
                                 check = "equal")
}