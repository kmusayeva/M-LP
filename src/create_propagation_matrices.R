#' @author K.Musayeva 2023 <khmusayeva@gmail.com>
#'
#' Computes the Laplacian matrix
#' 
#' @param dist_matrix_squared squared euclidean distance matrix
#' @param sigma width of the Gaussian kernel
#' 
#' @return Laplacian matrix
#' 
Laplacian_matrix <- function(dist_matrix_squared, sigma) {
  
  W <- weight_matrix(dist_matrix_squared, sigma, diagonal=F)
  
  D <- diag((rowSums(W)))
  
  return(D-W)
  
}




#' Computes the normalized Laplacian propagation matrix
#' 
#' @param dist_matrix_squared squared euclidean distance matrix
#' @param sigma width of the Gaussian kernel
#' 
#' @return Normalized Laplacian propagation matrix
#' 
propagation_matrix_normalized_Laplacian <- function(dist_matrix_squared, sigma) {
  ### returns similarity matrix with diagonal elements equal to zero
  P <- weight_matrix(dist_matrix_squared, sigma, diagonal=F)
  
  D <- diag((rowSums(P))^(-1/2))
  
  return(D%*%P%*%D)
  
}



#' Computes the transition matrix based on Gaussian kernel
#' 
#' @param dist_matrix_squared squared euclidean distance matrix
#' @param label training indices
#' @param nn number of nearest neighbors
#' @param sigma width of the Gaussian kernel
#' @param row_norm whether to row-normalize the matrix (default: FALSE)
#' 
#' @return Transition matrix based on Gaussian kernel
#' 
transition_matrix_gaussian <- function(dist_matrix_squared, label, nn, sigma, row_norm=FALSE) {
  
  n <- nrow(dist_matrix_squared)
  
  P <- round(as.matrix(weight_matrix(dist_matrix_squared, sigma, diagonal=TRUE)),3)
  
  P[label,] <- 0
  
  diag(P) <- 1
  
  if(nn<n) {
    
    neighbours_matrix <- t(apply(P, 1, order, decreasing=TRUE))[,1:as.integer(nn)]
    
    for(i in rownames(neighbours_matrix)) {
      
      P[i,-neighbours_matrix[i,]] <- 0
      
    }
    
  }
  
  if(row_norm) return(row_normalize(P))
  
  return(column_row_normalize(P))
  
}



#' Implements the propagation matrix of Wang and Zang, 2008
#' 
#' @param X input data
#' @param X_dist_matrix_squared squared euclidean distance matrix
#' @param nn number of nearest neighbors
#' 
#' @return Locally linear neighborhood propagation matrix
#' 
transition_matrix_linear_patch <- function(X, X_dist_matrix_squared, nn) {
  
  dist.ranked <- t(apply(X_dist_matrix_squared, 1, rank, ties.method= "first"))
  
  nns <- (dist.ranked<=nn+1 & dist.ranked>1)
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  W <- matrix(0, nrow=n, ncol=n)
  
  q <- rep(0,nn)
  
  A <- rbind(rep(1,nn),diag(nrow=nn))
  
  l <- c(1., rep(0,nn))
  
  u <- rep(1., nn+1)
  
  for(i in 1:n) {
    # k <- sum(nns[i,])
    # X <- X.scaled ### CHANGE
    # calculate the differences  between xi and its neighbours
    Z <- matrix(c(t(X)) - c(t(X[i,])), nrow=nrow(X), byrow = TRUE)
    Z <- matrix(Z[nns[i,],], ncol=d, nrow=nn )
    
    #gram-matrix
    G <- Z%*%t(Z)
    
    # Define problem data
    P <- 2*G
    
    settings <- osqpSettings(verbose = FALSE, eps_abs = 1e-10, eps_rel = 1e-10, eps_prim_inf = 1e-10, eps_dual_inf = 1e-10)
    # settings <- osqpSettings(verbose = FALSE)
    # Change alpha parameter and setup workspace
    model <- osqp(P, q, A, l, u, settings)
    
    # Solve problem
    res <- model$Solve()
    W[i,nns[i,]] <- res$x
    
  }
  
  return(W)
  
}



#' k-nearest neighbour similarity matrix
#' @param dist_matrix_squared squared euclidean distance matrix
#' @param nn number of nearest neighbors
#' @param sigma width of the Gaussian kernel
similarity_matrix <- function(dist_matrix_squared, nn, sigma){
  
  P <- as.matrix(weight_matrix(dist_matrix_squared, sigma, diagonal=FALSE))
  
  neighbours_matrix <- t(apply(P, 1, order, decreasing=TRUE))[,2:(as.integer(nn)+1)]
  
  for(i in rownames(neighbours_matrix)) {
    
    P[i,-neighbours_matrix[i,]] <- 0
    
    }
  
  return(P)
  
  }


#' implements chi-square transformation on 0/1 matrix
#' @param Y the label matrix
chi_square_transform <- function(Y) {

  total <- sum(Y)
  rows <- rowSums(Y)
  cols <- colSums(Y)

  step1 <- Y/ifelse(rows<1,1,rows)

  return(sqrt(total)*t(t(step1)/sqrt(ifelse(cols<1,1,cols))))

}



#' computes a gaussian similarity matrix
#' @param dist_matrix_squared squared euclidean distance matrix
#' @param sigma width of the Gaussian kernel
#'@param diagonal by default the diagonal of the matrix is zero
weight_matrix <- function(dist_matrix_squared, sigma, diagonal=F) {

  W <- exp(-dist_matrix_squared/sigma^2)

  if(diagonal==F) diag(W) <- 0

  return(W)

}

#' row normalizes the given matrix
#' @param P propagation matrix
row_normalize <- function(P) {
  return(sweep(P, MARGIN=1, FUN="/", STATS=ifelse(rowSums(P)<1,1,rowSums(P))))
}



#' column normalizes the given matrix
#' @param P propagation matrix
column_row_normalize <- function(P) {
  P.col.norm <- sweep(P, MARGIN=2, FUN="/", STATS=ifelse(colSums(P)<1,1,colSums(P)))
  P.row.norm <- sweep(P.col.norm, MARGIN=1, FUN="/", STATS=ifelse(rowSums(P.col.norm)<1e-15,1,rowSums(P.col.norm)))
  return(P.row.norm)
}
