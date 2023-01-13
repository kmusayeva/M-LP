#################################################################
# copyright (C) 2023 K.Musayeva
# This file contains the implementation of the following 
# transition/propagation matrices:
# 1. Normalized Laplacian (Zhou et al., 2004)
# 2. Locally linear neighbourhood (Wang and Zang, 2008)
# 3. Gaussian transition matrix (Musayeva and Binois, 2023)
#################################################################


### computes the normalized Laplacian propagation matrix
propagation_matrix_normalized_Laplacian <- function(dist_matrix, sigma) {

  ### returns similarity matrix with diagonal elements equal to zero
  P <- weight_matrix(dist_matrix, sigma, diagonal=F)

  D <- diag((rowSums(P))^(-1/2))
  
  return(D%*%P%*%D)

    }


### computes the transition matrix based on Gaussian kernel
transition_matrix_gaussian <- function(dist_matrix, label, nn, sigma, row_norm=FALSE) {

  n <- nrow(dist_matrix)

  P <- as.matrix(weight_matrix(dist_matrix, sigma, diagonal=TRUE))

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


### implements the propagation matrix of Wang and Zang, 2008
transition_matrix_linear_patch <- function(X, X.dist_matrix.squared, nn) {

  dist.ranked <- t(apply(X.dist_matrix.squared, 1, rank, ties.method= "first"))

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

    settings <- osqpSettings(verbose = FALSE, eps_abs = 0.0001, eps_rel = 0.0001, eps_prim_inf = 1e-10, eps_dual_inf = 1e-10)
    # settings <- osqpSettings(verbose = FALSE)
    # Change alpha parameter and setup workspace
    model <- osqp(P, q, A, l, u, settings)

    # Solve problem
    res <- model$Solve()
    W[i,nns[i,]] <- res$x

    }

  return(W)

  }


### does chi-square transformation on 0/1 matrix
chi_square_transform <- function(Y) {

  total <- sum(Y)
  rows <- rowSums(Y)
  cols <- colSums(Y)

  step1 <- Y/ifelse(rows<1,1,rows)

  return(sqrt(total)*t(t(step1)/sqrt(ifelse(cols<1,1,cols))))

}



### computes a similarity matrix
weight_matrix <- function(dist_matrix, sigma, diagonal=F) {

  W <- exp(-dist_matrix/sigma^2)

  if(diagonal==F) diag(W) <- 0

  return(W)

}

### row normalizes the given matrix
row_normalize <- function(P) {
  return(sweep(P, MARGIN=1, FUN="/", STATS=ifelse(colSums(P)<1,1,rowSums(P))))
}



### column normalizes the given matrix
column_row_normalize <- function(P) {
  P.col.norm <- sweep(P, MARGIN=2, FUN="/", STATS=ifelse(colSums(P)<1,1,colSums(P)))
  P.row.norm <- sweep(P.col.norm, MARGIN=1, FUN="/", STATS=ifelse(rowSums(P.col.norm)<1e-15,1,rowSums(P.col.norm)))
  return(P.row.norm)
}
