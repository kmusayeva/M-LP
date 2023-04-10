#' @author K.Musayeva 2023 <khmusayeva@gmail.com>
#' 
#' Label Propagation Algorithms/Instance-based logistic regression
#' 
#' This file contains the implementation of the following label
#' propagation algorithms:
#' 
#' 
#' \code{hf} the Harmonic Function and Iterative Label Propagation (Musayeva and Binois, 2023)
#' 
#' \code{lln} the Locally Linear Neighborhood method of Wang and Zang (2008).
#' 
#' \code{cm} the Consistency Method of Zhou et al. (2004).
#' 
#' \code{dlp} the Dynamic Label Propagation method of Wang et al. (2013).
#' 
#' \code{tram} the TRAM method of Kong et al. (2011).
#' 
#' \code{mlp_stacking} the harmonic function as a stacking method (Musayeva and Binois, 2023).
#' 
#' \code{iblr} the Instance-Based Logistic Regression method (Cheng and Hullermeier, 2009).
#' 
#' 
#' 
#' 
#' Multi-label propagation (delegates the query to the corresponding method)
#'
#' @param method a list containing the name of the propagation method (hf, cm, lln, dlp, tram) and possibly other parameters
#' @param X input data
#' @param Y label matrix
#' @param label list of training indices
#' @param hyperparam list containing the hyperparameters
#' @param dist_matrix_squared squared euclidean distance matrix of scaled X
#' 
#' @return predictions for unlabeled points
#'
#' @export
mlp <- function(method, X, Y, label, hyperparam, dist_matrix_squared) {
  # if((method$name %in% c("hf", "cm", "lln", "dlp", "tram")) == FALSE) stop("Available methods: hf, cm, lln, dlp, tram", call. = FALSE)
  
  if (method$name == "tram") {
    method$df <- "lhf"
  } else if (is.null(method$df) || (method$df %in% c("lco", "cmn")) == FALSE) {
    method$df <- "basic"
  }
  
  if (method$name == "lln") {
    method$df <- "atzero"
    preds <- lln(X, Y, label, dist_matrix_squared, hyperparam)
  } 
  else {
    preds <- do.call(method$name, list(Y, label, dist_matrix_squared, hyperparam))
  }
  
  do.call(method$df, list(preds, Y[label,]))
}



#' Harmonic function and iterative label propagation
#'
#' @param Y label matrix
#' @param label list of training indices
#' @param hyperparam list containing the hyperparameters
#' @param dist_matrix_squared  squared euclidean distance matrix of scaled X
#' 
#' @return raw predictions for unlabeled points (before applying a thresholding function)
#'
hf <- function(Y, label, dist_matrix_squared, hyperparam) {
  
  if (is.null(hyperparam$nn)) stop("Number of neighbours is missing.")
  if (is.null(hyperparam$sigma)) stop("Please specify sigma.")
  
  transition_matrix <- transition_matrix_gaussian(dist_matrix_squared, label, hyperparam$nn, hyperparam$sigma, row_norm=F)
  
  P_UU <- transition_matrix[-label,-label]
  P_UL <- transition_matrix[-label,label]
  
  u <- nrow(dist_matrix_squared) - length(label)
  I <- diag(u)
  
  if (is.null(hyperparam$iter)) {
    preds <- ginv(I - P_UU) %*% P_UL %*% as.matrix(Y[label,])
  } else {
    preds <- matrix(0, ncol = ncol(Y), nrow = u)
    
    for (i in 1:hyperparam$iter) {
      preds <- P_UL %*% as.matrix(Y[label,]) + P_UU %*% preds
    }
  }
  
  preds
  
  }


#' Locally linear neighbourhood (Wang and Zang, 2008)
#' 
#' @param X input data
#' @param Y label matrix
#' @param label list of training indices
#' @param dist_matrix_squared squared euclidean distance matrix of scaled X
#' @param hyperparam list containing the hyperparameters
#'
#'  @return raw predictions for unlabeled points
#'  
lln <- function(X, Y, label, dist_matrix_squared, hyperparam) {
  
  if(is.null(hyperparam$nn)) stop("Number of neighbours is missing.")
  if(is.null(hyperparam$reg)) stop("Regularization parameter is missing.")
  
  n <- nrow(X)
  
  transition_matrix <- transition_matrix_linear_patch(X, dist_matrix_squared, hyperparam$nn)
  
  Yprime <- as.matrix(Y)
  
  Yprime[which(Yprime==0)] <- -1
  
  Yprime[-label,] <- 0
  
  I <- diag(n)
  
  preds <- ginv(I-hyperparam$reg*transition_matrix) %*% Yprime
  
  preds[-label,]

}



#' Consistency method of Zhou et al., 2004
#' 
#' @param Y label matrix
#' @param label list of training indices
#' @param dist_matrix_squared squared euclidean distance matrix of scaled X
#' @param hyperparam list containing the hyperparameters
#'
#' @return raw predictions for unlabeled points
#' 
cm <- function(Y, label, dist_matrix_squared, hyperparam) {

  if(is.null(hyperparam$reg)) stop("Regularization parameter is missing.")
  if(is.null(hyperparam$sigma)) stop("Please specify sigma.")

  transition_matrix <- propagation_matrix_normalized_Laplacian(dist_matrix_squared, hyperparam$sigma)

  if(sum(is_nan_data(transition_matrix))) stop("Please choose different sigma.")
  
  Yprime <- as.matrix(Y)

  Yprime[-label,] <- 0

  I <- diag(nrow(Y))

  preds <- ginv(I-hyperparam$reg*transition_matrix) %*% Yprime

  preds[-label,]

}



#' Dynamic label propagation of Wang et al. 2013
#' 
#' @param Y label matrix
#' @param label list of training indices
#' @param dist_matrix_squared squared euclidean distance matrix of scaled X
#' @param hyperparam list containing the hyperparameters
#'
#' @return  raw predictions for unlabeled points
#' 
dlp <- function(Y, label, dist_matrix, hyperparam) {

  if(is.null(hyperparam$nn)) stop("Number of neighbours is missing.")
  if(is.null(hyperparam$sigma)) stop("Please specify sigma.")
  if(is.null(hyperparam$iter)) stop("Please specify the number of iterations.")

  if(is.null(hyperparam$lambda)) {hyperparam$lambda <- 0.01}
  if(is.null(hyperparam$alpha)) {hyperparam$alpha <- 0.0001}


  transition_matrix <- transition_matrix_gaussian(dist_matrix, label, hyperparam$nn, hyperparam$sigma)

  full_transition_matrix <- transition_matrix_gaussian(dist_matrix, label, nrow(Y), hyperparam$sigma)

  Yprime <- as.matrix(Y)
  Yprime[-label,] <- 0

  I <- diag(nrow(Y))

  for(i in 1:hyperparam$iter) {

    Yprime <- full_transition_matrix %*% Yprime
    Yprime[label,] <- as.matrix(Y[label,])
    full_transition_matrix <- transition_matrix %*% (full_transition_matrix + hyperparam$alpha * Yprime %*% t(Yprime)) %*% t(transition_matrix) + hyperparam$lambda * I

    }

  Yprime[-label,]

  }



#' TRAM (Kong et al. 2011): it is the harmonic function applied to the row-normalized transition matrix.
#' The decision function used is also based on the harmonic function.
#' @param Y label matrix
#' @param label list of training indices
#' @param dist_matrix_squared squared euclidean distance matrix of scaled X
#' @param hyperparam list containing the hyperparameters
#'
#' @return predictions for unlabeled points
tram <- function(Y, label, dist_matrix, hyperparam) {

  if(is.null(hyperparam$nn)) stop("Number of neighbours is missing.")
  if(is.null(hyperparam$sigma)) stop("Please specify sigma.")

  transition_matrix <- transition_matrix_gaussian(dist_matrix, label, hyperparam$nn, hyperparam$sigma, row_norm=T)

  P_UU <- transition_matrix[-label,-label]

  P_UL <- transition_matrix[-label,label]

  u <- nrow(Y)-length(label)

  I <- diag(u)

  # preds <- ginv(I-P_UU) %*% P_UL %*% as.matrix(row_normalize(Y[label,]))

  preds <- ginv(I-P_UU) %*% P_UL %*% as.matrix(Y[label,])

  preds

}



#' Stacking for harmonic function and IBLR
#' @param Y label matrix
#' @param Y_preds predictions of the external method
#' @param label list of training indices
#' @param alpha regularization parameter for Y_preds
#' @param hyperparam list containing the hyperparameters
#' @param X_dist_matrix_squared squared euclidean distance matrix of scaled X
#'
#' @return predictions for unlabeled points
run_stacking <- function(Y, Y_preds, label, alpha, method, hyperparam, X_dist_matrix_squared) {
  
  if(method$name=="hf") {
    result <- mlp_stacking(Y, Y_preds, label, alpha, method, hyperparam, X_dist_matrix_squared)
    }
  
  else {
    preds <- iblr(Y, label, X_dist_matrix_squared, hyperparam, Y_preds)
    result <- do.call(method$df, list(preds, Y[label,]))
    }

    return(result)
}


#' Harmonic function as a stacking method
#' @param Y label matrix
#' @param Y_preds predictions of the external method
#' @param label list of training indices
#' @param alpha regularization parameter for Y_preds
#' @param hyperparam list containing the hyperparameters
#' @param X_dist_matrix_squared squared euclidean distance matrix of scaled X
#'
#' @return predictions for unlabeled points
mlp_stacking <- function(Y, Y_preds, label, alpha, method, hyperparam, X_dist_matrix_squared) {

  if(is.null(alpha)) stop("Please specify the regularization parameter.")

  nn <- hyperparam$nn

  sigma <- hyperparam$sigma

  iter <- hyperparam$iter

  Yprime <- as.matrix(Y)

  Yprime[-label,] <- as.matrix(Y_preds)

  u <- nrow(Y)-length(label)

  I <- diag(u)

  Yprime_chi <- chi_square_transform(Yprime)

  Y_dist_matrix <- (as.matrix(stats::dist(Yprime_chi, method="euclidean")))^2

  dist_matrix <- alpha * Y_dist_matrix + X_dist_matrix_squared

  preds <- hf(Y, label, dist_matrix, hyperparam)

  return(do.call(method$df, list(preds, Y[label,])))

}



#' Istance-based logistic regression  (Cheng and Hullermeier, 2009)
#' @param Y label matrix
#' @param Y_preds predictions of the external method
#' @param label list of training indices
#' @param hyperparam list containing the hyperparameters
#' @param X_dist_matrix_squared squared euclidean distance matrix of scaled X
#'
#' @return predictions for unlabeled points
iblr <- function(Y, label, X_dist_matrix_squared, hyperparam, Y_preds=NULL) {
  
  if(is.null(hyperparam$nn)) stop("Number of neighbours is missing.")
  
  if(is.null(hyperparam$sigma)) stop("Please specify sigma.")
  
  P <- similarity_matrix(X_dist_matrix_squared, hyperparam$nn, hyperparam$sigma)
  
  n <- nrow(Y)
  
  m <- ncol(Y)
  
  unlabel <- c(1:n)[-label]
  
  Y_prime <- Y
  
  if(!is.null(Y_preds)) {Y_prime[-label,] <- Y_preds; label <- c(1:n)}
  
  p <- colSums(Y_prime[label,])/length(label)
  
  omega_zero <- log(p)-log(1-p)
  
  omega_plus <- matrix(0, ncol=m, nrow=n)
  
  for(j in 1:m) {
    y <- ifelse(Y[label,j]==1,1,-1)
    omega_plus[,j] <- rowSums(sweep(P[,label], 2, y, `*`))
  }
  
  preds <- matrix(0, ncol=m, nrow=length(unlabel))
  
  for(j in 1:m) {
    
    omega_zero_rep <- rep(omega_zero[j], length(label))
    
    mod <- glm(Y_prime[label,j]~omega_plus[label,]+0+offset(omega_zero_rep),family=binomial)
    
    a <- exp(omega_zero[j]+rowSums(sweep(omega_plus[unlabel, ], 2, mod$coefficients, `*`)))
    
    pi0 <- a/(1+a)
    
    pi0[is.na(pi0)] <- 0
    
    preds[,j] <- pi0
    
  }
  
  # print(sum(is.na(preds)))
  
  return(preds)
  
}


