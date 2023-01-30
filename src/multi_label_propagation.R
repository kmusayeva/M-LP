#####################################################################################
# copyright (C) 2023 K.Musayeva
# This file contains the implementation of the following label
# propagation algorithms: 
# 1. Harmonic function and iterative label propagation (Musayeva and Binois, 2023)
# 2. Consistency method (Zhou et al., 2004)
# 2. Locally linear neighbourhood (Wang and Zang, 2008)
# 3. Dynamic label propagation (Wang et al. 2013)
# 4. TRAM (Kong et al. 2011)
# 5. Harmonic function as a stacking method (Musayeva and Binois, 2023)
#####################################################################################


### multi-label propagation (delegates the query to the corresponding method)
mlp <- function(method, X, Y, label, hyperparam, dist_matrix_squared) {

  if((method$name %in% c("hf", "cm", "lln", "dlp", "tram"))==F) stop("Available methods: hf, cm, lln, dlp, tram", call. = FALSE)

  ### lln and tram methods have their decision functions
   if(method$name  %in% c("lln", "tram")) {
    return(do.call(method$name, list(Y, label, dist_matrix_squared, hyperparam)))
    }

  ### for the rest we can use one of the following decision functions
  if((method$df %in% c("lco","cmn"))==F) stop("Available decision functions: lco, cmn", call. = FALSE)

  preds <- do.call(method$name, list(Y, label, dist_matrix_squared, hyperparam))

  return(do.call(method$df, list(Y[label,], preds)))

  }



### Harmonic function and iterative label propagation
hf <- function(Y, label, dist_matrix, hyperparam) {

  if(is.null(hyperparam$nn)) stop("Number of neighbours is missing.")
  if(is.null(hyperparam$sigma)) stop("Please specify sigma.")

  transition_matrix <- transition_matrix_gaussian(dist_matrix, label, hyperparam$nn, hyperparam$sigma, row_norm=F)

  P_UU <- transition_matrix[-label,-label]
  P_UL <- transition_matrix[-label,label]

  
  u <- nrow(dist_matrix)-length(label)
  I <- diag(u)

  if(is.null(hyperparam$iter)) {
      preds <- ginv(I-P_UU) %*% P_UL %*% as.matrix(Y[label,])
    }

  ### if the iter parameter is not null, does iterative label propagation
  else {
    preds <- matrix(0,ncol=ncol(Y),nrow=u)

    for(i in 1:hyperparam$iter) {
      preds <- P_UL %*% as.matrix(Y[label,]) + P_UU %*% preds
        }

      }

    return(preds)

}


### Locally linear neighbourhood (Wang and Zang, 2008)
lln <- function(X, Y, label, dist_matrix, hyperparam) {

  if(is.null(hyperparam$nn)) stop("Number of neighbours is missing.")
  if(is.null(hyperparam$reg)) stop("Regularization parameter is missing.")

  n <- nrow(X)

  transition_matrix <- transition_matrix_linear_patch(X, dist_matrix, hyperparam$nn)

  Yprime <- as.matrix(Y)

  Yprime[which(Yprime==0)] <- -1

  Yprime[-label,] <- 0

  I <- diag(n)

  preds <- ginv(I-hyperparam$reg*transition_matrix) %*% Yprime

  preds <- preds[-label,]

  preds[which(preds>0)] <- 1

  preds[which(preds<=0)] <- 0

  return(preds)

  }




### Consistency method of Zhou et al., 2004
cm <- function(Y, label, dist_matrix, hyperparam) {

  if(is.null(hyperparam$reg)) stop("Regularization parameter is missing.")
  if(is.null(hyperparam$sigma)) stop("Please specify sigma.")

  transition_matrix <- propagation_matrix_normalized_Laplacian(dist_matrix, hyperparam$sigma)
  
  if(sum(is_nan_data(transition_matrix))) stop("Please choose different sigma.")
  
  Yprime <- as.matrix(Y)

  Yprime[-label,] <- 0

  I <- diag(nrow(Y))
  
  preds <- ginv(I-hyperparam$reg*transition_matrix) %*% Yprime

  return(preds[-label,])

}



### Dynamic label propagation of Wang et al. 2013
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

  return(Yprime[-label,])

  }


### TRAM (Kong et al. 2011): it is the harmonic function applied to the row-normalized transition matrix.
### The decision function used is also based on the harmonic function.
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


  return(lhf(preds))

}



### Harmonic function as a stacking method
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

  return(do.call(method$df, list(Y[label,], preds)))

}


