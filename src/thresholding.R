#' @author K.Musayeva 2023 <khmusayeva@gmail.com>
#' 
#' Label Propagation Algorithms/Instance-based logistic regression
#' 
#' This file contains the implementation of the following thresholding strategies:
#' 
#' \code{lco} Thresholding based on label cardinality (Read et al., 2011)
#' 
#' \code{cmn} Class-mass normalization (Zhu et al., 2003)
#' 
#' \code{lhf} Thresholding based on the harmonic function (Kong et al., 2011)
#' 
#' \code{basic} Thresholding at 0.5
#' 
#' \code{atzero} Thresholding at 0 
#'
#'
#' Implements the thresholding strategy of Read et al., 2011: we refer to it as the label cardinality optimizer (lco)
#'
#' @param preds raw predictions
#' @param Y_train training label matrix
#' @return thresholded predictions
lco <- function(preds, Y_train) {
  
  #index of the threshold vector that minimizes the absolute difference between the predicted and true label cardinalities
  which_threshold <- function(preds, Y_train) {
    
    lcard <- sum(Y_train)/(nrow(Y_train))
    l <- length(thresholds)
    res <- rep(0,l)
    
    for(i in 1:l) {
      preds_lcard <- sum(preds>=thresholds[i])/(nrow(preds))
      res[i] <- abs(lcard - preds_lcard)
    }
    
     which.min(res)
    
    }
  
  
  v <- which_threshold(preds, Y_train)
  
  fixed_threshold(preds, threshold = thresholds[v], probability = FALSE)

  }



#' Implements class-mass normalization method of Zhu et al., 2003
#'
#' @param preds raw predictions
#' @param Y_train training label matrix
#' @return thresholded predictions (in the mldr format)
cmn <- function(preds, Y_train) {

  Y_train.col <- ncol(Y_train)
  preds_n <- nrow(preds)

  A <- matrix(0, nrow=preds_n, ncol=Y_train.col)

  for(j in 1:Y_train.col) {
    Y_train.p <- mean(Y_train[,j])
    classj <- preds[,j]
    fu <- sum(preds[,j])
    fuc <- preds_n-fu
    A[,j] <- ifelse(sign(((Y_train.p*classj)/fu)-(1-Y_train.p)*((1-classj)/fuc)) < 0, 0, 1)
    A[is.na(A)] <- 0
  }

  convert_to_mldr_format(A, preds)

  }



#' The number of labels is computed based on the harmonic function (Kong et al. 2011).
#' We do not row-normalize the label matrix, thus we do not need the second parameter; it is for testing purposes.
#' @param preds raw predictions
#' @param preds_num_for_labels the number of relevant labels
#' @return thresholded predictions (in the mldr format)
lhf <- function(preds, ...) {

  A <- matrix(0, nrow=nrow(preds), ncol=ncol(preds))

  # if(preds_num_for_labels) {
  #   number_rel <- matrix(floor(rowSums(preds_num_for_labels)))
  # }
  # 
  # else {
    number_rel <- matrix(floor(rowSums(preds)))
        # }

  number_rel[which(number_rel<=0)] <- 1

  for(i in 1:length(number_rel)) {
    indices <- order(preds[i,],decreasing=T)[1:number_rel[i]]
    A[i,indices] <- 1
    }

  convert_to_mldr_format(A, preds)

}


#' standard thresholding at 0.5
#' @param  preds raw predictions 
basic <- function(preds, ...) {
  
  fixed_threshold(preds, threshold = 0.5, probability = FALSE)
  
  }


#' thresholding at 0
#' @param  preds raw predictions 
atzero <- function(preds, ...) {
  
  fixed_threshold(preds, threshold = 0, probability = FALSE)
  
}


#' convert the predictions into the mldr data format to use utiml evaluate method.
#' @param A thresholded predictions
#' @param preds raw predictions
convert_to_mldr_format <- function(A, preds){

  dimnames(A) <- dimnames(preds)

  ##if the thresholding strategy yields all zero labels for an instance,
  ## then convert to 1 the label with the highest score
  label <- apply(preds, 1, which.max)
  instance <- rowSums(A) == 0
  A[cbind(which(instance), label[instance])] <- 1


  ### this part is taken from mldr package
  preds <- as.matrix(preds)
  only.bipartitions <- A
  only.probabilities <- preds
  bipartitions <- A
  attr(preds, "classes") <- only.bipartitions
  attr(preds, "type") <- "probability"
  attr(bipartitions, "probs") <- only.probabilities
  attr(bipartitions, "type") <- "bipartition"
  class(preds) <- class(bipartitions) <- "mlresult"
  
  list(preds, bipartitions)[c(F, T)][[1]]

}



