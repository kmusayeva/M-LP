library(matrixcalc)
library(glmnet)

#' Exact-F1-plugin classifier (Dembczyński et al., 2013) 
#' Implemented K.Musayeva 2023 <khmusayeva@gmail.com>
#' 
#' @param X input data
#' @param Y label matrix
#' @param measures vector of evaluation measures
#' @param alpha alpha parameter in glmnet
#' @param lambda lambda parameter in glmnet
#' @return average and std of cross-validated performances wrt specified measures
#' 
#' @examples
#' alpha <- 0.1
#' lambda <- 0.05
#' exact_F1_plugin_evaluate(X, Y, trainlist, measures, alpha, lambda)
#' 

exact_F1_plugin_evaluate <- function(X, Y, trainlist, measures, alpha, lambda) {
  
  # Ystr <- as.factor(apply(Y, 1, function(x) paste(x, collapse="")))
  # trainlist <- createTrainList(Ystr, 5, 1, 456)
  
  Y <- as.matrix(Y)
  m <- ncol(Y)
  SY <- rowSums(Y)
  nsy <- max(SY)
  YS <- hadamard.prod(matrix(rep(SY,m),ncol=m), Y)+1
  Yzeroes <- ifelse(SY>0,0,1) #all zeroes
  
  l <- length(trainlist)  
  performance.eval <- matrix(rep(0, l*length(measures)), ncol=length(measures))
  
  for(k in 1:l) {
    
    label <- trainlist[[k]]
    X_train <- X[label,]
    Y.train <- Y[label,]
    X_test <- X[-label,]
    Y.test <- Y[-label,]
    
    YS_train <- YS[label,]
    
    Yzeroes.train <- Yzeroes[label]
    
    preds <- exact_F1_plugin(X_train,YS_train,Yzeroes.train,X_test, m, alpha, lambda) 
    
    Ytest.mldr <- mldr_from_dataframe(as.data.frame(Y.test), labelIndices=c(1:m))
    
    result <- multilabel_evaluate(Ytest.mldr, preds, measures, labels=TRUE)
    
    result <- result$multilabel
    
    performance.eval[k,] <- result
    
        }
  
  colnames(performance.eval) <- names(result)
  cv.performance <- data.frame(ave=colMeans(performance.eval),sd=apply(performance.eval, 2, sd))
  
  cv.performance <- cv.performance[measures,]
  cv.performance
  
  }




exact_F1_plugin <- function(X_train, YS_train, Yzeroes_train, X_test, m, alpha, lambda) {

  PY <- matrix(0,nrow=nrow(X_test),ncol=(m+1)*m)
  Pzeroes <- matrix(0,nrow=nrow(X_test),ncol=1)
  
  ### Step 1: compute the probability matrix by glm
  for(i in 1:m) {
    # check the frequency of classes for the i-th label
    YS.i.freq <- table(YS_train[,i])
    
    # is there a class that appears only once? If it is true then remove the corresponding rows.
    tmp <- which(YS.i.freq==1)
    
    
    if(length(tmp)>0) {
      
      row_ind <- which(YS_train[,i] %in% as.numeric(names(YS.i.freq[tmp])))
      
      X_train <- X_train[-row_ind,]
      
      YS_train_prime <- YS_train[-row_ind, i]
    
       }  

    else {

      YS_train_prime <- YS_train[,i]

      }
   
    # prepare labels 
    indices <- sort(unique(YS_train_prime))
   
    L <- sapply((YS_train_prime),function(x) which(indices==x)) ### for instance (1,3,4) ---> (1,2,3)

    print
    
    model <- glmnet(X_train, y=L, family="multinomial", lambda=lambda, alpha=alpha) 
    
    res <- predict(model, newx=as.matrix(X_test), type="response")

    start <- (m+1)*(i-1)
    
    PY[,c(start+indices)] <- res
    
   }
  
  
  ### The case of zero vector: binary logistic regression
  # model.zeroes <- glmnet(X_train, y=Yzeroes.train, family="multinomial", lambda=c(0.05), alpha=alpha2)
  # Pzeroes <- predict(model.zeroes, newx=as.matrix(X_test), type="response")
  # h0 <-Pzeroes[,,1][,2]

  ### Step 2: compute Delta matrix and formula (8)
  
  # s <- sort(unique(SY),decreasing = FALSE)
  # if(length(s)<m) s <- c(s,rep(0,m-length(s)))
  
  test.n <- nrow(X_test)
  
  hk <- matrix(0, nrow=test.n, ncol=m)
  
  sm <- matrix(rep(c(1:m),each=m),nrow=m,byrow=TRUE)
  
  W <- 1/(sapply(1:m, function(i) sm[,i]+i))
  
  
  for(i in 1:test.n) {
    
    P <- matrix(PY[i,], nrow=m, byrow = TRUE)[,2:(m+1)]
    
    delta <- P %*% W
    
    result.best <- 0
    
    for(k in 1:m) {
      
      ind <- order(delta[,k],decreasing = TRUE)[1:k]
      
      res <- sum(delta[ind,k])
      
      if(result.best < res) { result.best <- res; labels.best<-rep(0,m); labels.best[ind]<-1 }
      
      }
    
    hk[i,] <- labels.best
    
  }
  
  return(hk)
  
}




