### this file contains some auxilary functions


avgImb <- function(Y) {
  n <- nrow(Y)
  num <- function(a) {return(max(a,n-a))}
  den <- function(a) {return(min(a,n-a))}

  Fj <- colSums(Y)
  return(mean(sapply(Fj,num)/sapply(Fj,den)))
}



min_max_norm <- function(x) { return((x-min(x)) / (max(x)-min(x))) }


is_nan_data <- function(x)  do.call(cbind, lapply(x, is.nan))


createTrainList <- function(Ystr,nfolds,nrcv,seed) {
  set.seed(seed)
  return(caret::createMultiFolds(Ystr, k=nfolds, times=nrcv))
}



best_perf_init <- function(params, eval_measures=NULL) {

  if(is.null(eval_measures)) eval_measures <- c("hamming-loss", "coverage","ranking-loss", "F1", "macro-F1","micro-F1","macro-AUC","micro-AUC", "subset-accuracy", "average-precision")

  best_perf <- matrix(0,ncol=(length(params)+2),nrow=length(eval_measures))

  colnames(best_perf) <- c("val", "sd", params)

  rownames(best_perf) <- eval_measures

  if("hamming-loss" %in% eval_measures) best_perf["hamming-loss",1] <- Inf
  if("coverage" %in% eval_measures) best_perf["coverage",1] <- Inf
  if("ranking-loss" %in% eval_measures) best_perf["ranking-loss",1] <- Inf

  return(best_perf)

 }



best_perf_update <- function(res, best_perf, params) {

  eval_measures <- rownames(best_perf)

  for(i in 1:length(eval_measures)) {

    b <- best_perf[i,1]

    current <- unlist(res[eval_measures[i],])

        # if(is.na(r)) r <- 0

    if(((eval_measures[i] %in% loss) & b>current[1]) |  (!(eval_measures[i] %in% loss) & b<current[1]))
      {
        # if((i>3 & b<r[1]) | (i<=3 & b>r[1])) {

      best_perf[i,1:2] <- current[1:2]

      best_perf[i, 3:(2+length(params))] <- unlist(params)

         }

     }

  return(best_perf)

  }


