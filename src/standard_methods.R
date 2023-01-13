
standard_methods_cross_validate <- function(X, Y, method, hyperparams, data_name, verbose=F) {

  if(method$name=="mlknn") { X <- as.data.frame(apply(X, 2, min_max_norm)) } ### mlknn is based on the Euclidean distance, thus we need to normalize

  df_tmp <- as.data.frame(cbind(X, Y))

  df_mldr <- mldr_from_dataframe(df_tmp, labelIndices=c((ncol(X)+1):ncol(df_tmp)))

  Ystr <- as.factor(apply(Y, 1, function(x) paste(x, collapse="")))

  trainlist <- createTrainList(Ystr, nfolds, nrcv, 456)

  hyperparams_grid <- as.data.frame(expand.grid(hyperparams))
  
  perf_list <- list()

  best_perf <- best_perf_init(names(hyperparams))

  for(i in 1:nrow(hyperparams_grid)) {
    
    hyperparam <- hyperparams_grid[i,,drop=FALSE]

    if(verbose) print(hyperparam)

    res <- cross_validate_mldr(df_mldr, Y, trainlist, method, hyperparam)

    if(is.null(res)) stop("An error occurred.")

    best_perf <- best_perf_update(res, best_perf, hyperparam)

    if(verbose) print(best_perf)

    A <- list(params=hyperparam, objectives=res)

    perf_list[[i]] <- A

    }


  file <- paste0(results_dir, data_name,"_nfolds_",nfolds,"_nrcv_",nrcv,"_", paste0(method,collapse="_"))

  suppressWarnings(write.table(best_perf, paste0(file, ".txt"), append=TRUE))

  saveRDS(perf_list, paste0(file, ".obj"))

  return(perf_list)

  }



cross_validate_mldr <- function(df_mldr, Y, trainlist, method, hyperparam) {

  l <- length(trainlist)

  perf_eval <- matrix(rep(0, l*19), ncol=19)

  for(k in 1:l) {

    label <- trainlist[[k]]

    train <- df_mldr[label]

    test <- df_mldr[-label]

    param <- as.list(hyperparam)

    param[["mdata"]] <- train

    model <- do.call(method$name, param)

    preds <- predict(model, test)

    preds_thresh <- do.call(method$df, list(Y[label,], preds))

    result <- multilabel_evaluate(test, preds_thresh, c("example-based", "label-based","ranking"), labels=TRUE)$multilabel

    perf_eval[k,]  <- result

  }

  colnames(perf_eval) <-  names(result)

  cv_perf <- data.frame(ave=round(colMeans(perf_eval),4),sd=round(apply(perf_eval, 2, sd),4))

  return(cv_perf)

}

