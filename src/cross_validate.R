####################################################
# copyright (c) 2023 K.Musayeva
####################################################

### Creates the Euclidean distance matrix of the input data.
### Creates a training list.
### Cross validates.
### Saves the performances over the hyperparameter grid (as well as the best performance) to a file 
### (this is then used to find the KS solution).

mlp_cross_validate <- function(X, Y, method, hyperparams, stacking=NULL, what, verbose=F) {
  ### creates an Euclidean distance matrix of the input data
  X_scaled <- as.data.frame(apply(X, 2, min_max_norm)) ## normalize all features to [0,1]
  X_dist_matrix <- as.matrix(stats::dist(X_scaled, method="euclidean"))
  X_dist_matrix_squared <- (X_dist_matrix)^2

  ### create a training set
  Ystr <- as.factor(apply(Y, 1, function(x) paste(x, collapse="")))
  trainlist <- createTrainList(Ystr, nfolds, nrcv, 456)

  ### creates a data frame containing all combinations of hyper-parameters
  hyperparams_grid <- as.data.frame(expand.grid(hyperparams))

  ### a list for book-keeping hyper parameters and corresponding performance as measured by multiple evalution metrics
  perf_list <- list()

  best_perf <- best_perf_init(names(hyperparams))

  df_tmp <- as.data.frame(cbind(X, Y))
  df_mldr <- mldr_from_dataframe(df_tmp, labelIndices=c((ncol(X)+1):ncol(df_tmp)))


  for(i in 1:nrow(hyperparams_grid)) {

    hyperparam <- if(ncol(hyperparams_grid)>1) hyperparams_grid[i,] else hyperparam <- hyperparams_grid[i]

    if(verbose) print(hyperparam)

    res <- cross_validate(mlp, df_mldr, trainlist, X_dist_matrix_squared, method, hyperparam, stacking)

    if(is.null(res)) stop("An error occurred.")

    best_perf <- best_perf_update(res, best_perf, hyperparam)

    if(verbose) print(best_perf)

    A <- list(params=hyperparam, objectives=res)

    perf_list[[i]] <- A

  }

  #### saves the results
  file <- paste0(results_dir, what,"_nfolds_",nfolds,"_nrcv_",nrcv,"_", paste0(method,collapse="_"))

  if(!is.null(stacking)) file <- paste0(file, "_",paste0(stacking$method,collapse="_"))

  suppressWarnings(write.table(best_perf, paste0(file, ".txt"), append=TRUE))

  saveRDS(perf_list, paste0(file, ".obj"))

  return(perf_list)

  }


### Implements cross validation
cross_validate <- function(f, df_mldr, trainlist, dist_matrix_squared, method, hyperparam, stacking) {

  d <- length(df_mldr$attributesIndexes)
  X <- df_mldr$dataset[, 1:d]
  Y <- df_mldr$dataset[, df_mldr$labels$index]

  l <- length(trainlist)

  performance_eval <- matrix(rep(0, l*19), ncol=19)

  for(k in 1:l) {

    label <- trainlist[[k]]


    if(!is.null(stacking)) {

      if(stacking$method$name == "hf") {
        
        preds <- hf(Y, label, dist_matrix_squared, as.list(stacking$params$conf))
        
        }

      else {

        params <- as.list(stacking$params$conf)
        params[["mdata"]] <- df_mldr[label]
        model <- do.call(stacking$method$name, params)
        preds <- predict(model, df_mldr[-label])

        }

      preds_thresh <- do.call(stacking$method$df, list(Y[label,], preds))

      preds <- do.call(paste0(as.character(substitute(f)),"_stacking"),list(Y, preds_thresh, label, stacking$params$alpha, method, hyperparam, dist_matrix_squared))

      }

    else preds <- f(method, X, Y, label, hyperparam, dist_matrix_squared)

    result <- multilabel_evaluate(df_mldr[-label], preds, c("example-based", "label-based","ranking"), labels=TRUE)

    result <- result$multilabel

    performance_eval[k,] <- result

    }


  colnames(performance_eval) <- names(result)

  cv_performance <- data.frame(ave=round(colMeans(performance_eval),4),sd=round(apply(performance_eval, 2, sd),4))

  return(cv_performance)

}
