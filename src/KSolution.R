##################################################################################
# copyright (C) 2023 K.Musayeva, M. Binois
# This file computes the KS solution from the given list of performances evaluated 
# according to multiple evaluation measures over a hyperparameter grid.
##################################################################################

library(GPGame)

### computes the compromise solution for multiple evaluation metrics
ks_solution <- function(perf, measures, verbose=FALSE) {
  
  obj <- if(is.list(perf)) perf  else obj <- readRDS(file = perf)
  
  # get variables and objective matrices
  Xs <- Reduce(rbind, lapply(obj, function(alist) alist$params))
  Ys <- as.data.frame(Reduce(rbind,lapply(obj, function(alist) t(alist$objectives[measures,])[1,])))
  
  Ys[,measures[!(measures %in% loss)]] <- -Ys[,measures[!(measures %in% loss)]]

  ks_solution <- getEquilibrium(as.matrix(Ys), equilibrium = "KSE", nobj = length(measures), return.design = TRUE)

  params <- Xs[ks_solution$NE, ,drop = FALSE]
  
  ks_solution_values <- get_objective_values(obj, params)[measures,]

  # print the compromise solution
  if(verbose) {
    for(i in 1:nrow(ks_solution_values)) {
        cat(paste(colnames(Ys)[i], ":", round(ks_solution_values[i,1],4),"±",round(ks_solution_values[i,2],4)))
        cat("\n")
      }
    }

  return(list(params=params, values=ks_solution_values))

  }


### returns the values of the evaluation measures corresponding to the given hyperparameter configuration
get_objective_values <- function(obj, params) {

  for(i in 1:length(obj)) {
      if(isTRUE(all.equal(obj[[i]]$params,params))) return(obj[[i]]$objectives)
     }

  return(0)
  
  }


