#' @author K.Musayeva, M. Binois, 2023.
#' Computes the KS solution from the given list of performances evaluated
#' according to multiple evaluation measures over a hyperparameter grid.
#'
#' @param perf list containing the performance metrics evaluated over a hyperparameter grid
#' @param measures_to_optimize vector of evaluation measures to be optimized
#' @param verbose logical indicating whether to print the compromise solution.
#'
#' @return A list containing the hyperparameter configuration and the corresponding values of the evaluation measures.
ks_solution <- function(perf, measures_to_optimize, verbose=FALSE) {
  
  obj <- if(is.list(perf)) perf  else obj <- readRDS(file = perf)
  
  # get variables and objective matrices
  Xs <- Reduce(rbind, lapply(obj, function(alist) alist$params))
  Ys <- as.data.frame(Reduce(rbind,lapply(obj, function(alist) t(alist$objectives[measures_to_optimize,])[1,])))
  
  Ys[,measures_to_optimize[!(measures_to_optimize %in% loss)]] <- -Ys[,measures_to_optimize[!(measures_to_optimize %in% loss)]]

  ks_solution <- getEquilibrium(as.matrix(Ys), equilibrium = "KSE", nobj = length(measures_to_optimize), return.design = TRUE)

  params <- Xs[ks_solution$NE, ,drop = FALSE]
   
  ks_solution_values <- get_objective_values(obj, params)[measures,]

  # print the compromise solution
  if(verbose) {
    for(i in 1:nrow(ks_solution_values)) {
        cat(paste(round(ks_solution_values[i,1],3),"±",round(ks_solution_values[i,2],3)))
        cat("\n")
      }
    }

  list(params=params, values=ks_solution_values)

  }

#' Returns the values of the evaluation measures corresponding to the given hyperparameter configuration.
#'
#' @param obj list containing the performance metrics evaluated over a hyperparameter grid.
#' @param params hyperparameter configuration.
#'
#' @return matrix containing the values of the evaluation measures.
#'
get_objective_values <- function(obj, params) {

  for(i in 1:length(obj)) {
      if(isTRUE(all.equal(obj[[i]]$params,params))) return(obj[[i]]$objectives)
     }

  return(0)
  
  }


