#' reads csv data file and the training list
#' @param data_name name of the csv file
#' @return list of the input data, the label matrix and the training list
#' 
read_data <- function(data_name) {
  
  dat <- read.csv(paste0("data/",data_name,".csv"))
  indices <- 1:nrow(dat)
  # dat_trainlist <- readRDS(paste0("data/",data_name,".trainlist.rds"))

  ### the dimensions of the input space is entered manually here
  dim_X <- list(fungi=9, emotions=72, scene=294, yeast=103)
  
  d <- as.numeric(dim_X[data_name])
  
  trainlist <- readRDS(paste0("data/",data_name,"_trainlist.rds"))
  ### we used a subsample of scene and yeast datasets
  if(data_name=="scene" || data_name=="yeast") {
      indices <-readRDS(paste0("data/",data_name,".subsample.rds"))
        }

  list(X=dat[indices,c(1:d)], Y=dat[indices,(1+d):ncol(dat)], trainlist=trainlist)

  }

  
  

