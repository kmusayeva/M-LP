
read_data <- function(data_name) {
  
  dat <- read.csv(paste0("data/",data_name,".csv"))
  
  ### the dimensions of the input space is entered manually here
  dim_X <- list(fungi=9, emotions=72, scene=294, yeast=103)
  
  d <- as.numeric(dim_X[data_name])
  
  ### we used a subsample of scene and yeast datasets
  if(data_name=="scene" | data_name=="yeast") {
      indices <-readRDS(paste0("data/",data_name,".subsample"))
      return(list(X=dat[indices,c(1:d)], Y=dat[indices,(1+d):ncol(dat)]))
      }
  
  return(list(X=dat[,c(1:d)], Y=dat[,(1+d):ncol(dat)]))
  
  }

  
  

