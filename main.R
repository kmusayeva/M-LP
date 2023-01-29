##################################################################################################
# copyright (C) 2023 K.Musayeva
# This file contains examples on how to run the label propagation 
# algorithms, the harmonic function as a stacking method, as well as the 
# inductive methods used in the paper
##################################################################################################

rm(list=ls()) 

### change to the right directory
root_dir <- "~/M-LP/"
source_file_dir <- "src/"
results_dir <- "results/"


### load all source files
setwd(dir=root_dir) 
file_sources = list.files(path=source_file_dir, pattern="*.R")
sapply(paste0(source_file_dir, file_sources), source, .GlobalEnv)

########################## DATA #############################

data_name <- "yeast" # can be {fungi, emotions, yeast, scene}

dat <- read_data(data_name)
X <- dat$X
Y <- dat$Y

###################### LABEL PROPAGATION ####################

# method <- list(name=name_of_method, df=decision_function) 
# available methods: {hf, tram, dlp, cm, lln}
# available decision functions: {cmn, lco}

# examples

method <- list(name="hf", df="lco") 

#method <- list(name="tram") #the TRAM and LLN methods has their own decision functions, no need to enter "df"

#method <- list(name="dlp", df="cmn")

#method <- list(name="cm", df="cmn")

#method <- list(name="lln", df="lco")

### hyperparameters (example for each method)

###hf/tram
hyperparams <- list(sigma=seq(0.5, 1, by=0.2), nn=seq(25, 65, 5))

###iterative hf
#hyperparams <- list(iter=seq(1, 20, by=2),sigma=seq(1.2, 2, by=0.2), nn=seq(35, 65, 5))

###lln
#hyperparams <- list(nn=seq(5, 5, 2), reg=seq(0.2, 0.2,0.2))

###cm
#hyperparams <- list(sigma=seq(0.1, 1, by=0.2), reg=seq(0.05,0.5,0.05))

###dlp
#hyperparams <- list(iter=seq(1,5,1), sigma=seq(1, 2, by=0.2), nn=seq(15, 50, 5), lambda=c(0.01), alpha=c(0.0001))


### cross validates the method and saves the performance of the method on the hyperparameter grid to a file
result <- mlp_cross_validate(X, Y, method, hyperparams, stacking=NULL, data_name, verbose=T)


###################### STACKING ######################

### Example: we want to use the outputs of BR in HF
### read the outputs of BR 
ks_file <- "results/scene_nfolds_2_nrcv_5_br_cmn.obj"

### find the compromise solution and the corresponding hyperparameters
ks <- ks_solution(ks_file, measures_to_optimize, verbose=TRUE)

reg <- 0.5 #regularization parameter to penalize the label matrix

# in the methods field, enter the method whose outputs you want to use in the stacking:
# this should correspond to what is specified in the ks_file
stacking <- list(method=list(name="br", df="cmn"), params=list(conf=ks$params, alpha=reg)) 

method <- list(name="hf", df="lco") 

hyperparams <- list(sigma=seq(1.2, 2, by=0.2), nn=seq(25, 60, 5))

result <- mlp_cross_validate(X, Y, method, hyperparams, stacking=stacking, data_name, verbose=T)


###################### INDUCTIVE METHODS ######################
# ecc-svm
#hyperparams <- list(base.algorithm = getOption("utiml.base.algorithm","SVM"), gamma=c(0.1,0.3,0.5,0.7,1), cost=c(0.5,0.7,1), subsample=c(0.5,0.7,1))

# rakel
# hyperparams <- list(base.algorithm = getOption("utiml.base.algorithm","SVM"), gamma=c(0.3, 0.5, 0.7, 1, 10), cost=c(0.5, 0.7, 1), k=c(3:ncol(Y)))

# mlknn
# hyperparams <- list(k=seq(5, 10, by=2))

# br-svm
hyperparams <- list(base.algorithm = getOption("utiml.base.algorithm","SVM"), gamma=c(0.01,0.1,0.3,0.5,0.7,1), cost=c(0.5,0.7,1))
method <- list(name = "br", df = "cmn")

standard_methods_cross_validate(X, Y, method, hyperparams, data_name, verbose=T)



