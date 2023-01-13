########################
library(dplyr)
library(tidyverse)
library(tidyr)
library(mldr)
library(stats)
library(caret)
library(MASS)
library(utiml)
library(vegan)
library(osqp)
########################

### global variables

# cross validation
nfolds <- 2
nrcv <- 5

# evaluation metrics for KS solution
measures_to_optimize <- c("hamming-loss", "F1", "macro-F1", "micro-F1", "macro-AUC", "micro-AUC", "accuracy", "subset-accuracy")

# loss functions
loss <- c("hamming-loss","coverage","ranking-loss")
