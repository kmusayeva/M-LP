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
library(Matrix)
library(GPGame)


# cross validation
nfolds <-5
nrcv <- 2

# the evaluation metrcis we are interested in
measures <- c("hamming-loss", "subset-accuracy", "F1", "macro-F1", "micro-F1", "macro-AUC", "micro-AUC", "average-precision")

# evaluation metrics for KS solution; any subset of the measures above
measures_to_optimize <- c("hamming-loss", "F1", "macro-F1", "micro-F1", "macro-AUC", "micro-AUC", "subset-accuracy", "average-precision")

# loss functions
loss <- c("hamming-loss","coverage","ranking-loss")

# thresholds for lco
thresholds <- seq(0.01, 0.9, by=0.01)

