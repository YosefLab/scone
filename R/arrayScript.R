#!/bin/env Rscript
# arrayScript.R
library(somelibrary) # We don't want to reload the libraries every time that increases the individual compute time a lot
library(scone)
library(scran)


id <- Sys.getenv("SGE_TASK_ID")
subsample <- readRDS(toString(id)) # Might have to add directory to this
scored_subsample = scone(
  subsample,
  scaling = subsample@scaling_fn,
  k_qc = 8,
  k_ruv = 8,
  run = TRUE,
  zero = "postadjust",
  stratified_pam = FALSE,
  stratified_cor = FALSE,
  stratified_rle = FALSE,
  eval_kclust = 2:10,
  eval_pcs = 10,
  verbose = FALSE
)

saveRDS(scored_subsample, file =toString(id)) # save over it 

q(save="no") # End task