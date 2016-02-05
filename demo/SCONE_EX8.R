source("~/YosefCode/packages/RNASeq/OFBIT/SCONE/SCONE.R")

## ----- Standardize Quality Matrix -----
# Standardize quality metric
# q = quality metric matrix (columns = named features, rows = samples)
# ... = lists for specific transformations (see below)

PPQual = function(q, to.log = c("NREADS", "NALIGNED"), 
                  to.abs.log = c("MEDIAN_5PRIME_TO_3PRIME_BIAS","MEDIAN_5PRIME_BIAS","MEDIAN_3PRIME_BIAS"),
                  to.logit.one = c("PCT_RIBOSOMAL_BASES","PCT_CODING_BASES","PCT_UTR_BASES",
                                   "PCT_INTRONIC_BASES","PCT_INTERGENIC_BASES","PCT_MRNA_BASES"),
                  to.logit.hund = c("RALIGN")){
  
  ## ===== SD_EPSILON: Constant for purpose of correlation computation =====
  SD_EPSILON = 1e10 * .Machine$double.eps #~2.2e-6
  
  if(any(is.null(colnames(q)))){
    stop("No quality parameter names.")
  }
  
  quality.features = q
  
  # Convert non-numeric data fields
  for (i in 1:dim(quality.features)[2]){
    if(!is.numeric(quality.features[,i])){
      quality.features[,i] = as.numeric(as.character(quality.features[,i]))
    }
  }
  
  ## ----- Special Transformations
  # 0 to Infinity -> log
  quality.features[,to.log] = log(quality.features[,to.log]+.01)
  # 0 to Infinity, Best at 1
  quality.features[,to.abs.log] = exp(abs(log(quality.features[,to.abs.log]+.01)))-.01
  # 0 to 1
  quality.features[,to.logit.one] = log((quality.features[,to.logit.one]+.01)/(1-quality.features[,to.logit.one]+.01))
  # 0 to 100
  quality.features[,to.logit.hund] = log((quality.features[,to.logit.hund]+1)/(100-quality.features[,to.logit.hund] + 1))
  
  ## ----- Remove NA (missing data), Constants, and scale
  quality.features = t(na.omit(t(quality.features)))
  quality.features = quality.features[,apply(quality.features,2,sd) > SD_EPSILON]
  quality.features = scale(quality.features,center = T,scale = T)
  return(quality.features)
}

batch = read.table("/data/yosef/users/mbcole/BRAIN_DAVIDE/batch.txt")$V2
condition = read.table("/data/yosef/users/mbcole/BRAIN_DAVIDE/cell_type.txt")$V2
qual = PPQual(read.table("/data/yosef/users/mbcole/BRAIN_DAVIDE/filtered_qc.txt"))
e = as.matrix(read.table("/data/yosef/users/mbcole/BRAIN_DAVIDE/filtered_counts.txt"))
cc = read.table("/data/yosef/users/mbcole/WTP63/First/FP/cell_cycle_Tirosh.txt")$V2
hk = read.table("~/YosefCode/packages/RNASeq/summary/EXAMPLE/reference_files/house_keeping_human_names.txt")$V1
de = read.table("/data/yosef/users/mbcole/BRAIN_DAVIDE/macklis_markers.txt")$V1

batch = batch[condition == "L5"]
qual = qual[condition == "L5",]
e = e[,condition == "L5"]
condition = condition[condition == "L5"]

scone_out = SCONE(e = e,factor_free_only = F,design = c("nested"), nested_model = c("fixed"),
                  condition = condition,
                  batch = batch,
                  qual = qual,
                  dim_UV = 5,
                  is_HK = toupper(rownames(e)) %in% as.matrix(hk),
                  is_CC = toupper(rownames(e)) %in% as.matrix(cc),
                  is_DE = toupper(rownames(e)) %in% as.matrix(cc),
                  out_dir = "/data/yosef/users/mbcole/BRAIN_DAVIDE/SCONE_out_DE_L5")