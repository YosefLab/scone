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

data = tf.sc.eSet

# Load technical batch data
techbatch = c();
for(i in seq_along(colnames(data)))
{
  colname = colnames(data)[i];
  if(grepl("B04", colname)){
    techbatch = append(techbatch, "Week0");
  }
  else if(grepl("B17", colname)){
    techbatch = append(techbatch, "Week150");
  }
  else{
    stop("Didn't find B04 or B17 in sample label: ", colname)
  }
}
techbatch = as.factor(techbatch);

# Load biological batch data
biobatch = c();
for(i in seq_along(colnames(data)))
{
  colname = colnames(data)[i];
  if(grepl("A3001", colname)){
    biobatch = append(biobatch, "CMV");
  }
  else if(grepl("B3701", colname)) {
    biobatch = append(biobatch, "HIV");
  }
  else {
    stop("Didn't find A3001 or B3701 in sample label: ", colname);
  }
}
biobatch = as.factor(biobatch);

batch = as.factor(paste(techbatch,biobatch))
condition = biobatch
qual = PPQual(pData(protocolData(data)))
e = as.matrix(exprs(data)[gf.vec,])
cc = read.table("/data/yosef/users/mbcole/WTP63/First/FP/cell_cycle_Tirosh.txt")$V2
hk = read.table("~/YosefCode/packages/RNASeq/summary/EXAMPLE/reference_files/house_keeping_human_names.txt")$V1

scone_out = SCONE(e = e,factor_free_only = F,design = c("nested"), nested_model = c("fixed"),
                  condition = condition,
                  batch = batch,
                  qual = qual,
                  dim_UV = 5,
                  is_HK = rownames(e) %in% as.matrix(hk),
                  is_CC = rownames(e) %in% as.matrix(cc),
                  out_dir = "/data/yosef2/HIV/dat/CMV_HIV_150wk/SCONE_out_rerun")