generate_subsamples <- function(base=2, start=1, end=6, per=16, scone_object, dir = '~/Desktop/subsampled_scone_objects/'){
  dir.create(dir)
  d = 1
  for(exp in start:end){
    percent = base^exp
    for(i in 1:per){
      sub = subsample_scone(scone_obj, seed = i, subsample_cell_level = percent,
                            at_bio = FALSE, keep_all_control = TRUE,
                            verbose = FALSE)
      saveRDS(sub, file = paste(dir, toString(d), sep=""))
      d = d+1
    }
  }
}

compare_scores <- function(original_scone, scone_list){
  original = get_score_ranks(original_scone)
  scone_list = unlist(scone_list)
  
  scores_df <- lapply(scone_list,get_score_ranks)
  cor_by_order <- function(row){
    return(cor(original, row[names(original)]))
  }
  corrs <- lapply(scores_df, cor_by_order)
  corrs <- unlist(corrs)
  return(corrs)
}



compare_subsample <- function(original_scone, 
                              min_percent_power, max_percent_power, 
                              min_seed, max_seed
){
  original_scores = get_score_ranks(original_scone)
  percents = 2 ^ seq(from=min_percent_power, to = max_percent_power, by=1)
  seeds = runif(2^ max_seed,0,100)
  list_percents <- c()
  number_scored = 0
  start_time = format(Sys.time(), "%r")
  
  for(percent in percents){
    # generate a random seed array
    list_of_seeds <- c()
    for(seed in seeds){
      
      # Generate the subsample
      print("Starting subsample")
      print(format(Sys.time(), "%r"))
      subsample <- subsample_scone(original_scone, subsample_cell_level = percent, 
                                   at_bio = FALSE, verbose = FALSE, seed = seed)
      print('Subsampled')
      print(format(Sys.time(), "%r"))
      # Score that Subsample
      scored_subsample = scone(
        subsample,
        scaling = original_scone@scaling_fn,
        k_qc = 8,
        k_ruv = 8,
        run = TRUE,
        zero = "postadjust",
        stratified_pam = FALSE,
        stratified_cor = FALSE,
        stratified_rle = FALSE,
        eval_kclust = 2:10,
        eval_pcs = 6,
        verbose = FALSE
        
      )
      number_scored = number_scored + 1
      print("Number Scored")
      print(toString(number_scored))
      print('Scored')
      print(format(Sys.time(), "%r"))
      # Add the metadata
      #scored_subsample$metadata$'percent' = percent
      #scored_subsample$metadata$'seed' = seed
      # Append to the list of seeds
      list_of_seeds <- c(list_of_seeds, scored_subsample)
      
    }
    list_percents <- c(list_percents, list_of_seeds)
  }
  
  # subsampled_scores = get_score_ranks(subsampled_scone_scored)
  # return(plot(original_scores, subsampled_scores[names(original_scores)]))
  
  
  # Reshape list
  matrixed <- matrix(list_percents, ncol = max_percent_power - min_percent_power +1, byrow = FALSE)
  df <- as.data.frame(matrixed)
  names(df) <- percents
  print("Started At")
  print(start_time)
  print("Ended At")
  print(format(Sys.time(), "%r"))
  return(df)
}


compare_subsample_pca<- function(original_scone, subsampled_scone, method){
  # Get the normalizations
  subsampled_normalized <- get_normalized(subsampled_scone, method = method, log= TRUE)
  original_normalized <- get_normalized(original_scone, method = method, log= TRUE)
  
  # Calculate PCA of subsample and apply it to original
  pca_sub <- prcomp(subsampled_normalized, scale=TRUE, center=TRUE)
  pca_original <- predict(pca_sub, original_normalized)[,1:10]
  
  # This is useless because it jsut shows the whole point of principal components
  
  # Get proportions of Variance for original
  variance_orignals <- apply(pca_original,2, sd) ^ 2
  sum_variances_original <- sum(variance_orignals)
  proportions_of_variance_original <- variance_orignals / sum_variances_original
  
  # Get proportions of Variance for subsample
  variance_subbed <- pca_sub$sdev[1:10] ^ 2
  sum_variances_sub <- sum(variance_subbed)
  proportions_of_variance_subbed <- variance_subbed / sum_variances_sub
  
  # Plot the linreg of PC's
  df <- data.frame(proportions_of_variance_subbed, proportions_of_variance_original)
  #plot(lm(data = df))
  
  pca_original_2 <- prcomp(original_normalized, scale=TRUE, center=TRUE)
  
  # Return the proprtional difference between principal components
  # return(c(proportions_of_variance_original - proportions_of_variance_subbed 
  #   / proportions_of_variance_subbed, proportions_of_variance_original, proportions_of_variance_subbed, pca_sub, pca_original))
  return(pca_original_2)
  
}