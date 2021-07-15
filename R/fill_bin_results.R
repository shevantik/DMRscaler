#' fill_bin_results
#'
#'
#' @export


fill_bin_results<-function(results, sub_data, start_pos, stop_pos, FDR_scaler){
  results$start_pos <- start_pos
  results$stop_pos <- stop_pos

  results$numCG <- nrow(sub_data)
  results$numCG_positive <- length(which(   (sub_data$sign > 0) ) )
  results$numCG_negative <- length(which( (sub_data$sign < 0) ) )



  if(results$numCG <1){
    results$unsigned_bin_score <- 0
    results$pos_bin_score <- 0
    results$neg_bin_score <- 0
    next
  } else {
    unsigned_results <- Score_CpG_bin(numCG = results$numCG, CG_scores = sub_data$scoring_value, FDR_scaler = FDR_scaler)
  }
  results$unsigned_bin_score <- unsigned_results
  ### bin score stuff for pos signed significance peaks
  which_pos<-which(sub_data$sign > 0 )
  pos_sub_data <- sub_data[which_pos,]
  if(results$numCG_positive < 1){
    pos_results <- Score_CpG_bin_shortcircuit()
  } else {
    pos_results <- Score_CpG_bin(numCG = results$numCG_positive, CG_scores = pos_sub_data$scoring_value, FDR_scaler = FDR_scaler)
  }
  results$pos_bin_score <- pos_results
  ### bin score stuff for neg signed significance peaks
  which_neg <- which(sub_data$sign < 0)
  neg_sub_data <- sub_data[which_neg,]
  if(results$numCG_negative <1 ){
    neg_results <- Score_CpG_bin_shortcircuit()
  } else {
    neg_results<-Score_CpG_bin(numCG = results$numCG_negative, CG_scores = neg_sub_data$scoring_value, FDR_scaler = FDR_scaler)
  }
  results$neg_bin_score <- neg_results
  return(results)
}
