#' initialize_single_chr_CpG_df
#'
#'
#' @param chr chromosome name
#' @param num_spots number of spots
#' @export


initialize_single_chr_CpG_df <- function(chr, num_spots ){
  results<-data.frame(as.matrix(rep(chr,  num_spots, nrow=num_spots)))
  colnames(results)[1] <- "chr"
  results$start_pos <-rep(NA, nrow(results))
  results$stop_pos <-rep(NA, nrow(results))
  results$numCG <- rep(NA, nrow(results))
  results$numCG_positive <- rep(NA, nrow(results))
  results$numCG_negative <- rep(NA, nrow(results))
  
  results$unsigned_bin_score <- rep(NA, nrow(results))
  results$pos_bin_score <- rep(NA, nrow(results))
  results$neg_bin_score <- rep(NA, nrow(results))
  
  results$unsigned_bin_sig <- rep(NA, nrow(results))
  results$pos_bin_sig <- rep(NA, nrow(results))
  results$neg_bin_sig <- rep(NA, nrow(results))    
  
  return(results)
}