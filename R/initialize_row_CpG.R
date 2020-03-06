#' initialize_row_CpG_df
#'
#'
#' @export




initialize_row_CpG_df <- function(chr = 0, start_pos = 0, stop_pos = 0, numCG = 0, numCG_positive = 0, 
                                  numCG_negative = 0 ,
                                  unsigned_bin_score = 0, 
                                  pos_bin_score = 0, neg_bin_score = 0, unsigned_bin_sig =0,
                                  pos_bin_sig = 0, neg_bin_sig = 0 ){
  results<-data.frame(chr=chr)
  results$start_pos <- start_pos
  results$stop_pos <- stop_pos
  results$numCG <- numCG
  results$numCG_positive <- numCG_positive
  results$numCG_negative <- numCG_negative
  
  results$unsigned_bin_score <- unsigned_bin_score
  results$pos_bin_score <- pos_bin_score
  results$neg_bin_score <- neg_bin_score
  
  results$unsigned_bin_sig <- unsigned_bin_sig
  results$pos_bin_sig <- pos_bin_sig
  results$neg_bin_sig <- neg_bin_sig      
  
  return(results)
}
