#' merge_adjacent
#'
#'
#' @param
#' @param
#' @export

merge_adjacent<- function(left,right, CG_table, FDR_scaler, lookup_table){
  start_index<-left$start_index
  stop_index<-right$stop_index
  CG_scores <- CG_table$scoring_value[seq(from= start_index,to=stop_index)]
  num_cgs <- length(CG_scores)
  unsigned_bin_score<-score_CpG_bin_layers(CG_scores = CG_scores , FDR_scaler = FDR_scaler)
  unsigned_bin_sig<-bnl_get_significance(score=unsigned_bin_score, lookup_table = lookup_table, num_cgs = num_cgs )
  ## repack left
  left$stop_pos<-right$stop_pos
  left$numCG<-num_cgs
  left$unsigned_bin_score<-unsigned_bin_score
  left$unsigned_bin_sig<-unsigned_bin_sig
  left$stop_index<-right$stop_index
  return(left)
}
