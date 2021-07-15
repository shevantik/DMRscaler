#' bnl_select_max_scoring_overlap
#'
#'
#' @export
#'
bnl_select_max_scoring_overlap<-function(overlaps, CG_table, FDR_scaler){
  max_score <- 0
  max_at <- 0
  if(length(overlaps)==0){return(overlaps)} ## in case overlaps is empty list (no overlapping)
  for( i in 1:length(overlaps)){
    score <- raw_score_CpG_bin(CG_scores =  CG_table$scoring_value[overlaps[[i]]], FDR_scaler=FDR_scaler )
    if(score >= max_score){max_score<-score; max_at <- i;}
  }
  return(overlaps[[max_at]])
}
