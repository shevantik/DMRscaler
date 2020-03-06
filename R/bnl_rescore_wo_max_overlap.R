#' bnl_rescore_wo_max_overlap
#'
#'
#' @param
#' @param
#' @export
#'
###bnl_rescore_wo_max_overlap(windows_to_layer[i,], max_overlap, CG_table)
## window = df 1 row, max_overlap = vector, CG_table = as is.
bnl_rescore_wo_max_overlap<-function(window, max_overlap, CG_table, FDR_scaler){
  window_indices <- seq(from=window$start_index, to=window$stop_index)
  indices_remaining <- setdiff(window_indices, max_overlap)
  if(length(indices_remaining)==0){return(0)}
  score <- score_CpG_bin_layers(CG_scores = CG_table$scoring_value[indices_remaining], FDR_scaler = FDR_scaler)
  return(score)
}
