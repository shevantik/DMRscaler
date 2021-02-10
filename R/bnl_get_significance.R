#' bnl_get_significance
#'
#'
#' @param
#' @param
#' @export
bnl_get_significance <- function(score, lookup_table, num_cgs){
  low <- max(which(lookup_table$draw_count <= num_cgs)) ## lower CG bound
  high <- min(which(lookup_table$draw_count >= num_cgs )) ## upper CG bound
  if(length(which(lookup_table[high,-ncol(lookup_table)] < score ) )==0 ){ return(0)}
  if(length(which(lookup_table[low,-ncol(lookup_table)] < score ) )==0 ){ return(0)}
  sig_index<-min(max(which(lookup_table[low,-ncol(lookup_table)] < score )),max(which(lookup_table[high,-ncol(lookup_table)] < score )))
  sig<-as.numeric(colnames(lookup_table)[sig_index])
  return(sig)
}
