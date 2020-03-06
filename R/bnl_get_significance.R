#' bnl_get_significance
#'
#'
#' @param
#' @param
#' @export
bnl_get_significance <- function(score, lookup_table, num_cgs){
  low <- max(which(lookup_table$draw_count <= num_cgs)) ## lower CG bound
  high <- min(which(lookup_table$draw_count >= num_cgs )) ## upper CG bound
  sig_index<-min(max(which(lookup_table[low,] < score )),max(which(lookup_table[high,] < score )))
  sig<-as.numeric(colnames(lookup_table)[sig_index])
  return(sig)
}
