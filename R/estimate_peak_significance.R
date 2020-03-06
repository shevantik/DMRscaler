#' estimate_peak_significance
#'
#'
#' @param 
#' @param 
#' @export

estimate_peak_significance <- function(num_CG, score, lookup_table){
  if(num_CG < 2){return(0)}
  lower_row<-max(which(as.numeric(row.names(lookup_table)) <= num_CG  ))
  upper_row<-lower_row+1
  lower_row_sig <- as.numeric(colnames(lookup_table)[max(which(lookup_table[lower_row,] <= score))])
  upper_row_sig <- as.numeric(colnames(lookup_table)[max(which(lookup_table[upper_row,] <= score))])
  return(min(lower_row_sig,upper_row_sig))
}