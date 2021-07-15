#' bnl_check_significance
#'
#'
#' @export

bnl_check_significance <- function(score, lookup_table, num_cgs, quants_significance_cutoff="0.999"){
### CODE NEEDS A DEEP CLEAN ### Lookup significance in table, is better than cutoff?
  quants<-colnames(lookup_table)
  quants_index<-which(quants==quants_significance_cutoff)
  significance_cutoff <-max(lookup_table[,quants_index][c(max(which(lookup_table$draw_count <= num_cgs)), min(which(lookup_table$draw_count > num_cgs )))])
  if(score >= significance_cutoff){return(TRUE)}
  else{return(FALSE)}
  #
}
