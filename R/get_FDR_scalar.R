#' Get_FDR_scalar
#'
#' This function is a cludgy fix for a poorly designed FDR utility.
#' (NOTE THIS METHOD IS RATHER CRUDE, COULD STAND TO BE IMPROVED... SUGGESTIONS WELCOME :)
#'
#' @param MWW_FDR_table very specifically formated table == BAD DESIGN
#' @param MWW_FDR_threshold e.g. 0.2 for 20% FDR
#' @export

### Accessing the fdr scalar is non trivial, would write
get_FDR_scalar <- function(MWW_FDR_table, MWW_FDR_threshold){
  MWW_FDR_scalar <- as.numeric(strsplit(rownames(as.data.frame(MWW_FDR_table)[min(which(as.data.frame(MWW_FDR_table)$FDR < MWW_FDR_threshold)),]), "=")[[1]][2])
  return(MWW_FDR_scalar)
}
