#' Score_CpG_bin 
#'
#' This function is a cludgy fix for a poorly designed FDR utility.
#' (NOTE THIS METHOD IS RATHER CRUDE, COULD STAND TO BE IMPROVED... SUGGESTIONS WELCOME :)
#'
#' @param CG_scores vector of scores e.g. -log10pvalue
#' @param FDR_scaler score at which FDR is achieved 
#' @export

Score_CpG_bin <- function(numCG, CG_scores, FDR_scaler=1, numCG_cutoff=1 ){
  
  if(numCG <= numCG_cutoff){
    bin_score <- 0
  } else { ## numCG > 1 
    bin_score <- (sum( (CG_scores/FDR_scaler)^2 ) - (max(CG_scores)/FDR_scaler)^2) / (numCG-1)      
    #bin_score <- (sum( (CG_scores/FDR_scaler)^2 ) - (max(CG_scores)/FDR_scaler)^2) / (numCG/(FDR_scaler)^2)
    if(is.infinite(bin_score) | is.na(bin_score)){print("Inf or NA error in score"); print("CG_scores are: "); print(CG_scores) } 
  }
  return(bin_score)
}