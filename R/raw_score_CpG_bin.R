#' raw_score_CpG_bin
#'
#'
#' @param
#' @param
#' @export
#'

raw_score_CpG_bin <- function(CG_scores, FDR_scaler=1){
  bin_score <- sum( (CG_scores/FDR_scaler)^2 )
  if(is.infinite(bin_score) | is.na(bin_score)){
    print("Inf or NA error in score"); print("CG_scores are: "); print(CG_scores)
  }
  return(bin_score)
}
