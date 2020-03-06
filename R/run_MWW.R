#' Run Mann Whitney Test
#'
#' This function runs the Mann Whitney U test
#' to generate significance values for difference
#' between case and control groups for each CG in the
#' Beta matrix. (SHOULD DEFINE HOW TO HANDLE TIES)
#'
#' @param control_indices indices of controls in Beta matrix
#' @param case_indices indices of cases in Beta matrix
#' @param Beta matrix of beta values, proportion methylated, with rows as CGs and cols as samples
#' @export

run_MWW <- function(control_indices, case_indices, Beta ){
  p_val<-vector(mode = "numeric", length = nrow(Beta))

  for(i in 1:nrow(Beta)){
    temp<-wilcox.test(
      Beta[i,case_indices],
      y = Beta[i,control_indices],
      paired = FALSE
    )
    p_val[i]<-temp$p.value
  }
  MWW<-data.frame(p_val)
  rownames(MWW)<-rownames(Beta)
  return(MWW)
}
