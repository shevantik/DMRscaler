#' bnl_build_overlaps
#'
#'
#' @param
#' @param
#' @export
#'
#  returns list of CG index sets for each prev_layer DMR object that overlaps the windows DMR object
bnl_build_overlaps <- function(window, prev_layer){
  which_overlap <- which(
    (prev_layer$start_index <= window$start_index & prev_layer$stop_index >= window$stop_index)   ## prev layer completely covers window
    | (prev_layer$stop_index >= window$start_index & prev_layer$stop_index <= window$stop_index )  ## prev layer overlaps left window
    | (prev_layer$start_index >= window$start_index & prev_layer$stop_index <= window$stop_index ) ## window completely covers prev layer
    | (prev_layer$start_index >= window$start_index & prev_layer$start_index <= window$stop_index)
  )
  window_indices <- seq(from=window$start_index, to=window$stop_index)
  overlap_indices <- list()
  j=1
  for(i in which_overlap){
    overlap_indices[[j]] <- seq(from=prev_layer[i,]$start_index, to=prev_layer[i,]$stop_index)
    overlap_indices[[j]] <- intersect(window_indices, overlap_indices[[j]])
    j<-j+1
  }
  return(overlap_indices)
}
