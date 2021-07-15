#' bnl_which_overlap
#'
#'
#' @export
#'
#  returns list which rows in prev_layer overlap those in window layer

bnl_which_overlap <- function(window, prev_layer){
  which_overlap <- rep(0, nrow(prev_layer))
    for( i in 1:nrow(prev_layer)){
      overlaps_any <- which(
        (as.character(prev_layer$chr[i]) == as.character(window$chr))
        & ((prev_layer$start_index[i] <= window$start_index & prev_layer$stop_index[i] >= window$stop_index )   ## prev layer completely covers window
        | (prev_layer$stop_index[i] >= window$start_index & prev_layer$stop_index[i] <= window$stop_index )  ## prev layer overlaps left window
        | (prev_layer$start_index[i] >= window$start_index & prev_layer$stop_index[i] <= window$stop_index ) ## window completely covers prev layer
        | (prev_layer$start_index[i] >= window$start_index & prev_layer$start_index[i] <= window$stop_index) )
      )
      if(length(overlaps_any) != 0 ){ which_overlap[i] <- 1}
    }
  return(which(which_overlap==1))
}
