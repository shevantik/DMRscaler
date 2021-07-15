#' in_layer_does_overlap
#'
#'
#' @export


in_layer_does_overlap<-function(left, right){
  does_overlap<-((left$start_index <= right$start_index & left$stop_index >= right$stop_index)   ## left completely covers right
                 | (left$stop_index >= right$start_index & left$stop_index <= right$stop_index )  ## left overlaps left right
                 | (left$start_index >= right$start_index & left$stop_index <= right$stop_index ) ## right completely covers left
                 | (left$start_index >= right$start_index & left$start_index <= right$stop_index)
  )
  return(does_overlap)
}

in_layer_is_adjacent<-function(left,right){
  return(left$stop_index+1 == right$start_index)
}
