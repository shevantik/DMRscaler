#' in_layer_is_adjacent
#'
#'
#' @param
#' @param
#' @export

in_layer_is_adjacent<-function(left,right){
  return(left$stop_index+1 == right$start_index)
}
