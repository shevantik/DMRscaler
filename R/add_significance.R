#' add_significance
#'
#'
#' @param 
#' @param 
#' @export


add_significance<-function(three_window_list, lookup_table){
  for(i in 1:length(three_window_list)){
      if(is.null(nrow(three_window_list[[i]]))){next} ## bug fix for empty list
      if(nrow(three_window_list[[i]]) ==0){next} ## bug fix for empty list
    for(j in 1:nrow(three_window_list[[i]])){
      #print(j)
      three_window_list[[i]]$unsigned_bin_sig[j] <- estimate_peak_significance(num_CG=three_window_list[[i]]$numCG[j],
                                                                                    score=three_window_list[[i]]$unsigned_bin_score[j],
                                                                                    lookup_table = lookup_table )
    }
  }
  return(three_window_list)
}

