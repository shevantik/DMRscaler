#' expand_bounds
#'
#'
#' @export
#'


expand_bounds<-function(lower, upper, CG_table, FDR_scaler, lookup_table ){
  ## expand upper around lower
  ## should be able to assume significance doesn't drop, but may want to test this
  for(i in 1:nrow(upper)){
    which_lower_overlap_or_adj <- which( (lower$start_index < upper$start_index[i] & lower$stop_index >= upper$start_index[i])
                                         | (lower$start_index <= upper$stop_index[i] & lower$stop_index > upper$stop_index[i]))
    if(length(which_lower_overlap_or_adj)==0){next;}
    overlaps<-lower[which_lower_overlap_or_adj,]
    for(j in 1:nrow(overlaps)){
      if((overlaps[j,]$start_index < upper$start_index[i] & overlaps[j,]$stop_index >= upper$start_index[i])){  ## lower left overlap
        upper[i,] <- merge_adjacent(left = overlaps[j,], right = upper[i,] ,
                                    CG_table = CG_table, FDR_scaler = FDR_scaler,
                                    lookup_table = lookup_table)
      }
      if((overlaps[j,]$start_index <= upper$stop_index[i] & overlaps[j,]$stop_index > upper$stop_index[i])){  ## lower right overlap
        upper[i,] <- merge_adjacent(left = upper[i,], right = overlaps[j,] ,
                                    CG_table = CG_table, FDR_scaler = FDR_scaler,
                                    lookup_table = lookup_table)
      }
    }
  }
  return(upper)
}
