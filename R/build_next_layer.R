#' build_next_layer
#'
#'
#' @param
#' @param
#' @export
#'

build_next_layer<-function(prev_layer, windows_to_layer, CG_table, lookup_table, FDR_scaler){
  ## assume layer below is non_overlapping non_adjacent windows
  ## for each DMR object (i.e. windows_to_layer)
  ## nrow(windows_to_layer)  == number of DMR objects
  keep <- rep(FALSE,nrow(windows_to_layer))
  for(i in 1:nrow(windows_to_layer)){
    #print(i)
    window<-windows_to_layer[i,]
    overlaps <- bnl_build_overlaps(window = window, prev_layer=prev_layer) ## returns list of overlaps
    max_overlap <- bnl_select_max_scoring_overlap(overlaps=overlaps, CG_table=CG_table, FDR_scaler = FDR_scaler)
    rescore <- bnl_rescore_wo_max_overlap(window=window, max_overlap = max_overlap, CG_table=CG_table, FDR_scaler = FDR_scaler)
    rescore
    is_sig <- bnl_check_significance(score = rescore, lookup_table = lookup_table, num_cg = (window$numCG-length(max_overlap)))
    if(is_sig){
      keep[i] <- TRUE
    }
  }
  kept <- windows_to_layer[which(keep),]
  if(nrow(kept)>0){
    kept_merged <- in_layer_merge(dmrs=kept, CG_table = CG_table, FDR_scaler = FDR_scaler, lookup_table = lookup_table)
  } else { kept_merged <- kept}


  kept_expanded <- expand_bounds(lower=prev_layer, upper=kept_merged, CG_table = CG_table,
                                 FDR_scaler = FDR_scaler, lookup_table = lookup_table)


  prev_unchanged <- prev_layer[-bnl_which_overlap(window = kept_expanded, prev_layer=prev_layer ),]  ## which in the previous layer don't overlap anything, note minus sign
  if(length(bnl_which_overlap(window = kept_expanded, prev_layer=prev_layer ))==0 ){ prev_unchanged <- prev_layer} ## error correction for handling cases where nothing overlaps

  kept_expanded <- rbind(kept_expanded, prev_unchanged)
  kept_expanded <- kept_expanded[with(kept_expanded, order(chr, start_pos)),]

  return(kept_expanded)
}
