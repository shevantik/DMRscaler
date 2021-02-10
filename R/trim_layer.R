#' trim_layer
#'
#'
#' @param
#' @param
#' @export

## merge adjacent and overlapping windows within the same layer
trim_layer<-function(dmrs,CG_table, FDR_scaler,lookup_table){
  ## scan through and merge as we go
  ## for i and i+1
  ## test if rows are equivalent, if yes next
  ## if overlapping then merge else next
  ## if merge pairwise merge
  ## pairwise merge: test setdiffs and intersect. Keep all above significance.
  ## replace original rows with the solution. At end of merge, remove all redundant rows and return
  CG_table$above_fdr <- CG_table$scoring_values >= FDR_scaler
  keep <- rep(1,nrow(dmrs))
  for(i in 1:(nrow(dmrs))){
    #print(i)
    chr <- dmrs[i,]$chr
    start <- dmrs[i,]$start_pos
    stop <- dmrs[i,]$stop_pos
    which <- which(CG_table$chr == chr & CG_table$pos >= start & CG_table$pos <= stop & CG_table$above_fdr)  ## above_fdr used to trim edges
    if(length(which)<=1){keep[i]<-0; next;}
    which <- min(which):max(which)
    
    
    
    cgs <- CG_table$scoring_value[which]
    num_cgs <- length(cgs)
    unsigned_bin_score<-score_CpG_bin_layers(CG_scores = cgs , FDR_scaler = FDR_scaler)
    unsigned_bin_sig<-bnl_get_significance(score=unsigned_bin_score, lookup_table = lookup_table, num_cgs = num_cgs )
    sig <- bnl_check_significance(score=unsigned_bin_score, lookup_table = lookup_table, num_cgs = num_cgs)
    if(!sig){keep[i]<-0; next}
    
    dmrs[i,]$start_pos <- min(CG_table[which,]$pos)
    dmrs[i,]$start_index <- min(which)
    dmrs[i,]$stop_pos <- max(CG_table[which,]$pos)
    dmrs[i,]$stop_index <- max(which)
    dmrs[i,]$numCG <- length(which)
    dmrs[i,]$unsigned_bin_score <- unsigned_bin_score
    dmrs[i,]$unsigned_bin_sig <- unsigned_bin_sig
    
    
  }
  return(dmrs[which(keep==1),])
}
