#' in_layer_merge
#'
#'
#' @param
#' @param
#' @export

## merge adjacent and overlapping windows within the same layer
in_layer_merge<-function(dmrs,CG_table, FDR_scaler,lookup_table){
  ## scan through and merge as we go
  ## for i and i+1
  ## test if rows are equivalent, if yes next
  ## if overlapping then merge else next
  ## if merge pairwise merge
  ## pairwise merge: test setdiffs and intersect. Keep all above significance.
  ## replace original rows with the solution. At end of merge, remove all redundant rows and return

  if(nrow(dmrs)<=1){return(dmrs)} ## trivial don't need to merge anythign
  for(i in 1:(nrow(dmrs)-1)){
    #print(i)
    left<-dmrs[i,]
    right<-dmrs[i+1,]
    if(left$start_index == right$start_index & left$stop_index == right$stop_index){next;} ## skip if equivalent
    if(!in_layer_does_overlap(left,right) & !in_layer_is_adjacent(left,right)){next;} ## skip if non-overlapping and non-adjacent
    if(in_layer_is_adjacent(left,right)){
      merged <- merge_adjacent(left = left,right = right,CG_table = CG_table, FDR_scaler=FDR_scaler,lookup_table = lookup_table)
      dmrs[i,]<-merged
      dmrs[i,]$numCG<-0
      dmrs[i+1,]<-merged
      next
    }
    ### partial overlap case (left set diff ...---___ right set diff)
    ### if score of left set diff AND right set diff are both significant then merge to union
    ### else return max scoring windows (left or right)

    left_seq<-seq(from=left$start_index, to = left$stop_index)
    right_seq<-seq(from=right$start_index, to = right$stop_index)


    left_diff<-setdiff(left_seq, right_seq)
    right_diff<-setdiff(right_seq, left_seq)

    left_cgs <- CG_table$scoring_value[left_diff]
    left_num_cgs <- length(left_cgs)
    left_unsigned_bin_score<-score_CpG_bin_layers(CG_scores = left_cgs , FDR_scaler = FDR_scaler)
    left_unsigned_bin_sig<-bnl_get_significance(score=left_unsigned_bin_score, lookup_table = lookup_table, num_cgs = left_num_cgs )
    left_sig <- bnl_check_significance(score=left_unsigned_bin_score, lookup_table = lookup_table, num_cgs = left_num_cgs)

    right_cgs <- CG_table$scoring_value[right_diff]
    right_num_cgs <- length(right_cgs)
    right_unsigned_bin_score<-score_CpG_bin_layers(CG_scores = right_cgs , FDR_scaler = FDR_scaler)
    right_unsigned_bin_sig<-bnl_get_significance(score=right_unsigned_bin_score, lookup_table = lookup_table, num_cgs = right_num_cgs)
    right_sig <- bnl_check_significance(score=left_unsigned_bin_score, lookup_table = lookup_table, num_cgs = left_num_cgs )

    ## for cases
    if(right_sig & left_sig){
      merged<-merge_adjacent(left = left,right = right,CG_table = CG_table, FDR_scaler=FDR_scaler,lookup_table = lookup_table)
      dmrs[i,]<-merged
      dmrs[i,]$numCG<-0
      dmrs[i+1,]<-merged
      next;
    }else{
      if(left_sig){dmrs[i+1,]<-left;next;}
      if(right_sig){dmrs[i,]$numCG<-0;next;}
      if(right$unsigned_bin_score > left$unsigned_bin_score){ dmrs[i,]$numCG<-0;next}
      else{ dmrs[i,]$numCG<-0; dmrs[i+1,]<-left;next }
    }

  }
  dmrs<-dmrs[which(dmrs$numCG!=0),]


  return(dmrs)
}
