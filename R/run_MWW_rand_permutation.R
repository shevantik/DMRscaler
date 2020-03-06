#' Run Case Control Permutation
#'
#' This function runs the a series of
#' case-control label permutations. Recording
#' the result from each run
#' @param index_matrix matrix of permuted case-control labels
#' @param Beta matrix of methylation proportions unmethylated/methylated
#' @param num_permutations number of permutations
#' @export


run_MWW_rand_permutation<-function(index_matrix, Beta, num_permutations){
  export<-c("MWW_rand")
  noexport<-""
  MWW_rand<-foreach(i=1:num_permutations, .combine = cbind, .noexport = noexport, .export = export ) %dopar% {
    rand_0_index<-which(index_matrix[,i]==0)
    rand_1_index<-which(index_matrix[,i]==1)
    p_val_rand<-vector(mode = "numeric", length = nrow(Beta))
    for(j in 1:nrow(Beta)){
      temp<-wilcox.test( Beta[j,rand_0_index], y = Beta[j,rand_1_index],paired = FALSE)
      p_val_rand[j]<-temp$p.value
    }
    MWW_rand<-data.frame(a=p_val_rand)
    names(MWW_rand)[names(MWW_rand)=="a"]<-paste("p_val_rand_", i, sep = "")
    MWW_rand
  }
  return(MWW_rand)
}
