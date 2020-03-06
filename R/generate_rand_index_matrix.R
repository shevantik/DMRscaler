#' Case Control Permutation Description Matrix
#'
#' This function generates a matrix which describes
#' series of case-control label permutations.
#' @param num_controls number of controls
#' @param num_cases number of cases
#' @param num_permutations number of permutations
#' @export

generate_rand_index_matrix <- function(num_controls, num_cases, num_permutations ){
  ##### index matrix contains 0 for rand_0, 1 for rand_1 ####
  index_matrix<-matrix(nrow = (num_controls+num_cases), ncol = num_permutations)
  colnam<-vector()
  for(i in 1:num_permutations){
    colnam[i]<-paste("rand_", i, sep = "")
    index_matrix[,i]<-sample(c(rep(0,num_controls),rep(1,num_cases)),(num_controls+num_cases), replace = FALSE)
  }
  colnames(index_matrix)<-colnam
  #### demographics of randoms can be reconstructed by looking up the 0s and 1s in the index matrix
  return(index_matrix)
}
