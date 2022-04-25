#' get_loc_fdr_table
#'
#' @param mat matrix with samples on columns and locs on rows
#' @param cases case indices or column names
#' @param controls control indices or column names
#' @param stat_test statistical test used for signifance estimation
#' @param fdr false discover rate desired
#' @param resolution number of bins to use for table
#'
#' @importFrom foreach foreach
#' @importFrom foreach getDoParWorkers
#' @importFrom foreach %dopar%
#'
#' @return dmrscaler returns a list of dmrs at each layer
#'
#' @export

get_loc_fdr_table <- function(mat, cases, controls, stat_test, fdr=0.1, resolution=100){
  temp <- split(1:nrow(mat),cut(1:nrow(mat),max(getDoParWorkers(),2),labels=F))
  mat_subs_list <- list()
  for(i in 1:length(temp)){
    mat_subs_list[[i]] <- mat[temp[[i]],]
  }
  tp_pval <- foreach(mat_sub=mat_subs_list, .combine="c") %dopar% {
    tp_p <- numeric(length = nrow(mat_sub))
    for(i in 1:nrow(mat_sub)){
      tp_p[i] <- stat_test(mat_sub[i,cases],mat_sub[i,controls])$p.value
    }
    tp_p
  }



  num_permutations <- ceiling((1/fdr) * 2)
  if(num_permutations > choose(length(cases)+length(controls),length(cases) )){
    warning(paste("Warning: permutations required for accurate fdr estimation:",num_permutations,
                ",number of possible permutations:", choose(length(cases)+length(controls),length(cases) ),
                ",using ",choose(length(cases)+length(controls),length(cases) ), "permutations", sep=""))
    perms <- combn(length(cases)+length(controls), length(cases))
    perms <- t(perms[,1:floor(ncol(perms)/2)])
  } else {
    perms <- matrix(nrow=num_permutations, ncol=length(cases))
    for(i in 1:nrow(perms)){
      perms[i,] <- sample(1:length(c(cases,controls)),length(cases))
    }
  }

  print(paste("Running ",num_permutations, " permutations across ", getDoParWorkers()," worker nodes",sep=""))

  pvals <- foreach(i = 1:nrow(perms), .combine="cbind") %dopar% {
    g1 <- perms[i,]
    g2 <- which(!is.element(1:length(c(cases,controls)),perms[i,]))
    temp_p <- numeric(length = nrow(mat))
    for(j in 1:nrow(mat)){
      temp_p[j] <- stat_test(mat[j,g1],mat[j,g2])$p.value
    }
    temp_p
  }


  temp_seq <- seq(floor(min(log10(c(tp_pval, as.vector(pvals))))),0,length.out= resolution)
  tp_temp <- hist(log10(tp_pval), temp_seq, plot=F)
  tp_seq <- cumsum(tp_temp$counts / sum(tp_temp$counts))

  temp_seq <- seq(floor(min(log10(c(pvals, as.vector(pvals))))),0,length.out= resolution)
  temp <- hist(log10(pvals), temp_seq, plot=F)
  p_seq <- cumsum(temp$counts / sum(temp$counts))

  return(data.frame(log10pval_cutoff=temp_seq[-length(temp_seq)], fdr=p_seq/tp_seq))

}
#
