#' Write FDR Table
#'
#' This function uses the permutation and real data partition to
#' estimate the leve at which the desired FDR is achieved.
#' (NOTE THIS PART OF THE METHOD IS SUBOPTIMAL, SHOULD EVENTUALLY IMPROVE... SUGGESTIONS WELCOME :)
#'
#' @param real_table
#' @param rand_table
#' @param FDR_bin_count
#' @param FDR_table_row_names
#' @export

### Estimate the false discovery rate (FDR) from the prior permutation scheme
### Input: real table (vector) of results from the real patient-control partition of the data
###        rand table (matrix) of resutls from the permutation partition of the data, found
###        by the permutation function(s).
###        FDR_bin_count is just how many equal width bins FDR should be reported in (bad definition should clean up)
###        FDR_table_row_names is prefix to give row names (bad definition should clean up)
write_FDR_table <-function(real_table, rand_table, FDR_bin_count=50, FDR_table_row_names="MWW="){

  FDR_bins_max<-max(c(max(real_table),max(rand_table)))
  FDR_bins_step<-FDR_bins_max/FDR_bin_count
  FDR_bins<-seq(0, FDR_bins_max+FDR_bins_step, by=FDR_bins_step)

  binned_score_vals<-matrix()
  FDR_score_bins<-vector("numeric")
  for(i in 1:length(FDR_bins)){
    temp<- abs(real_table)
    FDR_score_bins[i]<-length(which( (temp >= (FDR_bins[i] ))  ))
  }
  binned_score_vals<-cbind(FDR_score_bins)
  row.names(binned_score_vals)<-FDR_bins
  colnames(binned_score_vals)<-"patient_v_control"

  noexport<-c("FDR_score_bins","temp", "j")
  export<-c("binned_score_vals")
  a<-foreach(i=1:length(rand_table[1,]) , .combine = cbind, .noexport = noexport, .export = export ) %dopar% {
    FDR_score_bins<-vector("numeric")
    for(j in 1:length(FDR_bins)){
      temp<- abs(rand_table[,i])
      FDR_score_bins[j]<-length(which( (temp >= (FDR_bins[j] ))  ))
    }
    FDR_score_bins
  }
  binned_score_vals<-cbind(binned_score_vals, a)

  rmedians<-rowMedians(binned_score_vals[,2:ncol(binned_score_vals)])
  binned_score_vals<-cbind(binned_score_vals, random_medians=rmedians)
  FirstQ<-vector()
  ThirdQ<-vector()
  for(i in 1:nrow(binned_score_vals) ){
    temp<-summary(binned_score_vals[i,2:ncol(binned_score_vals)])
    FirstQ[i]<-temp["1st Qu."]
    ThirdQ[i]<-temp["3rd Qu."]
  }
  binned_score_vals<-cbind(binned_score_vals, FirstQuartile=FirstQ, ThirdQuartile=ThirdQ)
  FDR<-binned_score_vals[1:nrow(binned_score_vals),"random_medians"]/binned_score_vals[1:nrow(binned_score_vals),"patient_v_control"]
  binned_score_vals<-cbind(binned_score_vals, FDR=FDR)

  FDR_table<-binned_score_vals[, c("patient_v_control", "FirstQuartile", "random_medians", "ThirdQuartile", "FDR")]
  rownames(FDR_table)<- paste(FDR_table_row_names, rownames(FDR_table), sep = "")
  return(FDR_table)
}
