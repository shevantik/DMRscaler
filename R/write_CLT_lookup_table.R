#' Write CLT lookup table
#'
#' This function generates a lookup table
#' to estimate the likelihood that a window
#' represents a DMR feature (i.e. appears enriched in differntially methylated CGs)
#' (NOTE THIS METHOD IS RATHER CRUDE, COULD STAND TO BE IMPROVED... SUGGESTIONS WELCOME :)
#'
#' @param
#' @param
#' @param
#' @export



### ADD ROW FOR INFINITY TO AVOID CONFUSING ERROR!!!
### (NOT DONE YET: ADD -1 as score for 0th percentile, ADD Inf as Score for 100th percentile)
write_CLT_lookup_table <- function( num_reps, data_to_sample, FDR_scaler, clt_numCGs = c(2,5,10,25,50,100, 150, 200, 250, 1000),
                                    quants=quants<-data.frame(row.names=c(0,0.5,0.75,0.9,0.95,0.99, 0.999, 0.9995, 0.9999, 0.99995, 0.99999, 1))
){
  # clt_numCGs <-c(2,5,10,25,50,100, 150, 200, 250, 1000) ### ADD INFINITY ROW!!!! COPY LAST ROW FOR INFINITY ROW!!!
  clt_time_start<-timestamp()
  numCores<-getDoParWorkers()
  perCore<-ceiling(num_reps/numCores)

  for(i in 1:length(clt_numCGs)){

    #print(clt_numCGs[i]) ######
    noexport=c("")
    export=c("clt")
    clt<-foreach(j=1:numCores, .noexport = noexport, .export = export, .combine = 'c' ) %dopar% {
      #print("hello")
      clt <- rep(-1,perCore)
      for(k in 1:perCore){
       # print(k)
        temp_sample<-base::sample(data_to_sample, clt_numCGs[i] )
        clt[k]<-Score_CpG_bin(numCG = length(temp_sample), CG_scores = temp_sample, FDR_scaler = FDR_scaler)
      }
      clt
    }
    clt<-clt[1:num_reps]
    # quants<-cbind(quants, as.data.frame(quantile(clt, probs = c(0,0.5,0.75,0.9,0.95,0.99, 0.999, 0.9995, 0.9999, 0.99995, 0.99999, 1))))
    quants<-cbind(quants, as.data.frame(quantile(clt, probs = as.numeric(row.names(quants)))))
  }


  colnames(quants)=clt_numCGs
  quants<-as.data.frame(t(quants))
  quants$draw_count<-as.numeric(rownames(quants))
  ## Add -1 as first row to make error obvious
  quants<-rbind("-1"=rep(-1,ncol(quants)), quants)
  quants<-rbind( quants,"Inf"=rep(0,ncol(quants))) ## attempted fix
  quants$draw_count[nrow(quants)]<-Inf
  quants$`0`<-rep(0,nrow(quants))
  quants$`1`<-rep(Inf,nrow(quants))
  clt_time_stop<-timestamp()
  clt_time_start_last<-clt_time_start
  return(quants)
}
