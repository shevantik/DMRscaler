#' n_nearest_window_scoring_func 
#'
#'
#' @param indat object that needs to be standardized and documented 
#' @param n_nearest number of adjacent cgs to consider
#' @param step_size step fraction = how many windows cover each cg
#' @param FDR_scaler might not be important, controls type 1 error sort of
#' @export


n_nearest_window_scoring_func <- function(indat, n_nearest, step_size, FDR_scaler){ 
  
  indat$chr
  noexport<-c("")
  export<-c("results")
  chrs <- levels(indat$chr)
  n_nearest_results<-foreach(i=chrs, .noexport = noexport, .export = export ) %dopar% {
    which_sub<-which( indat$chr == i)
    chr_indat <- indat[which_sub,]
    spots <- seq(from=1, to=length(chr_indat$scoring_value)-n_nearest, by=step_size)
    results <- initialize_single_chr_CpG_df(chr=i,
                                            num_spots=ceiling((length(chr_indat$scoring_value )-n_nearest)/step_size ) )
    
    j2=1
    for(j in 1:ceiling((length(chr_indat$scoring_value )-n_nearest)/step_size ) ){
      sub_data<-chr_indat[j2:(j2+n_nearest), ]
      j2=j2+step_size
      if(!is_empty(which(is.na(sub_data))) ){
        results[j,] <- initialize_row_CpG_df(chr=chrs[i], start_pos = -1, stop_pos = -1)
        next
      }
      results[j,]<-fill_bin_results(results = results[j,],
                                    sub_data = sub_data, start_pos =  min(as.numeric(as.character( sub_data$pos))),
                                    stop_pos = max(as.numeric(as.character( sub_data$pos))),
                                    FDR_scaler = FDR_scaler)
    }
    results<-results[which(results$numCG >0),]   
  }
}