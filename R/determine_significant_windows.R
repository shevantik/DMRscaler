#' determine_significant_windows
#'
#'
#' @param
#' @param
#' @export

determine_significant_windows<-function(window_results, indat, quants, quants_significance_cutoff = "0.9999"){
  chrs <- levels(indat$chr)
  quants_index<-which(colnames(quants)==quants_significance_cutoff)

  noexport<-c("")
  export<-c("sig_window_results")

  # Unsigned
  unsigned_significant_windows<-foreach(i=1:length(chrs), .noexport = noexport, .export = export ) %dopar% { #1:length(chrs)
    which_sub<-which( indat$chr == chrs[i])
    chr_indat <- indat[which_sub,]
    rm(which_sub)
    chr_window_results <- window_results[[i]]

    significance_cutoff <- rep(-1,  length(chr_window_results$unsigned_bin_score))
    for(j in 1:length(chr_window_results$unsigned_bin_score)){
      significance_cutoff[j] <- max(quants[,quants_index][c(max(which(quants$draw_count <= chr_window_results$numCG[j])),
                                                            min(which(quants$draw_count > chr_window_results$numCG[j] )))])
    }
    which_sig <- which(chr_window_results$unsigned_bin_score > significance_cutoff)

    sig_window_results<-chr_window_results[which_sig,]
  }
  return(unsigned_significant_windows)
}
