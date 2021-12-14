#' dmrscaler
#'
#'
#' @param locs dataframe of measured CpG loci with columns as "chr", "pos", "pval"
#' @param loc_signif_method one of "p-value" or "fdr" (false discovery rate)
#' @param loc_signif_cutoff p-value at which individual CpGs are considered significant or desired fdr is achieved for individual CpGs
#' @param region_signif_method one of "p-value" or "fdr" (false discovery rate) or "fwer" (family-wise error rate)
#' @param region_signif_cutoff p-value or fwer at which regions are considered significant or desired fdr is achieved for regions
#' @param window_type specifies  one of "k_nearest" or "genomic_width"
#' @param window_sizes vector of window size for each layer where window_type determines whether this represents the genomic_width of windows or the k_nearest neighbors that compose a window
#' @param dmr_constraint_list (ADD LATER)(optional) a list of rules that further constrains the definition of a differentially methylated region
#' @param output_type one of "simple" or "complete" where "simple" returns only the topmost layer as a dataframe of dmrs and "complete" returns a list of dataframes of dmrs that were called at and/or propegated up through each layer
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#'
#' @return dmrscaler returns a list of dmrs at each layer
#'

dmrscaler <- function(locs,
                      locs_pval_cutoff = 0.05,
                      region_signif_method = c("fwer","fdr","p-value"),
                      region_signif_cutoff = 0.01,
                      window_type = c("k_nearest", "genomic_width"),
                      window_sizes = c(2,4,8,16,32,64),
                      dmr_constraint_list = NULL,
                      output_type = c("simple", "complete")
){

  ### check input parameters are valid
  if( !all(is.element(c("chr","pos","pval"), colnames(locs))) ){
    stop("ERROR: locs must include column names: \"chr\",\"pos\",\"pval\". ")
  }
  stopifnot( !is.na(pmatch(region_signif_method, c("fwer","fdr","p-value")) ))
  region_signif_method <- match.arg(region_signif_method)
  stopifnot( !is.na(pmatch(window_type, c("k_nearest","genomic_width")) ))
  window_type <- match.arg(window_type)
  stopifnot( !is.na(pmatch(output_type, c("simple","complete")) ))
  output_type <- match.arg(output_type)
  stopifnot( locs_pval_cutoff > 0 & locs_pval_cutoff < 1 )
  stopifnot( region_signif_cutoff > 0 & region_signif_cutoff < 1 )
  stopifnot( all(window_sizes >= 2 ) )

  ## update locs to remove signal from locs with -log(p) < -log(cutoff) and set rank
  locs$pval[which(locs$pval > locs_pval_cutoff)] <- 1 ## set -log(p) to 0 if p is above cutoff
  locs$pval_rank <- rank(locs$pval, ties.method = "max") ## determine rank of each p value, used for region significance
  total_locs <- nrow(locs)

  ## organize locs into list of dataframes where each dataframe is from a unique chr
  locs_list <- list()
  for(chr in unique(locs$chr)){
    locs_list[[chr]] <- locs[which(locs$chr==chr),c("pos","pval","pval_rank")] ## separate locs by chromosome
    locs_list[[chr]] <- locs_list[[chr]][order(locs_list[[chr]]$pos),] ## order locs by position
    locs_list[[chr]]$start <- locs_list[[chr]]$pos  ## change pos to start and stop positions - generalized to work with meta_locs
    locs_list[[chr]]$stop <- locs_list[[chr]]$pos
    locs_list[[chr]]$pos <- NULL
    rownames(locs_list[[chr]]) <- NULL
  }


  ## build the primary output object i.e. dmr_layer_list
  dmr_layer_list <- list()
  for(window_index in 1:length(window_sizes)){
    window_size <- window_sizes[window_index]
    layer_name <- paste(window_size,"_loc_window_layer", sep="")
    dmr_layer_list[[layer_name]] <- foreach(chr_locs = locs_list, .final = function(x) setNames(x, names(locs_list))) %dopar% {
      ### first call features in layer using independently of all other layer
      chr_locs$in_dmr <- F
      which_signif <- which(chr_locs$pval < locs_pval_cutoff)
      which_signif_index <- 1

      ### paint all locs that are in a window that has significance < region_signif_cutoff
      while(TRUE){
        current_signif_index <- ifelse(which_signif_index > length(which_signif), -1, which_signif[which_signif_index])
        which_signif_index <- which_signif_index + 1
        ## test whether at end of array of significant locs, if yes break
        if(current_signif_index == -1 ){ break }
        right_index <- min( nrow(chr_locs), current_signif_index+window_size-1)
        right_signif_index <- max(intersect( which_signif, current_signif_index:right_index ))
        if(right_signif_index == current_signif_index){ next }

        ## use series of hypergeometric tests for significance
        window_locs <- chr_locs[current_signif_index:right_index,]
        window_locs <- window_locs[which(window_locs$pval < 1),] ## only retain significant locs, if pval==1, loop will multiple window_signif by 1 (doing nothing)
        window_loc_ranks <- window_locs$pval_rank[order(window_locs$pval_rank)][-1]  # [-1] drops most signif loc. Most signif loc serves as prior
        window_signif <- 1
        n <- total_locs
        k <- length(current_signif_index:right_index)-1  ## realized window_size minus most signif loc
        for(i in length(window_loc_ranks):1 ){
          window_signif = window_signif * dhyper(x=i, m=window_loc_ranks[i], n=max(0,n-window_loc_ranks[i]), k=k)
          n <- window_loc_ranks[i]-1
          k <- i-1
        }
        if(window_signif < region_signif_cutoff){
          chr_locs$in_dmr[current_signif_index:right_signif_index] <- TRUE
        }
      }

      ## build dmrs
      ## add dmrs ranges
      dmrs <- data.frame(start=numeric(),stop=numeric(),pval_region=numeric() )
      next_dmr <- data.frame(start=-1,stop=-1,pval_region=-1 )
      for(i in 1:nrow(chr_locs)){
        if(chr_locs$in_dmr[i]){
          if(next_dmr$start == -1){
            next_dmr$start <- chr_locs$start[i]
          }
          if( i+1 > nrow(chr_locs) | !chr_locs$in_dmr[i+1] ){
            next_dmr$stop <- chr_locs$stop[i]
            dmrs <- rbind(dmrs, next_dmr)
            next_dmr <- data.frame(start=-1,stop=-1,pval_region=-1 )
          }
        }
      }



      ## add dmrs significance
      if(nrow(dmrs) > 0 ){ ## need to test whether any dmrs were found
        for(i in 1:nrow(dmrs)){
          window_locs <- chr_locs[which(chr_locs$start==dmrs$start[i]):which(chr_locs$stop==dmrs$stop[i]),]
          k <- nrow(window_locs)   ## get current dmr size
          window_locs <- window_locs[which(window_locs$pval < 1),]
          window_loc_ranks <- window_locs$pval_rank[order(window_locs$pval_rank)][-1]
          window_signif <- 1

          n <- total_locs
          k <- max(window_size, k)  ## if current dmr size is not as large as window_size, use window_size for signif calc
          for(j in length(window_loc_ranks):1 ){
            window_signif = window_signif * dhyper(x=j, m=window_loc_ranks[j], n=max(0,n-window_loc_ranks[j]), k=k)
            n <- window_loc_ranks[j]-1
            k <- j-1
          }
          dmrs$pval_region[i] <- window_signif

        }
      }

      dmrs
    } ## end foreach

    ## update locs_list to replace locs added to dmrs with meta-locs
    for(chr in names(locs_list) ){
      dmrs <- dmr_layer_list[[layer_name]][[chr]]
      temp_locs <- locs_list[[chr]]
      if(nrow(dmrs) == 0){next}
      for(i in 1:nrow(dmrs)){
        # remove all locs contained in dmr and replace with the dmr in locs_list (e.g. a meta_loc)
        which <- which(temp_locs$start == dmrs$start[i]):which(temp_locs$stop==dmrs$stop[i])
        temp_locs[which[1],] <- data.frame("pval"=dmrs$pval_region[i], "pval_rank"=-1,"start"=dmrs$start[i],"stop"=dmrs$stop[i])
        temp_locs[which[-1],] <- NA
      }
      temp_locs <- temp_locs[complete.cases(temp_locs),]
      locs_list[[chr]] <- temp_locs
    }

    ## update pval_rank to reflect consolidation of locs into meta_locs
    temp_locs <- data.frame("pval"=numeric(),"pval_rank"=numeric(),"start"=numeric(),"stop"=numeric(),"chr"=character())
    for(chr in names(locs_list)){
      temp_temp_locs <- locs_list[[chr]]
      temp_temp_locs$chr <- chr
      temp_locs <- rbind(temp_locs, temp_temp_locs)
    }
    temp_locs$pval_rank <- rank(temp_locs$pval, ties.method = "max") ## determine rank of each p value, used for region significance
    ## update total_locs to reflect consolidation of locs into meta_locs
    total_locs <- nrow(temp_locs)

    ## update locs_list with updated pval_ranks
    locs_list <- list()
    for(chr in unique(temp_locs$chr)){
      locs_list[[chr]] <- temp_locs[which(temp_locs$chr==chr),] ## separate locs by chromosome
      locs_list[[chr]] <- locs_list[[chr]][order(locs_list[[chr]]$start),] ## order locs by position
      rownames(locs_list[[chr]]) <- NULL
    }

  }


  ## return output
  if(output_type == "simple"){
    return( dmr_layer_list[[length(dmr_layer_list)]])
  } else if(output_type == "complete"){
    return( dmr_layer_list)
  } else {return(-1)}
}

