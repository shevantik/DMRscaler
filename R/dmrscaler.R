#' dmrscaler
#'
#'
#' @param locs dataframe of measured CpG loci with columns as "chr", "pos", "pval"
#' @param loc_signif_method one of "p-value" or "fdr" (false discovery rate)
#' @param loc_signif_cutoff p-value at which individual CpGs are considered significant or desired fdr is achieved for individual CpGs
#' @param region_signif_method one of "p-value" or "benjamini-yekutieli" (false discovery rate control) or "bonferroni" (family-wise error rate)
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
#'
#'
#' @export

dmrscaler <- function(locs,
                      locs_pval_cutoff = 0.05,
                      region_signif_method = c("bonferroni","benjamini-yekutieli","p-value"),
                      region_signif_cutoff = 0.01,
                      window_type = c("k_nearest", "genomic_width"),
                      window_sizes = c(2,4,8,16,32,64),
                      dmr_constraint_list = NULL,
                      output_type = c("simple", "complete")
){

  # check input parameters are valid
  if( !all(is.element(c("chr","pos","pval"), colnames(locs))) ){
    stop("ERROR: locs must include column names: \"chr\",\"pos\",\"pval\". ")
  }
  stopifnot( !is.na(pmatch(region_signif_method, c("bonferroni","benjamini-yekutieli","p-value")) ))
  region_signif_method <- match.arg(region_signif_method)
  stopifnot( !is.na(pmatch(window_type, c("k_nearest","genomic_width")) ))
  window_type <- match.arg(window_type)
  stopifnot( !is.na(pmatch(output_type, c("simple","complete")) ))
  output_type <- match.arg(output_type)
  stopifnot( locs_pval_cutoff > 0 & locs_pval_cutoff < 1 )
  stopifnot( region_signif_cutoff > 0 & region_signif_cutoff < 1 )
  stopifnot( all(window_sizes >= 2 ) )

  # update locs to remove signal from locs with -log(p) < -log(cutoff) and set rank
  locs$pval[which(locs$pval > locs_pval_cutoff)] <- 1 ## set -log(p) to 0 if p is above cutoff
  locs$pval_rank <- rank(locs$pval, ties.method = "max") ## determine rank of each p value, used for region significance
  total_locs <- nrow(locs)

  if(region_signif_method == "bonferroni"){
    region_signif_cutoff <- region_signif_cutoff/nrow(locs)
  }

  # organize locs into list of dataframes where each dataframe is from a unique chr
  locs_list <- list()
  for(chr in unique(locs$chr)){
    locs_list[[chr]] <- locs[which(locs$chr==chr),c("pos","pval","pval_rank")] ## separate locs by chromosome
    locs_list[[chr]] <- locs_list[[chr]][order(locs_list[[chr]]$pos),] ## order locs by position
    rownames(locs_list[[chr]]) <- NULL
    locs_list[[chr]]$in_dmr <- FALSE
    locs_list[[chr]]$dmr_id <- NA ## keep track of dmr membership

  }


  # build the primary output object i.e. dmr_layers_list
  dmr_layers_list <- list()  ## main output object
  for(window_index in 1:length(window_sizes)){
    window_size <- window_sizes[window_index]
    layer_name <- paste(window_size,"_loc_window_layer", sep="")


    ### BEGIN: Set up benjamini-yekutielie ###
    if(region_signif_method=="benjamini-yekutieli"){
      benjimini_yekutieli_windows <- foreach(chr_locs = locs_list,
                                             .final = function(x) setNames(x, names(locs_list)),
                                             .errorhandling = "pass" ) %dopar% {
        benjimini_yekutieli_windows <- c()
        which_signif <- which(chr_locs$pval < locs_pval_cutoff)
        which_signif_index <- 1
        ## paint all locs that are in a window that has significance < region_signif_cutoff
        while(TRUE){
          current_signif_index <- ifelse(which_signif_index > length(which_signif), -1, which_signif[which_signif_index])
          which_signif_index <- which_signif_index + 1
          ### test whether at end of array of significant locs, if yes break
          if(current_signif_index == -1 ){ break }
          right_index <- min( nrow(chr_locs), current_signif_index+window_size-1)
          right_signif_index <- max(intersect( which_signif, current_signif_index:right_index ))
          if(right_signif_index == current_signif_index){ next }
          ### use series of hypergeometric tests for significance
          window_locs <- chr_locs[current_signif_index:right_index,]
          window_locs <- window_locs[which(window_locs$pval < 1),] ## only retain significant locs, if pval==1, loop will multiple window_signif by 1 (doing nothing)
          #### if all dmrs already belong to the same dmr, nothing to do
          if( length(unique(window_locs$dmr_id))==1 & !all(is.na(window_locs$dmr_id)) ){ next }
          #### build groups to sequentially mask
          if(any(is.na(window_locs$dmr_id))){
            which <- which.min( window_locs$pval[which(is.na(window_locs$dmr_id) )] )
            if(all(is.na(window_locs$dmr_id))){
              window_locs$dmr_id[which(is.na(window_locs$dmr_id) )][which] <- 1
            } else {
              window_locs$dmr_id[which(is.na(window_locs$dmr_id) )][which] <- max(window_locs$dmr_id)+1
            }
          }
          #### sequentially mask groups
          window_all_signif <- TRUE
          for(mask_id in unique(window_locs$dmr_id[!is.na(window_locs$dmr_id)]) ){
            which_mask <- which(window_locs$dmr_id==mask_id)
            window_loc_ranks <- window_locs$pval_rank[-which_mask]
            window_loc_ranks <- window_loc_ranks[order(window_loc_ranks)]
            window_signif <- 1
            n <- total_locs
            k <- length(current_signif_index:right_index)-length(which_mask)  ##  window_size minus length mask group
            for(i in length(window_loc_ranks):1 ){
              window_signif = window_signif * phyper(q=i-0.01, m=window_loc_ranks[i], n=max(0,n-window_loc_ranks[i]), k=k, lower.tail=F)
              n <- window_loc_ranks[i]-1
              k <- i-1
            }

            benjimini_yekutieli_windows <- c(benjimini_yekutieli_windows, window_signif)
          }
        }
        benjimini_yekutieli_windows
      }

      temp <- unname(unlist(benjimini_yekutieli_windows))
      temp <- temp[order(temp)]
      temp_df <- data.frame("benjamini-yekutieli"=temp,"rank"=1:length(temp))
      temp_df$crit <- (region_signif_cutoff * temp_df$rank) / (length(which(locs$pval<1)) * sum(1/(1:length(which(locs$pval<1)))))
      temp_df$bh_lt_crit <- temp_df$benjamini.yekutieli < temp_df$crit
      if(length(which(temp_df$bh_lt_crit))==0){
        benj_hoch_signif_cutoff <- 0
      } else {
        benj_hoch_signif_cutoff <- temp_df$benjamini.yekutieli[max(which(temp_df$bh_lt_crit))]
      }
    }
    ### END: Set up benjamini-yekutielie ###

    dmr_layers_list[[layer_name]] <- foreach(chr_locs = locs_list, .final = function(x) setNames(x, names(locs_list))) %dopar% {
      ## first call features in layer using independently of all other layer
      which_signif <- which(chr_locs$pval < locs_pval_cutoff)
      which_signif_index <- 1


      ## paint all locs that are in a window that has significance < region_signif_cutoff
      while(TRUE){
        current_signif_index <- ifelse(which_signif_index > length(which_signif), -1, which_signif[which_signif_index])
        which_signif_index <- which_signif_index + 1
        ### test whether at end of array of significant locs, if yes break
        if(current_signif_index == -1 ){ break }

        right_index <- min( nrow(chr_locs), current_signif_index+window_size-1)
        right_signif_index <- max(intersect( which_signif, current_signif_index:right_index ))
        if(right_signif_index == current_signif_index){ next }

        ### use series of hypergeometric tests for significance
        window_locs <- chr_locs[current_signif_index:right_index,]
        window_locs <- window_locs[which(window_locs$pval < 1),] ## only retain significant locs, if pval==1, loop will multiple window_signif by 1 (doing nothing)
        #### if all dmrs already belong to the same dmr, nothing to do
        if( length(unique(window_locs$dmr_id))==1 & !all(is.na(window_locs$dmr_id)) ){ next }
        #### build groups to sequentially mask
        if(any(is.na(window_locs$dmr_id))){
          which <- which.min( window_locs$pval[which(is.na(window_locs$dmr_id) )] )
          if(all(is.na(window_locs$dmr_id))){
            window_locs$dmr_id[which(is.na(window_locs$dmr_id) )][which] <- 1
          } else {
            window_locs$dmr_id[which(is.na(window_locs$dmr_id) )][which] <- max(window_locs$dmr_id)+1
          }
        }

        #### sequentially mask groups
        window_all_signif <- TRUE
        for(mask_id in unique(window_locs$dmr_id[!is.na(window_locs$dmr_id)]) ){
          which_mask <- which(window_locs$dmr_id==mask_id)
          window_loc_ranks <- window_locs$pval_rank[-which_mask]
          window_loc_ranks <-  window_loc_ranks[order(window_loc_ranks)]
          window_signif <- 1
          n <- total_locs
          k <- length(current_signif_index:right_index)-length(which_mask)  ##  window_size minus length mask group
          for(i in length(window_loc_ranks):1 ){
            window_signif = window_signif * phyper(q=i-0.01, m=window_loc_ranks[i], n=max(0,n-window_loc_ranks[i]), k=k, lower.tail=F)
            n <- window_loc_ranks[i]-1
            k <- i-1
          }
          if(region_signif_method=="benjamini-yekutieli"){
            if(window_signif > benj_hoch_signif_cutoff){
              window_all_signif <- FALSE
              break ## if window_signif > region_signif_cutoff for any masked group, region is not altered in current layer
            }
          } else {
            if(window_signif > region_signif_cutoff){
              window_all_signif <- FALSE
              break ## if window_signif > region_signif_cutoff for any masked group, region is not altered in current layer
            }
          }
        }
        if(window_all_signif){
          chr_locs$in_dmr[current_signif_index:right_signif_index] <- TRUE
        }
      } # end paint all in dmrs

      ### build dmrs
      #### add dmrs ranges
      dmrs <- data.frame(start=numeric(),stop=numeric(),pval_region=numeric(), dmr_id=numeric() )
      dmr_counter <- 1
      next_dmr <- data.frame(start=-1,stop=-1,pval_region=-1,dmr_id=dmr_counter)
      for(i in 1:nrow(chr_locs)){
        if(chr_locs$in_dmr[i]){
          if(next_dmr$start == -1){
            next_dmr$start <- chr_locs$pos[i]
          }
          if( i+1 > nrow(chr_locs) | !chr_locs$in_dmr[i+1] ){
            next_dmr$stop <- chr_locs$pos[i]
            dmrs <- rbind(dmrs, next_dmr)
            dmr_counter <- dmr_counter + 1
            next_dmr <- data.frame(start=-1,stop=-1,pval_region=-1,dmr_id=dmr_counter )
          }
        }
      }
      #### add dmrs significance
      if(nrow(dmrs) > 0 ){ ## need to test whether any dmrs were found
        for(i in 1:nrow(dmrs)){
          window_locs <- chr_locs[which(chr_locs$pos==dmrs$start[i]):which(chr_locs$pos==dmrs$stop[i]),]
          k <- nrow(window_locs)
          window_locs <- window_locs[which(window_locs$pval < 1),]
          window_loc_ranks <- window_locs$pval_rank[order(window_locs$pval_rank)]
          window_signif <- 1
          n <- total_locs
          for(j in length(window_loc_ranks):1 ){
            window_signif = window_signif *phyper(q=j-0.01, m=window_loc_ranks[j], n=max(0,n-window_loc_ranks[j]), k=k, lower.tail=F)
            n <- window_loc_ranks[j]-1
            k <- j-1
          }
          dmrs$pval_region[i] <- window_signif
          ### ADD ADJUSTED PVALUE AS  WELL
        }
      }
      dmrs
    } ## end foreach


    ## add unique dmr names to locs in dmr - will to form groups to mask
    for(chr in names(locs_list) ){
      dmrs <- dmr_layers_list[[layer_name]][[chr]]
      temp_locs <- locs_list[[chr]]
      if(nrow(dmrs) == 0){next}
      for(i in 1:nrow(dmrs)){
        # remove all locs contained in dmr and replace with the dmr in locs_list (e.g. a meta_loc)
        which <- which(temp_locs$pos == dmrs$start[i]):which(temp_locs$pos==dmrs$stop[i])
        temp_locs[which,]$in_dmr <- TRUE
        temp_locs[which,]$dmr_id <- i
      }
      locs_list[[chr]] <- temp_locs
    }
  } # end for window_size loop


  for(i in 1:length(dmr_layers_list)){
    temp <- data.frame(chr=character(), dmr_layers_list[[i]][[1]][0,])
    for(j in 1:length(dmr_layers_list[[i]])){
      if(nrow(dmr_layers_list[[i]][[j]])>0){
        temp <- rbind(temp, data.frame(chr=names(dmr_layers_list[[i]])[j], dmr_layers_list[[i]][[j]]) )
      }
    }
    dmr_layers_list[[i]] <- temp
  }

  # return output
  if(output_type == "simple"){
    return( dmr_layers_list[[length(dmr_layers_list)]])
  } else if(output_type == "complete"){
    return( dmr_layers_list)
  } else {return(-1)}
}

