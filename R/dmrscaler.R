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
                      window_sizes = c(1,2,4,8,16,32,64),
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
  stopifnot( all(window_sizes >= 1 ) )

  ## update locs to remove signal from locs with -log(p) < -log(cutoff) and set rank
  locs$pval[which(locs$pval > locs_pval_cutoff)] <- 1 ## set -log(p) to 0 if p is above cutoff
  locs$pval_rank <- rank(locs$pval, ties.method = "max") ## determine rank of each p value, used for region significance
  total_locs <- nrow(locs)

  ## organize locs into list of dataframes where each dataframe is from a unique chr
  locs_list <- list()
  for(chr in unique(locs$chr)){
    locs_list[[chr]] <- locs[which(locs$chr==chr),c("pos","pval","pval_rank")] ## separate locs by chromosome
    locs_list[[chr]] <- locs_list[[chr]][order(locs_list[[chr]]$pos),] ## order locs by position
    rownames(locs_list[[chr]]) <- NULL
  }

  ## build the primary output object i.e. dmr_layer_list
  dmr_layer_list <- foreach(chr_locs = locs_list, .final = function(x) setNames(x, names(locs_list))) %dopar% {
    ### first call features in layer using independently of all other layer
    for(window_index in 1:length(window_sizes)){
      window_size <- window_sizes[window_index]
      which_signif <- which(chr_locs$pval < locs_pval_cutoff)
      which_signif_index <- 1
      dmrs <- data.frame(start=numeric(),stop=numeric(),pval_region=numeric() )
      next_dmr <- data.frame(start=-1,stop=-1,pval_region=-1 )
      left_signif_index <- -1

      while(TRUE){
        current_signif_index <- ifelse(which_signif_index > length(which_signif), -1, which_signif[which_signif_index])
        left_signif_index <- ifelse(left_signif_index == -1, current_signif_index, left_signif_index)
        which_signif_index <- which_signif_index + 1
        right_signif_index <- ifelse(which_signif_index > length(which_signif), -1, which_signif[which_signif_index])
        ## test whether at end of array of significant locs, if yes record last significant region and break
        if(current_signif_index == -1 | right_signif_index == -1 ){
          if(nrow(next_dmr) > 0){
            dmrs <- rbind(dmrs, next_dmr)
          }
          break
        }

        ## stopping case for individual region: window does not contain any additional significant locs
        if(right_signif_index - current_signif_index > window_size){ ## extending region would not include any additional signif locs
          if(!all(next_dmr == -1)){
            dmrs <- rbind(dmrs, next_dmr)
            next_dmr <- data.frame(start=-1,stop=-1,pval_region=-1 )
          }
          left_signif_index <- -1
          next
        }

        ## update window_locs
        window_locs <- chr_locs[left_signif_index:right_signif_index,]

        ## use series of hypergeometric tests for significance
        window_loc_ranks <- window_locs$pval_rank[order(window_locs$pval_rank)][-1]  # [-1] drops most significant loc. Most sig serves as prior
        window_signif <- 1
        n <- total_locs
        for(i in length(window_loc_ranks):1 ){
          window_signif = window_signif * dhyper(x=i, m=window_loc_ranks[i], n=max(0,n-window_loc_ranks[i]), k=i)
          n <- window_loc_ranks[i]-1
        }

        if(window_signif < region_signif_cutoff){
          next_dmr$start <- chr_locs$pos[left_signif_index]
          next_dmr$stop <- chr_locs$pos[right_signif_index]
          next_dmr$pval_region <- window_signif
          next
        } else { next }


      }
    }
  }

  ## return output
  if(output_type == "simple"){
    return( dmr_layer_list[[length(dmr_layer_list)]])
  } else if(output_type == "complete"){
    return( dmr_layer_list)
  } else {return(-1)}
}


# dmrscaler <- function(locs, layer_sizes = c(4,8,16,32,64),  layer_step_fracs = NA, fdrscaler, cltable, signif_cutoff = "0.9999" ){
#
#   if(is.na(layer_step_fracs)){ layer_step_fracs <- rep(2, length(layer_sizes))}
#   layers<-list()
#   for(i in 1:length(layer_sizes)){
#     print(paste("layer", i, sep="_"))
#     nn_step_fraction = layer_step_fracs[i]  ##
#     nn_step_size <- max( floor(layer_sizes[i] / nn_step_fraction), 1)
#     nn <- n_nearest_window_scoring_func(indat = locs, n_nearest = layer_sizes[i], step_size = nn_step_size, FDR = fdrscaler)
#     signn <- determine_significant_windows(window_results=nn, indat=locs, quants=cltable , quants_significance_cutoff = signif_cutoff )
#     signn <- add_significance(three_window_list =  signn, lookup_table = cltable)
#
#     ## multiple chromosomes each a list, coerce to single dataframe
#     signn <- bind_rows(signn)
#     layers[[i]]<-signn
#   }
#
#
#   names(layers)<-paste("layer", layer_sizes, sep="_")
#
#   print("adding CG indices")
#   atomic_layer <- locs
#   for(i in 1:length(layers)){
#     for(k in 1:length(layers[[i]]$start_pos)){
#       layers[[i]]$start_index[k]<-which(atomic_layer$pos==layers[[i]]$start_pos[k] & as.character(atomic_layer$chr) == as.character(layers[[i]]$chr[k]))
#       layers[[i]]$stop_index[k]<-which(atomic_layer$pos==layers[[i]]$stop_pos[k] & as.character(atomic_layer$chr) == as.character(layers[[i]]$chr[k]))
#     }
#   }
#   print("merging layers")
#   built_layers <- list()
#   built_layers[[1]] <- in_layer_merge(dmrs = layers[[1]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
#   built_layers[[1]] <- in_layer_merge(dmrs = built_layers[[1]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
#   built_layers[[1]] <- trim_layer(dmrs = built_layers[[1]], CG_table = atomic_layer, FDR_scaler = 2, lookup_table = cltable)
#   for(i in 2:length(layers)){
#     print(paste("merging layer", i))
#     built_layers[[i]] <- build_next_layer(prev_layer = built_layers[[i-1]],
#                                                      windows_to_layer = layers[[i]],
#                                                      CG_table = atomic_layer,
#                                                      FDR_scaler=fdrscaler,
#                                                      lookup_table = cltable)
#     built_layers[[i]] <- in_layer_merge(dmrs = built_layers[[i]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
#     temp <- in_layer_merge(dmrs = built_layers[[i]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
#     while(!identical(built_layers[[i]], temp)){
#       built_layers[[i]] <- temp
#       temp <- in_layer_merge(dmrs = built_layers[[i]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
#     }
#
#     built_layers[[i]] <- trim_layer(dmrs = built_layers[[i]], CG_table = atomic_layer, FDR_scaler = 2, lookup_table = cltable)
#     # print("done with one")
#   }
#
#   for(i in 1:length(built_layers)){
#     built_layers[[i]] <- built_layers[[i]][,which(is.element(colnames(built_layers[[i]]),
#                                                                       c("chr","start_pos","stop_pos","numCG",
#                                                                         "unsigned_bin_score", "unsigned_bin_sig",
#                                                                         "start_index","stop_index")))]
#   }
#   return(built_layers)
#
# }
