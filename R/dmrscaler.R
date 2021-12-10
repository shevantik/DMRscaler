#' dmrscaler
#'
#'
#' @param locs dataframe of measured CpG loci with columns as "chr", "pos", "pval"
#' @param loc_signif_method one of "p-value" or "fdr" (false discovery rate)
#' @param loc_signif_cutoff p-value at which individual CpGs are considered significant or desired fdr is achieved for individual CpGs
#' @param region_signif_method one of "p-value" or "fdr" (false discovery rate) or "fwer" (family-wise error rate)
#' @param region_signif_cutoff p-value or fwer at which regions are considered significant or desired fdr is achieved for regions
#' @param layer_type specifies  one of "k_nearest" or "genomic_width"
#' @param layer_sizes vector of window size for each layer where layer_type determines whether this represents the genomic_width of windows or the k_nearest neighbors that compose a window
#' @param dmr_constraint_list (ADD LATER)(optional) a list of rules that further constrains the definition of a differentially methylated region
#' @export
#'



dmrscaler <- function(locs,
                      loc_signif_method = c("fdr","p-value"),
                      loc_signif_cutoff = 0.10,
                      region_signif_method = c("fwer","fdr","p-value"),
                      region_signif_cutoff = 0.05,
                      layer_type = c("k_nearest", "genomic_width"),
                      layer_sizes = c(1,2,4,8,16,32,64),
                      dmr_constraint_list = NULL
                      ){

  ### check input parameters are valid
  if( !all(is.element(c("chr","pos","pval"), colnames(locs))) ){
    stop("ERROR: locs must include column names: \"chr\",\"pos\",\"pval\". ")
  }
  stopifnot( !is.na(pmatch(loc_signif_method, c("fdr","p-value") )) )
  loc_signif_method <- match.arg(loc_signif_method)
  stopifnot( !is.na(pmatch(region_signif_method, c("fwer","fdr","p-value")) ))
  region_signif_method <- match.arg(region_signif_method)
  stopifnot( !is.na(pmatch(layer_type, c("k_nearest","genomic_width")) ))
  layer_type <- match.arg(layer_type)
  stopifnot( loc_signif_cutoff > 0 & loc_signif_cutoff < 1 )
  stopifnot( region_signif_cutoff > 0 & region_signif_cutoff < 1 )
  stopifnot( all(layer_sizes > 0) )

  ##

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
