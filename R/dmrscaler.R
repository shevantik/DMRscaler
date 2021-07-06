#' dmrscaler
#'
#'
#' @param
#' @param
#' @export
#'

dmrscaler <- function(locs, layer_sizes = c(4,8,16,32,64),  layer_step_fracs = NA, fdrscaler, cltable, signif_cutoff = "0.9999" ){
  if(is.na(layer_step_fracs)){ layer_step_fracs <- rep(2, length(layer_sizes))}
  layers<-list()
  for(i in 1:length(layer_sizes)){
    print(paste("layer", i, sep="_"))
    nn_step_fraction = layer_step_fracs[i]  ##
    nn_step_size <- max( floor(layer_sizes[i] / nn_step_fraction), 1)
    nn <- n_nearest_window_scoring_func(indat = locs, n_nearest = layer_sizes[i], step_size = nn_step_size, FDR = fdrscaler)
    signn <- determine_significant_windows(window_results=nn, indat=locs, quants=cltable , quants_significance_cutoff = signif_cutoff )
    signn <- add_significance(three_window_list =  signn, lookup_table = cltable)

    ## multiple chromosomes each a list, coerce to single dataframe
    signn <- bind_rows(signn)
    layers[[i]]<-signn
  }


  names(layers)<-paste("layer", layer_sizes, sep="_")


  atomic_layer <- locs
  for(i in 1:length(layers)){
    for(k in 1:length(layers[[i]]$start_pos)){
      layers[[i]]$start_index[k]<-which(atomic_layer$pos==layers[[i]]$start_pos[k] & as.character(atomic_layer$chr) == as.character(layers[[i]]$chr[k]))
      layers[[i]]$stop_index[k]<-which(atomic_layer$pos==layers[[i]]$stop_pos[k] & as.character(atomic_layer$chr) == as.character(layers[[i]]$chr[k]))
    }
  }
  built_layers <- list()
  built_layers[[1]] <- in_layer_merge(dmrs = layers[[1]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
  built_layers[[1]] <- in_layer_merge(dmrs = built_layers[[1]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
  built_layers[[1]] <- trim_layer(dmrs = built_layers[[1]], CG_table = atomic_layer, FDR_scaler = 2, lookup_table = cltable)
  for(i in 2:length(layers)){
    #print(i)
    built_layers[[i]] <- build_next_layer(prev_layer = built_layers[[i-1]],
                                                     windows_to_layer = layers[[i]],
                                                     CG_table = atomic_layer,
                                                     FDR_scaler=fdrscaler,
                                                     lookup_table = cltable)
    built_layers[[i]] <- in_layer_merge(dmrs = built_layers[[i]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
    temp <- in_layer_merge(dmrs = built_layers[[i]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
    while(!identical(built_layers[[i]], temp)){
      built_layers[[i]] <- temp
      temp <- in_layer_merge(dmrs = built_layers[[i]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
    }

    built_layers[[i]] <- trim_layer(dmrs = built_layers[[i]], CG_table = atomic_layer, FDR_scaler = 2, lookup_table = cltable)
    # print("done with one")
  }
  return(built_layers)

}
