#' example_generate_dmr_tree
#'
#'
#' @param
#' @param
#' @export
#'

example_generate_dmr_tree <- function(dmrscaler_result, layer, chr, start, stop ){

  strip_grs <- function(node){
    node <- list("name" = node$name, "children" = node$children)
    if(length(node$children)==0){return(node)}
    for(i in 1:length(node$children)){
      node$children[[i]]<-strip_grs(node$children[[i]])
    }
    return(node)
  }

  x_full_node_set <- list()
  for(j in 1:length(dmrscaler_result)){
    if(nrow(dmrscaler_result[[j]])==0){next;}
    x_which <- which(dmrscaler_result[[j]]$chr==chr & dmrscaler_result[[j]]$start>=start & dmrscaler_result[[j]]$stop<=stop )
    x_base <- GenomicRanges::GRanges(seqnames = dmrscaler_result[[j]]$chr[x_which],
                                     IRanges(start = dmrscaler_result[[j]]$start[x_which], end = dmrscaler_result[[j]]$stop[x_which]  ))
    x_base_nodes <- list()
    for(i in 1:length(x_base)){
      temp_name <- paste(as.character(seqnames(x_base[i,])), start(x_base[i,]), end(x_base[i,]), sep = "_")
      temp_child <- list("name"= temp_name, "alt_name"= temp_name, "grange"=x_base[i,], "children"=list() )
      x_base_nodes[[i]] <-  list("name"= temp_name, "alt_name"= temp_name, "grange"=x_base[i,], "children"=list())
    }
    x_full_node_set[[j]]<-x_base_nodes
  }
  x_full_node_set[sapply(x_full_node_set, is.null)] <- NULL

  for(j in 2:length(x_full_node_set)){
    #print(j)

    lower <- x_full_node_set[[j-1]]
    upper <- x_full_node_set[[j]]

    lower_gr <- lower[[1]]$grange
    if(length(lower) < 2){next}
    for(i in 2:length(lower)){
      lower_gr <- c(lower_gr, lower[[i]]$grange)
    }
    upper_gr <- upper[[1]]$grange
    if(length(upper) >= 2){
      for(i in 2:length(upper)){
        upper_gr <- c(upper_gr, upper[[i]]$grange)
      }
    }

    hits <- GenomicRanges::findOverlaps(lower_gr,upper_gr,type = "within")
    for( i in 1:length(hits@from)){
      upper[[ hits@to[i] ]]$children <- c(upper[[ hits@to[i] ]]$children, list(lower[[ hits@from[i] ]]))
      if(upper[[ hits@to[i] ]]$name == lower[[ hits@from[i] ]]$name |
         upper[[ hits@to[i] ]]$alt_name == lower[[ hits@from[i] ]]$name  |
         upper[[ hits@to[i] ]]$name == lower[[ hits@from[i] ]]$alt_name
      ){  ### avoid redundant child-parent pairs
        upper[[ hits@to[i] ]]$name <- "..."
      }
    }
    x_full_node_set[[j]]<-upper
  }

  x_dmr_tree <- list("name"="root","children"=x_full_node_set[[length(x_full_node_set)]])

  x_dmr_tree_clean <- list()

  node <- x_dmr_tree




  x_dmr_tree_clean <- strip_grs(x_dmr_tree)

  return(x_dmr_tree_clean)

}
