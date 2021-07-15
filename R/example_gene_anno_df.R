#' example_gene_anno_df
#'
#'
#' @export
#'



example_gene_anno_df <- function(chr, start, stop){
  # get gene annotations   , 'ensembl_exon_id'
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  genes <- getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol', 'ensembl_gene_id','gene_biotype'),
                 filters = list('biotype'='protein_coding'),
                 mart = ensembl, useCache = F) ## update BiocManager to avoid need to include useCache false flag
  genes <- genes[which(is.element(genes$chromosome_name, c(1:22, "X", "Y")) & genes$hgnc_symbol != "" ) ,]
  genes$chromosome_name <- paste("chr", genes$chromosome_name, sep = "")



  which_genes <- which(genes$chromosome_name == chr & ((genes$start_position >= start & genes$start_position <= stop) | (genes$end_position >= start & genes$end_position <= stop)  ))
  which_genes_ids <- genes$ensembl_gene_id[which_genes]

  if(is.null(which_genes_ids)){ genesub <- NA
  }else{
    # attribute list without hgnc_symbol
    attributes.1 = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "ensembl_gene_id",  "ensembl_transcript_id","ensembl_exon_id")
    # attribute list with hgnc_symbol & ensembl_gene_id
    attributes.2 = c("hgnc_symbol","ensembl_gene_id")
    # get results for each query
    results.1 = getBM(attributes = attributes.1, filters = list('ensembl_gene_id'=which_genes_ids), mart = ensembl,  useCache = F)
    results.2 = getBM(attributes = attributes.2, filters = list('ensembl_gene_id'=which_genes_ids), mart = ensembl,  useCache = F)
    # merge the results for both queries
    genesub = merge(results.1,results.2,by='ensembl_gene_id',all.x=T)

    genesub$chromosome_name <- paste("chr", genesub$chromosome_name, sep = "")
    colnames(genesub)[colnames(genesub)=="chromosome_name"] <- "chromosome"
    colnames(genesub)[colnames(genesub)=="exon_chrom_start"] <- "start"
    colnames(genesub)[colnames(genesub)=="exon_chrom_end"] <- "end"
    colnames(genesub)[colnames(genesub)=="ensembl_gene_id"] <- "gene"
    colnames(genesub)[colnames(genesub)=="ensembl_transcript_id"] <- "transcript"
    colnames(genesub)[colnames(genesub)=="ensembl_exon_id"] <- "exon"
    colnames(genesub)[colnames(genesub)=="hgnc_symbol"] <- "symbol"
  }
  return(genesub)
}
