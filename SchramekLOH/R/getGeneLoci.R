#' getGeneLoci
#'
#' @param gene Gene ID in HUGO format, returns entrez ids
#'
#' @return
#' @export
#'
#' @examples
getGeneLoci <- function(gene){
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- genes(txdb, single.strand.genes.only=FALSE)
  map <- as.data.frame(org.Hs.egSYMBOL)

  # Map Hugo genes to entrez ids
  eid <- sapply(gene, function(i) grep(paste0("^", i, "$"), map$symbol))
  gene.id <- map[eid,]$gene_id
  # Lookup loci
  lookup <- sapply(gene.id, function(i) genes[[i]][1,])


  data.frame("gene"=gene,
             "chr"=sapply(lookup, function(i) as.character(seqnames(i))),
             "start"=sapply(lookup, function(i) start(i)),
             "end"=sapply(lookup, function(i) end(i)))
}
