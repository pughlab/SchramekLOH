#' .getAltType
#'
#' @param x
#'
#' @return
#'
#' @examples
.getAltType <- function(x) {
  alt.x <- rep("NA", length(x))
  non.idx <- unique(sort(unlist(sapply(c("fs", "del", "\\*",
                                         "splice"), grep, x = x))))
  mis.idx <- grep("^[A-Z][0-9]+[A-Z]$", x = x)
  alt.x[non.idx] <- "TRUNC"
  alt.x[mis.idx] <- "MISSENSE"
  alt.x
}

#' .cbioMutations
#'
#' @param i
#' @param genes
#'
#' @return
#'
#' @examples
.cbioMutations <- function(i, genes){
  mut.idx <- match(paste0(genes, "_MUT"), names(i))
  muts <- sapply(i[mut.idx], strsplit, split="\t")
  mut.type <- sapply(muts, .getAltType)

  mut.gene <- rep(genes, sapply(muts, length))
  muts <- as.character(unlist(muts))
  mut.type <- as.character(unlist(mut.type))

  keep.idx <- which(!mut.type=='NA')

  if(length(keep.idx) > 0){
    data.frame('Sample'=rep(as.character(i['TRACK_ID']), length(keep.idx)),
               'Gene'=as.character(mut.gene)[keep.idx],
               'mut'=as.character(muts)[keep.idx],
               'mut.type'=as.character(mut.type)[keep.idx])
  } else { NULL }
}


#' .cbioFusions
#'
#' @param i
#' @param genes
#'
#' @return
#'
#' @examples
.cbioFusions <- function(i, genes){
  fus.idx <- match(paste0(genes, "_FUSION"), names(i))
  fus.type <- rep("NA", length(fus.idx))
  fus.type[which(i[fus.idx] == 'FUSION')] <- 'FUSION'
  keep.idx <- which(!fus.type=='NA')

  if(length(keep.idx) > 0){
    data.frame('Sample'=rep(as.character(i['TRACK_ID']), length(keep.idx)),
               'Gene'=genes[keep.idx],
               'mut'=as.character(i[fus.idx])[keep.idx],
               'mut.type'=as.character(fus.type)[keep.idx])
  } else { NULL }
}

#' .cbioCNA
#'
#' @param i
#' @param genes
#' @param use.absolute
#' @param sample.stdres
#'
#' @return
#'
#' @examples
.cbioCNA <- function(i, genes, use.absolute=FALSE, sample.stdres=NULL){
  if(!use.absolute){
    cna.idx <- match(paste0(genes, "_CNA"), names(i))
    keep.idx <- which(i[cna.idx] != 'no_alteration')
    cnas <- as.character(i[cna.idx])[keep.idx]
  } else {
    id <- mapIds(as.character(i['TRACK_ID']), mapping.cov, in.type='tcga', out.type='affy')
    seg.mean <- sample.stdres[[id]][genes, 'seg.mean']
    cn <- sapply(seg.mean, function(x) switch(x,
                                              '-1'='HETLOSS',
                                              '-2'='HOMDEL',
                                              'NA'))
    cn <- as.character(cn)
    keep.idx <- which(cn != 'NA')
    cnas <- cn[keep.idx]
  }


  if(length(keep.idx) > 0){
    data.frame('Sample'=rep(as.character(i['TRACK_ID']), length(keep.idx)),
               'Gene'=genes[keep.idx],
               'mut'=cnas,
               'mut.type'=rep("CNA", length(keep.idx)))
  } else { NULL }
}

#' cbioConv
#'
#' @param i
#' @param genes
#' @param ... params passed to .cbioCNA; use.absolute and sample.stdres
#'
#' @return
#' @export
#'
#' @examples
cbioConv <- function(i, genes, ...){
  fus <- .cbioFusions(i, genes)
  mut <- .cbioMutations(i, genes)
  cna <- .cbioCNA(i, genes, ...)

  Reduce(function(x,y) rbind(x,y), list(fus, mut, cna))
}
