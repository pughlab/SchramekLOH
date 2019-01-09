#' generateIgvAttributes
#'
#' @param sample.stdres
#' @param idx
#' @param mapping.cov
#' @param stdres.thresh
#'
#' @return
#' @export
#'
#' @examples
generateIgvAttributes <- function(idx, sample.stdres,
                                  mapping.cov, stdres.thresh=-5){
  i <- sample.stdres[[idx]]

  # Get SEG ID
  if(names(sample.stdres)[idx] == mapping.cov$V3[idx]){
    loh.df <- data.frame("TRACK_ID" = gsub("01[AB]-.*", "01", mapping.cov$V2[idx]),
                         stringsAsFactors = FALSE)
  }

  loh.stat <- i$stdres < stdres.thresh
  loh.val <- rep("Het", nrow(i))
  loh.val[which(loh.stat)] <- 'LOH'
  names(loh.val) <- i$gene
  cbind(loh.df, t(data.frame(loh.val)))
}

#' doTheChi
#'
#' @param ctbl
#' @param attributes
#' @param mut.attr
#' @param acceptable.terms
#' @param tbl.idx
#' @param gene
#' @param mut
#'
#' @return
#' @export
#'
#' @examples
doTheChi <- function(ctbl=NULL, attributes, mut.attr,
                     acceptable.terms=c("HOMDEL", "HETLOSS", "LOH",
                                        "Het", "no_alteration", 'Unknown'),
                     tbl.idx=c(1,3), gene='ADAM10', mut='ADAM10_CNA'){
  return.tbl <- FALSE
  if(is.null(ctbl)){
    return.tbl <- TRUE

    ## Generates the contigency table of Mutation against LOH
    require(taRifx)
    loh.attr <- remove.factors(attributes)


    # Combine LOH with Mutation
    m.attr <- remove.factors(merge(mut.attr, loh.attr, by="TRACK_ID", all=TRUE))

    # Clean the data; simplifies mutations just down to "Mut"
    m.attr[is.na(m.attr)] <- 'Unknown'
    ids <- m.attr$TRACK_ID
    m.attr <- as.matrix(m.attr)
    m.attr[!m.attr %in% acceptable.terms] <- 'Mut'  # e.g. V347G
    m.attr[,1] <- ids

    ctbl <- apply(m.attr[,-1], 2, function(i){
      apply(m.attr[,-1], 2, function(j){
        table(data.frame(i,j))
      })
    })
  }

  .chi <- function(tbl, idx=c(1,3)){
    ch <- chisq.test(tbl[,idx])
    list("Contigency"=tbl,
         "p"=ch$p.value,
         "Std.Res"=ch$stdres)
  }

  chitbl <- .chi(ctbl[[gene]][[mut]], idx=tbl.idx)
  if(return.tbl) chitbl <- list("ctbl"=ctbl, "chi"=chitbl)

  chitbl
}


#' plotIdeoGene
#'
#' @param ideo.ref
#' @param sample.stdres
#' @param gene
#' @param thresh
#'
#' @return
#' @export
#'
#' @examples
plotIdeoGene <- function(ideo.ref, sample.stdres,
                         gene, thresh=-0.1){
  .getGeneGr <- function(x, gene){
    GRanges(x[gene,  c('chr', 'seg.start', 'seg.end', 'seg.mean')])
  }

  gr.list <- lapply(sample.stdres, .getGeneGr, gene=gene)
  names(gr.list) <- gsub("01[AB]-.*", "01", mapping.cov$V2)
  gr.seg.mean <- sapply(gr.list, function(x) x$seg.mean)
  gr.list <- gr.list[which(as.numeric(gr.seg.mean) < thresh)]
  gr.list.ord <- order(sapply(gr.list, width))

  plotOnIdeo(chrom=as.character(seqnames(gr.list[[1]])),
             ideoTable=ideo_hg19,
             values_GR=gr.list[gr.list.ord], value_cols="seg.mean",
             plotType="seg_tracks",
             col=alpha("blue",0.6), vertical=F)
}
