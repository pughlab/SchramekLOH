#' visOneGene
#'
#' @param gene
#' @param rge
#' @param mut
#'
#' @return
#' @export
#'
#' @examples
visOneGene <- function(gene, rge, mut, use.absolute=FALSE, purity=NULL,
                       plot.cr=FALSE, use.seg=FALSE){
  ## Calculate the ranged Std.Res in context
  # Combines the StdRes of surrounding 1Mb bins to look for stable LOH regions
  range.stdres <- sapply(names(sample.stdres), aggregateStdRes,
                         range=rge, gene=gene,
                         sample.stdres=sample.stdres,
                         all.stdres=all.stdres,
                         use.absolute=use.absolute,
                         use.seg=use.seg)
  range.stdres <- data.frame(t(range.stdres), stringsAsFactors=FALSE)

  if(use.absolute & (!is.null(purity))){
    # Set a boolean tag for HETLOSS annotation
    range.stdres$HETLOSS <- FALSE
    range.stdres[which(range.stdres$seg < 0),]$HETLOSS <- TRUE

    range.stdres$GAIN <- FALSE
    range.stdres[which(range.stdres$seg > 0),]$GAIN <- TRUE

    range.stdres$sample <- rownames(range.stdres)
    range.stdres <- merge(range.stdres, purity[,c('sample', 'purity')], by='sample', all.x=TRUE)
    range.stdres <- range.stdres[,-grep("sd", colnames(range.stdres))]
    range.stdres$sd <- range.stdres$purity
    range.stdres$sd <- round(rescale((1-range.stdres$sd), to=c(1,5)), 0)

    if(plot.cr) {
      range.stdres$cn <- range.stdres$seg
      range.stdres$seg <- range.stdres$copy_ratio
    }
  } else {
    # Index all the HETLOSS samples
    hetloss.idx <- which(mut.attr[,mut]  == 'HETLOSS')
    homloss.idx <- which(mut.attr[,mut]  == 'HOMDEL')
    mut.ids <- mut.attr[c(hetloss.idx, homloss.idx), 'TRACK_ID', drop=TRUE]
    mapped.ids <- mapIds(mut.ids, mapping.cov, in.type='tcga', out.type='affy')

    # Set a boolean tag for HETLOSS annotation
    range.stdres$HETLOSS <- FALSE
    range.stdres[which(rownames(range.stdres) %in% mapped.ids),]$HETLOSS <- TRUE

    range.stdres$cn <- range.stdres$seg
  }


  nan.idx <- which(is.nan(range.stdres$mean))
  if(length(nan.idx) > 0) range.stdres <- range.stdres[-nan.idx, ]

  ## Visualize the plots
  with(range.stdres, plot(seg~mean, cex=sd,
                          ylim=if(plot.cr) c(0,1) else c(-1, 0.5),
                          col=alpha('grey', 1-rescale(sd, to=c(0,1))),
                          pch=16, xlab="StdRes", ylab="Seg"))
  axis(side=2, at = seq(-5, 4, by=0.1),
       labels = rep("", length(seq(-5, 4, by=0.1))), tick = TRUE)

  with(range.stdres[which(range.stdres$HETLOSS),],
       points(seg~mean, cex=sd,
              col=alpha('red', 1-rescale(sd, to=c(0,1))),
              pch=16))

  range.stdres
}
