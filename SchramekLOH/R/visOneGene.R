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
visOneGene <- function(gene, rge, mut){
  ## Calculate the ranged Std.Res in context
  # Combines the StdRes of surrounding 1Mb bins to look for stable LOH regions
  range.stdres <- sapply(names(sample.stdres), aggregateStdRes,
                         range=rge, gene=gene,
                         sample.stdres=sample.stdres,
                         all.stdres=all.stdres)
  range.stdres <- data.frame(t(range.stdres), stringsAsFactors=FALSE)

  # Index all the HETLOSS samples
  hetloss.idx <- which(mut.attr[,mut]  == 'HETLOSS')
  homloss.idx <- which(mut.attr[,mut]  == 'HOMDEL')
  mut.ids <- mut.attr[c(hetloss.idx, homloss.idx), 'TRACK_ID', drop=TRUE]
  mapped.ids <- mapIds(mut.ids, mapping.cov, in.type='tcga', out.type='affy')

  # Set a boolean tag for HETLOSS annotation
  range.stdres$HETLOSS <- FALSE
  range.stdres[which(rownames(range.stdres) %in% mapped.ids),]$HETLOSS <- TRUE

  ## Visualize the plots
  with(range.stdres, plot(mean~seg, cex=sd, xlim=c(-1, 0.5),
                          col=alpha('grey', 1-rescale(sd, to=c(0,1))),
                          pch=16, ylab="StdRes", xlab="Seg"))
  axis(side=1, at = seq(-5, 4, by=0.1),
       labels = rep("", length(seq(-5, 4, by=0.1))), tick = TRUE)

  with(range.stdres[which(range.stdres$HETLOSS),],
       points(mean~seg, cex=sd,
              col=alpha('red', 1-rescale(sd, to=c(0,1))),
              pch=16))
}
