#' .reduceId
#'
#' @param i
#' @param mut.i
#'
#' @return
#'
#' @examples
.reduceId <- function(i, mut.i){
  j <- paste(i[mut.i][which(i[mut.i] != 'NA')], collapse="-")
  if(j == '') j <- 'NA'
  j
}

#' .simplifyMutInfo
#'
#' @param mut.attr
#' @param mut.x
#' @param mut.y
#'
#' @return
#'
#' @examples
.simplifyMutInfo <- function(mut.attr, mut.x, mut.y){
  mut.tmp <- mut.attr[,c(mut.x, mut.y)]
  mut.tmp[mut.tmp == 'no_alteration'] <- 'NA'
  mut.tmp <- as.matrix(mut.tmp)
  mut.tmp[!mut.tmp %in% c('HETLOSS', 'HOMDEL', 'NA', 'GAIN')] <- 'SNV'
  mut.tmp <- data.frame(mut.tmp, stringsAsFactors = FALSE)
  mut.tmp$TCGA_ID <- mut.attr$TRACK_ID
  mut.tmp
}

#' compTwoGenes
#'
#' @param gene.x
#' @param gene.y
#' @param col.df
#' @param cex.type
#'
#' @return
#' @export
#'
#' @examples
compTwoGenes <- function(gene.x, gene.y, col.df, cex.type='xy', plot.legend=TRUE){
  ## Obtain the Seg + LOH info of both Gene X and Gene Y
  # Gene.x: AJUBA
  range.stdres.x <- sapply(names(sample.stdres), aggregateStdRes,
                           range=rge, gene=gene.x, sample.stdres=sample.stdres, all.stdres=all.stdres)
  range.stdres.x <- data.frame(t(range.stdres.x), stringsAsFactors=FALSE)
  range.stdres.x$unity <- apply(range.stdres.x, 1, function(x) mean(c(x['mean'], x['seg'])))
  range.stdres.x$TRACK_ID <- rownames(range.stdres.x)

  # Gene.y: ADAM10
  range.stdres.y <- sapply(names(sample.stdres), aggregateStdRes,
                           range=rge, gene=gene.y, sample.stdres=sample.stdres, all.stdres=all.stdres)
  range.stdres.y <- data.frame(t(range.stdres.y), stringsAsFactors=FALSE)
  range.stdres.y$unity <- apply(range.stdres.y, 1, function(x) mean(c(x['mean'], x['seg'])))
  range.stdres.y$TRACK_ID <- rownames(range.stdres.y)

  range.stdres.xy <- merge(range.stdres.x, range.stdres.y, by='TRACK_ID', all=TRUE)


  ## Obtain the mutation information both Gene X and Gene Y
  mut.x <- colnames(mut.attr)[grep(gene.x, colnames(mut.attr))]
  mut.y <- colnames(mut.attr)[grep(gene.y, colnames(mut.attr))]

  # Format the two gene mutations into one dataframe and simplify
  mut.tmp <- .simplifyMutInfo(mut.attr, mut.x, mut.y)
  x.ids <- apply(mut.tmp, 1, .reduceId, mut.i=mut.x)
  y.ids <- apply(mut.tmp, 1, .reduceId, mut.i=mut.y)
  mut.tmp$UID <- paste(x.ids, y.ids, sep="_")
  mut.tmp$UID <- gsub("NA-NA", "NA", mut.tmp$UID)

  ## Map the TCGA IDs to the LOH/Seg data frame
  mapped <- mapIds(range.stdres.xy$TRACK_ID, mapping.cov, in.type='affy', out.type='tcga')
  mapped <- gsub("01[AB]-.*", "01", mapped)
  range.stdres.xy$TCGA_ID <- mapped

  # Map the Unique Mutation ID to the LOH/Seg data frame
  range.stdres.xyz <- merge(range.stdres.xy, mut.tmp[,c("TCGA_ID", "UID")], by="TCGA_ID", all.x=TRUE)

  # Map the Unique Mutation ID colour schema
  range.stdres.xyz <- merge(range.stdres.xyz, col.df, by="UID", all.x=TRUE)

  ## Visualization
  #LEGEND
  if(plot.legend){
    plot(0, type='n', ylim=c(0, nrow(col.df)), xlim=c(0,5), xaxt='n', yaxt='n', xlab='', ylab='')
    for(i in seq_along(col.df$UID)){
      rect(xleft=0, ybottom=(i-1), xright = 1, ytop=i, col=alpha(col.df$col[i], 0.7))
      text(x=1.5, y=i-0.5, labels = paste(col.df$UID[i], col.df$col[i], sep="-"), adj=0)
    }
  }
  #SEG-LOH Plot
  cex.xy <- switch(cex.type,
                   "xy"=with(range.stdres.xyz, rescale(abs(mean.x * mean.y), to=c(1,4))),
                   "x"=with(range.stdres.xyz, rescale(abs(mean.x), to=c(1,4))),
                   "y"=with(range.stdres.xyz, rescale(abs(mean.y), to=c(1,4))))
  with(range.stdres.xyz, plot(seg.x, seg.y, col=alpha(col, 0.7), pch=16,
                              xlim=c(-1,0.5), ylim=c(-1, 0.5),
                              xlab=paste(gene.x, "seg"), ylab=paste(gene.y, "seg"),
                              cex=cex.xy))
  range.stdres.xyz
}

#' plotSegSizes
#'
#' @return
#' @export
#'
#' @examples
plotSegSizes <- function(){
  range.stdres.xyz$width.x <- with(range.stdres.xyz, seg.end.x - seg.start.x)
  range.stdres.xyz$width.y <- with(range.stdres.xyz, seg.end.y - seg.start.y)
  mb <- 1000000
  with(range.stdres.xyz[which(range.stdres.xyz$UID == 'HETLOSS_HETLOSS'),],
       plot((width.x / mb), (width.y / mb), pch=16,
            cex=rescale(abs(seg.x * seg.y), to=c(1,5)), col=alpha(col, 0.8),
            xlab=paste0(gene.x, " HetLoss Seg size"), ylab=paste0(gene.y, " HetLoss Seg size"),
            sub="Circle size proportional to purity"))
}
