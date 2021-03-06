#' combineSegLohPurity
#'
#' @param i
#' @param seg.ids
#' @param all.stdres
#'
#' @return
#' @export
#'
#' @examples
combineSegLohPurity <- function(i, seg.ids, all.stdres){
  # Map the IDs
  affy.id <- i
  tcga.id <- mapIds(i, mapping.cov, in.type='affy', out.type='tcga')

  if(any(tcga.id == names(seg.ids)) & any(affy.id == names(all.stdres))){
    # Obtain GRanges objects for the CN segments and StdRes estimations
    stdres.gr <- sort(makeGRangesFromDataFrame(all.stdres[[affy.id]], keep.extra.columns=TRUE))
    seg.gr <- sort(makeGRangesFromDataFrame(seg.ids[[tcga.id]], keep.extra.columns = TRUE))

    # Overlap the two GRanges object and make them comparable
    ov <- findOverlaps(seg.gr, stdres.gr)
    s.seg.df <- data.frame(seg.gr[queryHits(ov), c('ID', 'sample', 'seg.mean', 'copy_ratio', 'subclonal')],
                           stringsAsFactors=FALSE)
    s.stdres.df <- data.frame(stdres.gr[subjectHits(ov), c('stdres')], stringsAsFactors=FALSE)

    # Combine and add purity
    s.df <- cbind(s.seg.df, s.stdres.df)
    s.df$purity <- purity.ids[[tcga.id]]$purity[1]
    s.df
  } else {
    NULL
  }
}

#' getLohProb
#'
#' @param s sample specific dataframe, split from combineSegLohPurity()
#' @param mod logistic regression model containing predictors: 'mean', 'purity'
#' @param pos
#'
#' @return
#' @export
#'
#' @examples
getLohProb <- function(s, mod, pos=c("seqnames", "start", "end", "width")){
  dedup.sdf <- s

  # Remove duplicate indices arbitrarily
  dup.idx <- which(duplicated(dedup.sdf[,pos]))
  dedup.sdf <- dedup.sdf[-dup.idx,]

  # Calculates the LOH probabilities
  predict(mod, dedup.sdf, type='response')
}


#' plotLohProb
#'
#' @param goi Gene of interest (e.g. ADAM10)
#' @param i if diff=TRUE; returned smoothed LOESS from plotLohProb()
#' @param j if diff=TRUE; returned smoothed LOESS from plotLohProb()
#' @param diff compares two LOESS segments
#' @param main Title
#'
#' @return
#' @export
#'
#' @examples
plotLohProb <- function(goi=NULL, soi=NULL, i=NULL, j=NULL, diff=FALSE,
                        main=NULL, is.expr=FALSE, xlab='Probability',
                        minx=NULL){
  if(is.null(goi) & is.null(soi)){
    soi <- as.character(cbio.abs$Sample)
    soi <- soi[soi %in% colnames(loh.mat)]
    goi <- 'all'
  } else if (is.null(soi)) {
    soi <- as.character(cbio.abs[which(cbio.abs$Gene == goi),]$Sample)
    main <- goi
  }

  if(!is.expr){
    print("Processing LOH...")
    loh <- data.frame("prob"=apply(loh.mat[,soi], 1, mean),
                      "idx"=c(1:nrow(loh.mat)))
    rle.x <- rle(as.character(loh.pos$seqnames))
    na.chr <- which(!is.na(rle.x$values))
    chr.xpos <- c(1, cumsum(rle.x$length))
  } else {
    print("Processing expression...")
    soi <- gsub("-", ".", soi)

    gaf.spl <- split(gaf.ord, f=gaf.ord$chr)
    cum.x <- sapply(gaf.spl, function(x) c(min(x$start), diff(x$start)))
    cum.x <- cumsum(as.numeric(unlist(cum.x[paste0("chr", c(1:22, "X", "Y"))])))
    length(cum.x) <- nrow(gaf.ord)
    gaf.ord$cumx <- cum.x

    loh <- data.frame("prob"=apply(z.ex[,soi], 1, function(x) mean(as.numeric(x), na.rm=TRUE)),
                      "idx"=gaf.ord$cumx)

    rle.x <- rle(as.character(gaf.ord$chr))
    na.chr <- which(!is.na(rle.x$values))
    chg <- c(1, cumsum(rle.x$length))
    chr.xpos <- na.omit(gaf.ord$cumx[chg])
  }



  # Fit the Loess Model
  loessMod10 <- loess(prob ~ idx, data=loh, span=0.03) # 10% smoothing span
  smoothed10 <- predict(loessMod10, loh$idx)

  # Plot the loess model
  if(diff){
    smoothed10 <- (i - min(i, na.rm=TRUE)) - (j - min(j, na.rm=TRUE))
    min.x <- min(smoothed10, na.rm=TRUE) - 0.03
    max.x <- max(smoothed10, na.rm=TRUE) + 0.03
  } else {
    min.x <- min(loh$prob, na.rm=TRUE) - 0.03
    max.x <- max(loh$prob, na.rm=TRUE) + 0.03
  }

  with(loh, plot((prob), -1*idx, pch=16, col=alpha("grey", 0.3),
                 xlim=if(!is.null(minx)) minx else c(min.x, max.x), las=1,
                 main=main, yaxt='n', ylab='', xlab=xlab,
                 type=if(diff) 'n' else 'p'))
  rect(ytop=-1*chr.xpos[-length(chr.xpos)][c(TRUE,FALSE)],  xleft = -10,
       ybottom = -1*chr.xpos[-1][c(TRUE,FALSE)], xright = 10,
       border=FALSE, col=alpha("grey", 0.3))
  lines(smoothed10, -1*(loh$idx), col="blue")


  # Label chromosomes
  #abline(v=chr.xpos, lty=2, col="grey")
  if(diff){
    abline(v=quantile(smoothed10, 0.001), lty=2, col="black")
  } else {
    if(!is.expr){
      abline(v=quantile(loh$prob, 0.1, na.rm=TRUE), lty=2, col="black")
    } else {
      abline(v=0, lty=2, col="black")
    }
  }
  text(y=(-1*(chr.xpos[-length(chr.xpos)] + 10)),
       x=rep(c(min.x, min.x+0.01), length(na.chr)),
       labels=gsub("^chr", "", rle.x$values[na.chr]),
       adj=0, cex=0.5)
  smoothed10
}
