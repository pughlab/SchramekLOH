#' .intersectLohWithSeg
#'
#' @param i
#' @param chr
#' @param bin.left
#' @param bin.right
#' @param seg
#' @param stdres
#'
#' @return
#' @export
#'
#' @examples
.intersectLohWithSeg <- function(i, chr, bin.left, bin.right, seg, stdres){
  gr.i <- GRanges(i)
  gr.bin <- GRanges(seqnames = chr, ranges=IRanges(bin.left, bin.right))
  gr.seg.reduce <- cleanSeg(seg)
  gr.seg <- GRanges(seg)

  overlaps <- findOverlaps(gr.i, gr.bin)
  ov.idx <- subjectHits(overlaps)

  red.overlaps <- findOverlaps(gr.i, gr.seg.reduce)
  red.overlaps <- filterOverlaps(red.overlaps)
  red.ov.idx <- subjectHits(red.overlaps)

  seg.overlaps <- findOverlaps(gr.i, gr.seg)
  seg.overlaps <- filterOverlaps(seg.overlaps)
  seg.ov.idx <- subjectHits(seg.overlaps)

  list("genes"=data.frame("gene"=i$gene,
                          "chr"=i$chr,
                          "gene.start"=i$start,
                          "gene.end"=i$end,
                          "bin.start"=if(length(ov.idx) == 0) 0 else start(gr.bin[ov.idx]),
                          "bin.end"=if(length(ov.idx) == 0) 1 else end(gr.bin[ov.idx]),
                          "seg.start"=if(length(red.ov.idx) == 0) 0 else start(gr.seg.reduce[red.ov.idx]),
                          "seg.end"=if(length(red.ov.idx) == 0) 1 else end(gr.seg.reduce[red.ov.idx]),
                          "seg.mean"=if(length(seg.ov.idx) == 0) 0 else gr.seg[seg.ov.idx]$seg.mean,
                          "stdres"=if(length(ov.idx) == 0) 0 else stdres[ov.idx]),
       "stdres"=SchramekLOH:::.getStdResDf(chr, bin.left, bin.right,
                                           stdres))
}

#' .getStdResDf
#'
#' @param chr Chr ID
#' @param bin.left Left loci of each genomic bin
#' @param bin.right Right loci of each genomic bin
#' @param stdres The raw Standardized residuals of each genomic bin for AB
#'
#' @return
#'
#' @examples
.getStdResDf <- function(chr, bin.left, bin.right, stdres){
  gr.bin <- GRanges(seqnames = chr, ranges=IRanges(bin.left, bin.right))
  qval.df <- data.frame("chrom"=as.character(seqnames(gr.bin)),
                        "start"=start(gr.bin),
                        "end"=end(gr.bin),
                        "stdres"=stdres)
  qval.df
}

#' .chiAB
#'
#' @param cnts
#'
#' @return
#'
#' @examples
.chiAB <- function(cnts){
  qval <- apply(cnts, 1, function(i){
    p <- tryCatch({ suppressWarnings(chisq.test(i[c(1,3,5)])$stdres[2]) },
                  error=function(e){NA})
    p
  })
  qval.raw <- qval
  qval[qval <= -10] <- -10

  list("raw"=qval.raw,
       "thr"=qval)
}

#' .mapId
#'
#' @param id
#' @param mapping
#' @param use.affy
#'
#' @return
#'
#' @examples
.mapId <- function(id, mapping, use.affy=FALSE){
  tcga.id <- mapping[match(id, mapping$V3), 'V2']
  tcga.id <- gsub("01[AB]-.*", "01", tcga.id)
  pdf.id <- tcga.id
  if(use.affy){
    tcga.id <- id
  }

  list("id"=id, "tcga"=tcga.id, "pdf"=pdf.id)
}

#' .formatExprZ
#'
#' @param z A vector of z-values for gene expression
#'
#' @return
#'
#' @examples
.formatExprZ <- function(z){
  z.to.y <- z
  z.to.y[z.to.y > 2] <- 2  # Upper-thresh of expression z-score (Viz. only)
  z.to.y[z.to.y < -2] <- -2 # Lower-thresh of expression z-score (Viz. only)
  range <- c(-2.1, -1.1)  # Range of y-values to plot on to
  z.to.y <- rescale(z.to.y, to=range, from=c(-2, 2))
  mid.y <- (c(range[1] - range[2]) /2) + (range[2])
  z.col <- rep("red", length(z.to.y))  # Set high expr as red
  z.col[z.to.y < mid.y] <- 'blue'  # and low expr as blue

  list("col"=z.col,
       "z0"=mid.y,
       "z"=z.to.y)
}

#' mapAndPlotFeatures
#'
#' @param id
#' @param mapping.cov
#' @param use.affy
#' @param plotsdir
#' @param snp6.chr
#' @param seg.chr
#' @param bs.chr
#' @param z.chr
#' @param r
#' @param p
#' @param cn
#' @param gen.plot
#'
#' @return
#' @export
#'
#' @examples
mapAndPlotFeatures <- function(id, mapping.cov, use.affy, plotsdir,
                               snp6.chr, seg.chr, bs.chr, z.chr,
                               r, p, cn, gen.plot=TRUE){
  ids <- SchramekLOH:::.mapId(id, mapping.cov, use.affy=use.affy)

  if(gen.plot){
    pdf(file.path(plotsdir, paste0(ids[['pdf']], ".", use.affy, ".pdf")), height=50)
    split.screen(c(length(names(snp6.chr)), 1))
  }
  seg.id <- seg.chr[[ids[['tcga']]]]


  gene.stdres <- lapply(names(snp6.chr), function(chr){
    if(gen.plot){
      screen(match(chr, names(snp6.chr)))
      par(mar=c(3.1, 4.5, 0.5, 2.1))
    }

    df <- data.frame(x = snp6.chr[[chr]]$start,
                     y = bs.chr[[chr]][,id])
    chr.bins <- round(max(snp6.chr[[chr]]$start) / bin.size, 0)

    ## Visualize the density of AA, AB, BB SNPs in a 1Mb bin
    h2 <- hist2d(df, nbins=c(chr.bins, 5), col=r, ylim=c(-2.2,2),
                 xlim=c(0, max(snp6.chr[[chr]]$start)),
                 yaxt='n', xaxt='n', show=gen.plot,
                 xlab = "Genomic Loci (Mb)", ylab="")

    ## Labels, titles, formatting, etc.
    if(gen.plot){
      title(ylab=chr, line = 3.5)
      spacer <- c(TRUE,rep(FALSE, round(chr.bins / spacer.param, 0)))
      axis(side = 1, at = h2$x[spacer],
           labels = round(h2$x[spacer]/bin.size, 1),
           tick = FALSE, las=2, cex.axis=0.8)
      axis(side=2, at=c(-1.5, -0.8, -0.4, 0, 1, 2),
           labels = c("expr", "CN", "Std.Res", "AA", "AB", "BB"),
           cex.axis=0.8, las=2)
    }

    ## Visualize the Std.Res LOH track
    ab.stdres <- SchramekLOH:::.chiAB(h2$counts)  # Run the ChiSq to get the STD.RES of AB
    bin.left = h2$x.breaks[-length(h2$x.breaks)]
    bin.right = h2$x.breaks[-1]
    if(gen.plot){
      rect(xleft = bin.left, ybottom = -0.1,
           xright = bin.right, ytop = -0.5, border = FALSE,
           col = p[round(abs(ab.stdres[['thr']]), 2) * 100])
    }

    ## Visualize the CN track
    if(gen.plot){
      if(!is.null(seg.id[[chr]])){
        rect(xleft = seg.id[[chr]]$loc.start, ybottom = -1,
             xright = seg.id[[chr]]$loc.end, ytop = -0.6, border = FALSE,
             col = cn[as.character(round(seg.id[[chr]]$seg.mean,1))])
      }
      abline(h = c(0, -0.1, -0.5, -0.6, -1.0, -1.1), col="black", lty=1)
    }

    ## Visualize the Expr Track
    if(gen.plot){
      z.expr <- SchramekLOH:::.formatExprZ(z.chr[[chr]][,gsub("-", ".", ids[['pdf']])])
      rect(xleft = as.numeric(gaf.chr[[chr]]$start), ybottom = z.expr[['z0']],
           xright = as.numeric(gaf.chr[[chr]]$end), ytop = z.expr[['z']],
           col=alpha(z.expr[['col']], 0.6), border=FALSE)
      abline(h = -1.55, lty=3, col="black")
    }

    ## Combine Seg  with LOH status
    goi.tmp <- goi.chr[chr]
    if(!is.null(goi.tmp[[1]])){
      ## Visualize the Gene Of Interest is it is located in this region
      i <- goi.tmp[[1]]
      if(gen.plot){
        rect(xleft = i$end - (10000 * chr.bins), ybottom = 0,
             xright = i$end, ytop = 10,
             col=alpha("blue", 0.5), border=FALSE)
        text(x=i$end, y = 1.4, labels = i$gene, pos=2, cex=0.8)
      }

      SchramekLOH:::.intersectLohWithSeg(i, chr, bin.left, bin.right,
                                         seg.id[[chr]], ab.stdres[['raw']])
    } else {
      list("genes"=NULL,
           "stdres"=SchramekLOH:::.getStdResDf(chr, bin.left, bin.right,
                                               ab.stdres[['raw']]))
    }

  })

  if(gen.plot) dev.off()

  all.stdres <- do.call("rbind", lapply(gene.stdres, function(i) i[['stdres']]))
  stdres <- do.call("rbind", lapply(gene.stdres, function(i) i[['genes']]))
  rownames(stdres) <- NULL
  list("genes"=stdres,
       "all"=all.stdres)
}
