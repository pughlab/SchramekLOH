#' plotSeg
#'
#' @param seg
#' @param chr
#' @param cn
#' @param ylim
#' @param outdir
#'
#' @return
#' @export
#'
#' @examples
plotSeg <- function(seg, chr='Chromosome', cn=NULL,
                    ylim=c(0,1), outdir){
  col.li <- list('0'='darkblue',
                 '1'='blue',
                 '2'='grey',
                 '3'='red',
                 '4'='darkred')
  seg$modal_cn[which(seg$modal_cn > 4)] <- 4
  seg.chr <- split(seg, f=seg[,chr])

  if(all(grepl("chr", names(seg.chr)))){
    chr.order <- paste0("chr", c(1:22, "X", "Y"))
  } else {
    chr.order <- c(1:22, "X", "Y")
  }
  seg.chr <- seg.chr[chr.order]
  seg.chr <- seg.chr[!is.na(names(seg.chr))]

  dir.create(file.path(outdir, "abs_plots"))
  pdf(file.path(outdir, "abs_plots", paste0(unique(seg$sample), ".pdf")),
      width=8, height=30)
  split.screen(c(length(seg.chr), 1))

  lapply(names(seg.chr), function(s.id){
    s <- seg.chr[[s.id]]

    screen(grep(paste0("^", s.id, "$"), names(seg.chr)));
    par(mar=c(3.1, 4.1, 0.1, 2.1))

    plot(0, type='n', ylim=ylim, xlim=c(0, max(s$End.bp)),
         xlab='Genomic Coordinate', ylab=paste0("Chr", s.id))
    rect(xleft = s$Start.bp, ybottom = (s$copy_ratio - 0.05),
         xright = s$End.bp, ytop = (s$copy_ratio + 0.05),
         col=as.character(col.li[as.character(s$modal_cn)]),
         border=FALSE)
    abline(h = median(seg$copy_ratio), lty=2, col="black")
  })
  close.screen(all.screens=TRUE)
  dev.off()

}
