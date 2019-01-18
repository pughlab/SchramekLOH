

#' filterOverlaps
#'
#' @param i The returned object for GRanges overlap
#'
#' @return A filtered version of GRanges overlap with a
#' unique queryHits
#' @export
#'
#' @examples
filterOverlaps <- function(i){
  dups <- duplicated(queryHits(i))
  if(isTRUE(any(dups))){
    i <- i[-which(dups),]
  }
  i
}

#' cleanSeg
#'
#' @param seg.x A seg file in data frame, will be converted ot GRanges
#' @param min.size The minimum size of a segment to be allowed
#' @param loss.thresh Anything lower than this will be classified as a loss
#' @param gain.thresh Anything higher than this will be classified as a gain
#'
#' @return
#' @export
#'
#' @examples
cleanSeg <- function(seg.x, min.size=50000,
                     loss.thresh=-0.1, gain.thresh=0.1){
  x.gr <- GRanges(seg.x)
  x.gr <- x.gr[width(x.gr) > min.size,] # Remove small CNVs

  x.gr$gain <- x.gr$seg.mean > gain.thresh
  x.gr$loss <- x.gr$seg.mean < loss.thresh
  x.gr$neut <- !(x.gr$gain | x.gr$loss)

  x.gr.g <- reduce(x.gr[x.gr$gain,], min.gapwidth=min.size)
  x.gr.l <- reduce(x.gr[x.gr$loss,], min.gapwidth=min.size)
  x.gr.n <- reduce(x.gr[x.gr$neut,], min.gapwidth=min.size)

  x.gr <- c(x.gr.g, x.gr.l, x.gr.n)
  x.gr
}

#' mapIds
#'
#' @param id input ID (e.g. TCGA-CN-6011-01)
#' @param mapping Mapping data frame
#' @param in.type either 'tcga' (TCGA-CN-6011-01) or
#'   affy ('FISTS_p_TCGA_b107_121_SNP_N_GenomeWideSNP_6_D08_777884')
#' @param out.type either 'tcga' or 'affy'
#'
#' @return
#' @export
#'
#' @examples mapIds(id, mapping.cov, in.type='tcga', out.type='affy')
mapIds <- function(id, mapping, in.type='tcga',
                   out.type='affy', clean=TRUE){
  in.idx <- switch(in.type,
                   tcga=2,
                   affy=3)
  out.idx <- switch(out.type,
                    tcga=2,
                    affy=3)

  idx <- sapply(id, function(i) grep(i, mapping[,in.idx])[1])
  if(out.type == 'tcga' & clean){
    gsub("01[AB]-.*", "01", mapping[idx, out.idx])
  } else {
    mapping[idx, out.idx]
  }
}


#' aggregateStdRes
#'
#' @param sample.id
#' @param range
#' @param gene
#' @param sample.stdres
#' @param all.stdres
#'
#' @return
#' @export
#'
#' @examples
aggregateStdRes <- function(sample.id, range, gene, sample.stdres,
                            all.stdres, use.absolute=FALSE, use.seg=TRUE){
  gene.df <- sample.stdres[[sample.id]]
  seg.df <- all.stdres[[sample.id]]

  gene.info <- gene.df[match(gene, gene.df$gene),]
  if(use.seg){
    seg.range <- intersect(which(seg.df$chrom %in% gene.info$chr),
                           which(seg.df$start > as.numeric(gene.info$seg.start) &
                                   seg.df$end < as.numeric(gene.info$seg.end)))
  } else {
    bin.idx <- intersect(which(seg.df$chrom %in% gene.info$chr),
                         which(seg.df$start %in% gene.info$bin.start))
    seg.range <- c((bin.idx-range):(bin.idx+range))
  }
  rge.stdres <- seg.df[seg.range,]$stdres

  c("mean"=mean(rge.stdres, na.rm=TRUE),
    "sd"=sd(rge.stdres, na.rm=TRUE),
    "seg"=as.numeric(gene.info[,'seg.mean']),
    "copy_ratio"=if(use.absolute) as.numeric(gene.info[,'copy.ratio']) else 0,
    "seg.start"=as.numeric(gene.info[,'seg.start']),
    "seg.end"=as.numeric(gene.info[,'seg.end']))
}

#' getGeneExp
#'
#' @param gene.id i.e. AJUBA or JUB
#'
#' @return
#' @export
#'
#' @examples
getGeneExp <- function(gene.id, z.ex=z.ex){
  gene.ex <- data.frame(t(z.ex[grep(paste0("^", gene.id), rownames(z.ex)), ,drop=FALSE]),
                        stringsAsFactors=FALSE)
  gene.ex$TRACK_ID <- gsub("\\.", "-", rownames(gene.ex))
  gene.ex
}


#' getSegIQR
#'
#' @param segdf
#' @param lo.q
#' @param hi.q
#'
#' @return
#' @export
#'
#' @examples
getSegIQR <- function(segdf, lo.q=0.25, hi.q=0.75){
  lrr <- rep(segdf$seg.mean, round(segdf$num.mark/10,0))
  iqr <- quantile(lrr, hi.q) - quantile(lrr, lo.q)
  c("TRACK_ID" = unique(as.character(segdf$ID)), "IQR"=iqr)
}

#' parseIdsByMutation
#'
#' @param gene
#' @param gene.ex
#'
#' @return
#' @export
#'
#' @examples
parseIdsByMutation <- function(gene, gene.ex, seg.ids=seg.ids, lo.q=0.1, hi.q=0.9){
  .getAltType <- function(x){
    alt.x <- rep("NA", length(x))
    non.idx <- unique(sort(unlist(sapply(c("fs", "del", "\\*", "splice"), grep, x=x))))
    mis.idx <- grep("^[A-Z][0-9]+[A-Z]$", x=x)

    alt.x[non.idx] <- 'Nonsense'
    alt.x[mis.idx] <- 'Missense'
    alt.x
  }

  ## Get a list of all the mutation subtype (e.g. CNA, MUT, FUSION)
  gene.mut <- colnames(mut.attr)[grep(gene, colnames(mut.attr))]

  ## Get the IQR for seg LRR as a pseudo-approximation for purity
  seg.iqr <- data.frame(t(sapply(seg.ids, getSegIQR, lo.q=lo.q, hi.q=hi.q)),
                        stringsAsFactors=FALSE)

  ## Find Samples that are "no_alteration" across all mutation subtypes
  no_alt.samples <- lapply(gene.mut, function(m) mut.attr[which(mut.attr[,m] == 'no_alteration'), 'TRACK_ID'])
  no_alt.samples <- Reduce(function(x,y) intersect(x,y), no_alt.samples)

  ## Identify samples for each particular mutation
  # Samples with no_alteration (CNA, MUT, or FUSION)
  no_alt.samples <- lapply(gene.mut, function(m) mut.attr[which(mut.attr[,m] == 'no_alteration'), 'TRACK_ID'])
  no_alt.samples <- Reduce(function(x,y) intersect(x,y), no_alt.samples)
  no_alt.samples <- data.frame("TRACK_ID"=no_alt.samples, "Alt"="no_alteration")
  no_alt.samples$Alt.type <- .getAltType(no_alt.samples[,2])
  no_alt.samples <- merge(no_alt.samples, gene.ex, by='TRACK_ID', all.x=TRUE)
  no_alt.samples <- merge(no_alt.samples, seg.iqr, by='TRACK_ID', all.x=TRUE)

  # Samples with specific types of alterations
  mut.samples <- mut.attr[which(mut.attr[, paste0(gene, "_MUT")] != 'no_alteration'),
                          c('TRACK_ID',paste0(gene, "_MUT"))]
  mut.samples$Alt.type <- .getAltType(mut.samples[,2])
  mut.samples <- merge(mut.samples, gene.ex, by='TRACK_ID', all.x=TRUE)
  mut.samples <- merge(mut.samples, seg.iqr, by='TRACK_ID', all.x=TRUE)

  cna.samples <- mut.attr[which(mut.attr[, paste0(gene, "_CNA")] != 'no_alteration'),
                          c('TRACK_ID', paste0(gene, "_CNA"))]
  cna.samples$Alt.type <- .getAltType(cna.samples[,2])
  cna.samples <- merge(cna.samples, gene.ex, by='TRACK_ID', all.x=TRUE)
  cna.samples <- merge(cna.samples, seg.iqr, by='TRACK_ID', all.x=TRUE)

  fusion.samples <- mut.attr[which(mut.attr[, paste0(gene, "_FUSION")] != 'no_alteration'),
                             c('TRACK_ID', paste0(gene, "_FUSION"))]
  fusion.samples$Alt.type <- .getAltType(fusion.samples[,2])
  fusion.samples <- merge(fusion.samples, gene.ex, by='TRACK_ID', all.x=TRUE)
  fusion.samples <- merge(fusion.samples, seg.iqr, by='TRACK_ID', all.x=TRUE)

  ex.by.mut <- list("NA"=no_alt.samples,
                    "mut"=mut.samples,
                    "cna"=cna.samples,
                    "fusion"=fusion.samples)
  ex.by.mut
}


#' .mutBoxplot
#'
#' @param ex.by.mut
#' @param add.purity
#'
#' @return
#' @export
#'
#' @examples
.mutBoxplot <- function(ex.by.mut, add.purity=TRUE){
  cols <- list("NA"='grey',
               'Nonsense'='red',
               'Missense'='blue')

  ex.id <- colnames(ex.by.mut[[1]])[4]

  boxplot(lapply(ex.by.mut, function(x) x[,ex.id ]),
          ylab="Gene z-score", ylim=c(-2, 3), outline=FALSE)

  x <- sapply(seq_along(ex.by.mut), function(x){
    ex <- ex.by.mut[[x]]
    points(x=rep(x, nrow(ex)), y=ex[, ex.id],
           pch=16,
           col=alpha(cols[as.character(ex$Alt.type)], 0.6),
           cex=if(add.purity) rescale(as.numeric(ex[,5]), to=c(1,2)) else 1)
  })
  NULL
}
