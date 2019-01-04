library(Homo.sapiens)
library(taRifx) ## Remvoes factors
library(scales)
#### Functions ####
getGeneLoci <- function(gene){
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

filterOverlaps <- function(i){
  dups <- duplicated(queryHits(i))
  if(isTRUE(any(dups))){
    i <- i[-which(dups),]
  }
  i
}

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

#### Load in Mappings ####
# Bed file obtained form http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/GenomeWideSNP_6.hg19.bed.zip
pdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2'
ref <- file.path(pdir, "ref")
snp6 <- read.table(file.path(ref, "GenomeWideSNP_6.hg19.headerless.bed"),
                   header=FALSE, sep="\t", stringsAsFactors = FALSE,
                   check.names = FALSE)

mapping <- read.table(file.path(ref, "mappingIds.txt"),
                      header=FALSE, sep="\t", stringsAsFactors = FALSE,
                      check.names = FALSE)

expr.mapping <- read.table(file.path(pdir, "gene_expression", "ref", "mappingIds.txt"),
                           header=FALSE, sep="\t", stringsAsFactors = FALSE,
                           check.names = FALSE)

gaf <- read.table(file.path(pdir, "gene_expression", "ref", "TCGA.hg19.June2011.gaf"),
                  header=FALSE, sep="\t", stringsAsFactors = FALSE,
                  check.names = FALSE)

#### Read in Birdseed + Segs ####
use.affy <- TRUE

outdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/summary_af'
vcfdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/symlink'
segdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc'
exprdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/gene_expression/symlink'
affysegdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/cnv_seg/symlink'

goi <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
         "JAG1", "JAG2", "ADAM10", "AJUBA")

seg <- read.table(file.path(segdir, "hnsc_tcga_segments.seg"),
                  header=TRUE, stringsAsFactors = FALSE,
                  check.names=FALSE)
affyseg <- read.table(file.path(affysegdir, "TCGA_HNSC_affy6.seg"),
                      header=TRUE, stringsAsFactors = FALSE,
                      check.names = FALSE)
if(use.affy) seg <- affyseg

#### Load in all Birdseed files and aggregate to a matrix
birdseed <- list.files(vcfdir, pattern = "birdseed.data.txt")
all.bs <- lapply(mapping$V3, function(i){
  f <- file.path(vcfdir, birdseed[grep(i, birdseed)])
  if(file.exists(f)){
    bs <- read.table(f, header=FALSE, skip = 2, 
                     stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    warning(paste0("Could not find: ", i))
    bs <- data.frame(matrix(ncol=3))
  }
  bs[,1:2]
})

## Aggregate to a matrix
match.idx <- sapply(all.bs, nrow) > 1
ref.probe.ord <- all.bs[match.idx][[1]]$V1
df.bs <- sapply(all.bs[match.idx], function(e.bs){
  if(!all(e.bs$V1 == ref.probe.ord)){
    warning("Non-matching probe order")
    ord <- match(ref.probe.ord, e.bs$V1)
    e.bs <- e.bs[ord,]
  }
  e.bs$V2
})

mapping.cov <- mapping[which(match.idx),]
mapping.noncov <- mapping[which(!match.idx),]
colnames(df.bs) <- mapping.cov$V3
rownames(df.bs) <- ref.probe.ord


#### Load in all expression calls and aggreagte to a matrix
rpkmfiles <- list.files(exprdir, pattern = "trimmed.annotated.gene.quantification.txt")
all.ex <- lapply(mapping$V2, function(i){
  i <- gsub("01[AB]-.*", "01", i)
  f <- file.path(exprdir, rpkmfiles[grep(i, rpkmfiles)])
  if(isTRUE(file.exists(f))){
    ex <- read.table(f, header=TRUE,  
                     stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    warning(paste0("Could not find: ", i))
    ex <- data.frame(matrix(ncol=4))
    colnames(ex) <- c("gene", "raw_counts", 
                      "median_length_normalized", "RPKM")
  }
  ex[,c(1,4)]
})

## Aggregate to a matrix
ref.gene.ord <- all.ex[[1]]$gene
df.ex <- sapply(all.ex, function(e.ex){
  if(!isTRUE(all(e.ex$gene == ref.gene.ord))){
    warning("Non-matching probe order")
    ord <- match(ref.gene.ord, e.ex$gene)
    e.ex <- e.ex[ord,]
  }
  e.ex$RPKM
})
colnames(df.ex) <- gsub("01[AB]-.*", "01", mapping$V2)
rownames(df.ex) <- all.ex[[1]]$gene

## Generate z-score per gene
z <- function(x){ (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}
z.ex <- data.frame(t(apply(df.ex, 1, z)), stringsAsFactors=FALSE)


#### Map Probesets to Genomic Loci ####
snp6 <- snp6[match(ref.probe.ord, snp6$V4),]
snp6.ord <- snp6[,c(4, 1:3)]
colnames(snp6.ord) <- c("probeset_id", "chrom", "start", "end")

#### Map Genes to Genomic Loci ####
ord <- match(rownames(df.ex), gaf$V2)
gaf.ord <- gaf[ord, c("V2", "V17")]
gaf.ord$chr <- gsub(":.*", "", gaf.ord$V17)
gaf.ord$start <- as.numeric(gsub("^.*:", "", gsub("-.*", "", gaf.ord$V17)))
gaf.ord$end <- as.numeric(gsub(":.*", "", gsub("^.*?-", "", gaf.ord$V17)))

## Reorder all the matrices
chr.ord <- paste0("chr", c(1:22, "X", "Y"))
gaf.ord <- gaf.ord[order(gaf.ord$start),]
gaf.ord <- gaf.ord[order(match(gaf.ord$chr, chr.ord)), ]

ord <- match(gaf.ord$V2, rownames(df.ex))
df.ex <- df.ex[ord,]
z.ex <- z.ex[ord,]


#### Chromosome order datasets ####
seg$chrom <- gsub("(chr).*\\1", "\\1", paste0("chr", seg$chrom))
seg.ids <- split(seg, f=seg$ID)
seg.chr <- lapply(seg.ids, function(seg.i){
  seg.tmp <- split(seg.i, f=seg.i$chrom)
  chrom.ord <- match(paste0("chr", c(1:22)), names(seg.tmp))
  seg.tmp[chrom.ord]
})

snp6.chr <- split(snp6.ord, f=snp6.ord$chrom)
bs.chr <- split(as.data.frame(df.bs), snp6.ord$chrom)
goi.df <- getGeneLoci(goi)
goi.chr <- split(goi.df, f=goi.df$chr)
gaf.chr <- split(gaf.ord, f=gaf.ord$chr)
z.chr <- split(z.ex, f=gaf.ord$chr)

chrom.ord <- match(paste0("chr", c(1:22, "X", "Y")), names(snp6.chr))
snp6.chr <- snp6.chr[chrom.ord]
bs.chr <- bs.chr[chrom.ord]
chrom.ord <- match(paste0("chr", c(1:22, "X", "Y")), names(gaf.chr))
gaf.chr <- gaf.chr[chrom.ord]
z.chr <- z.chr[chrom.ord]

#### Visualization and LOH test ####
library(gplots)
library(RColorBrewer)
rf <- colorRampPalette(c("white", "black"))
pf <- colorRampPalette(c("white", "white", "red", "darkred"))
cf <- colorRampPalette(c("blue4", "blue4", "blue", "white", "red", "darkred",  "darkred"))
r <- rf(32); p <- pf(1000); cn <- cf(100)
names(cn) <- round(seq(-4.9, 5.0, by=0.1),1)

## Plot Legend
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
pdf(file.path(outdir, "legend.pdf"), width=6)
split.screen(c(1,3))
screen(1); color.bar(p, min=0, max=-10, title="AB-Std.Res")
screen(2); color.bar(r, min=0, max=1, title="SNP density")
screen(3); color.bar(cn, min=-5, max=5, title="CN-LogR")
close.screen(all.screens=TRUE)
dev.off()

bin.size <- 1000000
spacer.param <- 20
sample.stdres <- lapply(colnames(bs.chr[[1]]), function(id){
  tcga.id <- mapping.cov[match(id, mapping.cov$V3), 'V2']
  tcga.id <- gsub("01[AB]-.*", "01", tcga.id)
  pdf.id <- tcga.id
  if(use.affy){
    tcga.id <- id
  }
  
  pdf(file.path(outdir, paste0(pdf.id, ".", use.affy, ".pdf")), height=50)
  split.screen(c(length(names(snp6.chr)), 1))
  seg.id <- seg.chr[[tcga.id]]
  
  gene.stdres <- lapply(names(snp6.chr), function(chr){
    screen(match(chr, names(snp6.chr)))
    par(mar=c(3.1, 4.5, 0.5, 2.1))
    
    df <- data.frame(x = snp6.chr[[chr]]$start, 
                     y = bs.chr[[chr]][,id])
    chr.bins <- round(max(snp6.chr[[chr]]$start) / bin.size, 0)
    
    h2 <- hist2d(df, nbins=c(chr.bins, 5), col=r, ylim=c(-2.2,2), 
                 xlim=c(0, max(snp6.chr[[chr]]$start)), 
                 yaxt='n', xaxt='n',
                 xlab = "Genomic Loci (Mb)", ylab="")
    title(ylab=chr, line = 3.5)
    
    spacer <- c(TRUE,rep(FALSE, round(chr.bins / spacer.param, 0)))
    axis(side = 1, at = h2$x[spacer], 
         labels = round(h2$x[spacer]/bin.size, 1), 
         tick = FALSE, las=2, cex.axis=0.8)
    axis(side=2, at=c(-1.5, -0.8, -0.4, 0, 1, 2), 
         labels = c("expr", "CN", "Std.Res", "AA", "AB", "BB"), 
         cex.axis=0.8, las=2)
    
    qval <- apply(h2$counts, 1, function(i){ 
      p <- tryCatch({ chisq.test(i[c(1,3,5)])$stdres[2]}, 
                    error=function(e){NA})
      p
    })
    qval.raw <- qval
    qval[qval <= -10] <- -10
    
    bin.left = h2$x.breaks[-length(h2$x.breaks)]
    bin.right = h2$x.breaks[-1]
    rect(xleft = bin.left, ybottom = -0.1, 
         xright = bin.right, ytop = -0.5, border = FALSE,
             col = p[round(abs(qval), 2) * 100])
    if(!is.null(seg.id[[chr]])){
      rect(xleft = seg.id[[chr]]$loc.start, ybottom = -1, 
           xright = seg.id[[chr]]$loc.end, ytop = -0.6, border = FALSE,
           col = cn[as.character(round(seg.id[[chr]]$seg.mean,1))])
    }
    abline(h = c(0, -0.1, -0.5, -0.6, -1.0, -1.1), col="black", lty=1)
    
    ## Add in Expr Track
    z.val <- z.chr[[chr]][,gsub("-", ".", pdf.id)]
    z.to.y <- z.val
    z.to.y[z.to.y > 2] <- 2
    z.to.y[z.to.y < -2] <- -2
    range <- c(-2.1, -1.1)
    z.to.y <- rescale(z.to.y, to=range, from=c(-2, 2))
    mid.y <- (c(range[1] - range[2]) /2) + (range[2])
    z.col <- rep("red", length(z.to.y))
    z.col[z.to.y < mid.y] <- 'blue'
    
    rect(xleft = as.numeric(gaf.chr[[chr]]$start), ybottom = mid.y, 
         xright = as.numeric(gaf.chr[[chr]]$end), ytop = (z.to.y), 
         col=alpha(z.col, 0.6), border=FALSE)
    abline(h = -1.55, lty=3, col="black")
    
    ## Combine Seg  with LOH status
    gr.bin <- GRanges(seqnames = chr, ranges=IRanges(bin.left, bin.right))
    qval.df <- data.frame("chrom"=as.character(seqnames(gr.bin)),
                          "start"=start(gr.bin),
                          "end"=end(gr.bin),
                          "stdres"=qval.raw)
    
    goi.tmp <- goi.chr[chr]
    if(!is.null(goi.tmp[[1]])){
      i <- goi.tmp[[1]]
      rect(xleft = i$end - (10000 * chr.bins), ybottom = 0, 
           xright = i$end, ytop = 10, 
           col=alpha("blue", 0.5), border=FALSE)
      text(x=i$end, y = 1.4, labels = i$gene, pos=2, cex=0.8)  
      
      gr.i <- GRanges(i)
      gr.bin <- GRanges(seqnames = chr, ranges=IRanges(bin.left, bin.right))
      gr.seg.reduce <- cleanSeg(seg.id[[chr]])
      gr.seg <- GRanges(seg.id[[chr]])
      
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
                              "bin.start"=start(gr.bin[ov.idx]),
                              "bin.end"=end(gr.bin[ov.idx]),
                              "seg.start"=start(gr.seg.reduce[red.ov.idx]),
                              "seg.end"=end(gr.seg.reduce[red.ov.idx]),
                              "seg.mean"=gr.seg[seg.ov.idx]$seg.mean,
                              "stdres"=qval[ov.idx]),
           "stdres"=qval.df)
      
    } else {
      list("genes"=NULL,
           "stdres"=qval.df)
    }
    
  })
  dev.off()
  
  all.stdres <- do.call("rbind", lapply(gene.stdres, function(i) i[['stdres']]))
  stdres <- do.call("rbind", lapply(gene.stdres, function(i) i[['genes']]))
  rownames(stdres) <- NULL
  list("genes"=stdres,
       "all"=all.stdres)
})
all.stdres <- lapply(sample.stdres, function(i) i[['all']])
names(all.stdres) <- colnames(bs.chr[[1]])
## Reduce the gene to a single segment
sample.stdres.bkup <- sample.stdres
sample.stdres <- lapply(sample.stdres, function(i) {
  single.j <- sapply(split(i[['genes']], f=i[['genes']]$gene), function(j){
    uniq.j <- apply(j, 2, unique)
    if(any(sapply(uniq.j[c('seg.start', 'seg.end', 'seg.mean')], length) > 1)){
      uniq.j[['seg.start']] <- min(uniq.j[['seg.start']])
      uniq.j[['seg.end']] <- max(uniq.j[['seg.end']])
      uniq.j[['seg.mean']] <- mean(uniq.j[['seg.mean']])
    }
    
    sapply(uniq.j, function(x) x)
  })
  remove.factors(data.frame(t(single.j)))
})
names(sample.stdres) <- colnames(bs.chr[[1]])

#### Format Attributes file using SEG IDs####
stdres.thresh <- -5
attributes <- lapply(seq_along(sample.stdres), function(idx){
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
})
attributes <- do.call("rbind", attributes)

write.table(attributes, file.path(outdir, "LOH_attributes.txt"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
save(sample.stdres, attributes, mapping.cov, file=file.path(outdir, "gene_stdRes.R"))


#### EXTRA: Comparison of Attributes ####
loh.attr <- remove.factors(attributes)
mut.attr <- read.table(file.path(segdir, "hnsc_tcga_attributes.txt"),
                       sep="\t", header=TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)
m.attr <- remove.factors(merge(mut.attr, loh.attr, by="TRACK_ID", all=TRUE))
m.attr[is.na(m.attr)] <- 'Unknown'

acceptable.terms <- c("HOMDEL", "HETLOSS", "LOH", "Het", "no_alteration", 'Unknown')
ids <- m.attr$TRACK_ID
m.attr <- as.matrix(m.attr)
m.attr[!m.attr %in% acceptable.terms] <- 'Mut'
m.attr[,1] <- ids

ctbl <- apply(m.attr[,-1], 2, function(i){
  apply(m.attr[,-1], 2, function(j){
    table(data.frame(i,j))
  })
})

doTheChi <- function(tbl, idx=c(1,3)){
  ch <- chisq.test(tbl[,idx])
  list(tbl, ch$p.value, ch$stdres)
}

doTheChi(ctbl[['ADAM10']][['ADAM10_CNA']])
doTheChi(ctbl[['NOTCH1']][['NOTCH1_CNA']], idx=c(1,2))
doTheChi(ctbl[['AJUBA']][['AJUBA_CNA']])


#### EXTRA: Distribution of STDRES
pdf(file.path(outdir, "stdres.pdf"), width=20)
library(reshape)
m.stdres <- do.call("rbind", all.stdres)
m.stdres$ID <- gsub("\\..*", "", rownames(m.stdres))
m.stdres <- melt(m.stdres, id.vars = "ID", measure.vars = "stdres")
m.stdres$IDidx <- match(m.stdres$ID, unique(m.stdres$ID))
h2 <- hist2d(m.stdres[,c("IDidx", "value")], 
             nbins=c(length(unique(m.stdres$ID)), 30), 
             col=r, ylim=c(-20,10), 
             xlab = "Genomic Loci (Mb)", ylab=chr)
dev.off()

#### EXTRA: Ideogram a gene
require(IdeoViz)
ideo_hg19 <- getIdeo("hg19")
getGeneGr <- function(x, gene){
  GRanges(x[gene,  c('chr', 'seg.start', 'seg.end', 'seg.mean')])
}

plotIdeoGene <- function(ideo.ref, sample.stdres, 
                         gene, thresh=-0.1){
  gr.list <- lapply(sample.stdres, getGeneGr, gene=gene)
  names(gr.list) <- gsub("01[AB]-.*", "01", mapping.cov$V2)
  gr.seg.mean <- sapply(gr.list, function(x) x$seg.mean)
  gr.list <- gr.list[which(as.numeric(gr.seg.mean) < thresh)]
  gr.list.ord <- order(sapply(gr.list, width))
  #load("~/Desktop/tmp.RData")
  
  #save(ideo_hg19, gr.list, sample.stdres, mapping.cov, file="~/Desktop/tmp2.RData")
  #load("~/Desktop/tmp2.RData")
  
  pdf(file.path(outdir, paste0("ideo_", gene, ".pdf")), height=20)
  plotOnIdeo(chrom=as.character(seqnames(gr.list[[1]])),
             ideoTable=ideo_hg19,
             values_GR=gr.list[gr.list.ord], value_cols="seg.mean",
             plotType="seg_tracks",
             col=alpha("blue",0.6), vertical=F)
  dev.off()
}

for(gene in rownames(sample.stdres[[1]])){
  plotIdeoGene(ideo_hg19, sample.stdres, gene, thresh=-0.1)
}
x <- sapply(sample.stdres, function(x) x['ADAM10', c('seg.mean', 'stdres')])
x <- as.data.frame(t(x))
storage.mode(x) <- "numeric"
x <- x[which(x$seg.mean < -0.2),]
plot(x)
cor(as.numeric(x[,1]), as.numeric(x[,2]))
