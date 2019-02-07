
library(Homo.sapiens)
library(taRifx) ## Removes factors
library(scales)
library(SchramekLOH)
library(gplots)
library(IdeoViz)
library(reshape)

#detach("package:SchramekLOH", unload=TRUE)
#library(SchramekLOH)

data("birdseed")  # df.bs, mapping.cov, mapping, ref.probe.ord
data("expr.mapping")  # expr.mapping
data("gaf")       # gaf
data("geneExpr")  # df.ex
data("mapping")   # mapping
data("snp6")      # snp6
data("Affyseg")   # affyseg
data("TCGAseg")   # seg
data("segmaf")    # segmaf
data("purity")    # purity


segdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc'

mut.attr <- read.table(file.path(segdir, "hnsc_tcga_extended_attributes.ripk4.txt"),
                       sep="\t", header=TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)


#### Load in Mappings ####
# Bed file obtained form http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/GenomeWideSNP_6.hg19.bed.zip
pdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2'
#outdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/summary_af'
outdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/loh_analysis'
tmpdir <- file.path(outdir, "tmp")
plotsdir <- file.path(outdir, "plots")
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotsdir, recursive = TRUE, showWarnings = FALSE)

ref <- file.path(pdir, "ref")


goi <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
         "JAG1", "JAG2", "ADAM10", "AJUBA", "RIPK4")

#### Read in Birdseed + Segs ####
use.affy <- FALSE
use.absolute <- TRUE
if(use.affy) seg <- affyseg
if(use.absolute) seg <- segmaf
 

cn.center <- c('modal_cn', 'copy_ratio')
cn.center.choice <- 1

if(use.absolute){
    # Create mapping for Absoltue SEGMAF samples
    segmaf.ids <- unique(sort(seg$sample))
    segmaf.tmp <- data.frame("sample"=segmaf.ids[which(segmaf.ids %in% mapping$V3)],
                             "ID"=gsub("01[AB]-.*", "01", 
                                       mapIds(segmaf.ids[which(segmaf.ids %in% mapping$V3)], 
                                              mapping, in.type='affy', out.type='tcga')))
    
    seg <- merge(seg, segmaf.tmp, by='sample', all.x=TRUE)
    if(cn.center[cn.center.choice] == 'modal_cn') seg$seg.mean <- (seg$modal_cn - 2)
    if(cn.center[cn.center.choice] == 'copy_ratio') seg$seg.mean <- (seg$copy_ratio - 0.5)
    
    colnames(seg) <- c('sample', 'chrom','loc.start','loc.end','num.mark','length',
                       'seg_sigma','W','copy_ratio','modal_cn','expected_cn','subclonal',
                       'cancer_cell_frac','ccf_ci95_low','ccf_ci95_high','hz','ID', 'seg.mean')
}

if(use.absolute){
  purity.ids <- purity$sample
  purity.tmp <- data.frame("sample"=purity.ids[which(purity.ids %in% mapping$V3)],
                           "ID"=gsub("01[AB]-.*", "01", 
                                     mapIds(purity.ids[which(purity.ids %in% mapping$V3)], 
                                            mapping, in.type='affy', out.type='tcga')))
  purity <- merge(purity, purity.tmp, by="sample", all.x=TRUE)
}

seg$chrom <- gsub("(chr).*\\1", "\\1", paste0("chr", seg$chrom))
seg.ids <- split(seg, f=seg$ID)
seg.i <- seg.ids[[1]]
head(seg.i)

#### Chromosome order datasets ####
seg$chrom <- gsub("(chr).*\\1", "\\1", paste0("chr", seg$chrom))
seg.ids <- split(seg, f=seg$ID)
seg.chr <- lapply(seg.ids, function(seg.i){
  seg.i <- seg.i[order(seg.i$loc.start),]
  seg.tmp <- split(seg.i, f=seg.i$chrom)
    
  chr.ord <- paste0("chr", c(1:22, "X", "Y"))
  x.factor <- factor(names(seg.tmp), levels = chr.ord, ordered=TRUE)

  seg.tmp[order(x.factor)]
})

if(use.absolute)  {
    purity.ids <- split(purity, f=purity$ID)
    print(head(names(purity.ids)))
    head(purity.ids[[1]])
}


head(seg.chr[['BALMS_p_TCGAb54and67_SNP_N_GenomeWideSNP_6_A03_730402']][['chr3']])
head(seg.chr[['TCGA-4P-AA8J-01']][['chr3']])

#### Map Probesets to Genomic Loci ####
if(exists("ref.probe.ord")){
    snp6 <- snp6[match(ref.probe.ord, snp6$V4),]
    snp6.ord <- snp6[,c(4, 1:3)]
    colnames(snp6.ord) <- c("probeset_id", "chrom", "start", "end")

    snp6.chr <- split(snp6.ord, f=snp6.ord$chrom)
    bs.chr <- split(as.data.frame(df.bs), snp6.ord$chrom)
    goi.df <- getGeneLoci(goi)
    goi.chr <- split(goi.df, f=goi.df$chr)

    chrom.ord <- match(paste0("chr", c(1:22, "X", "Y")), names(snp6.chr))
    snp6.chr <- snp6.chr[chrom.ord]
    bs.chr <- bs.chr[chrom.ord]
}


#### Expression analysis
if(exists("df.ex")){
    ## Generate z-score per gene
    z <- function(x){ (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}
    z.ex <- data.frame(t(apply(df.ex, 1, z)), stringsAsFactors=FALSE)
    
    ## Map Genes to Genomic Loci ##
    ord <- match(rownames(df.ex), gaf$V2)
    gaf.ord <- gaf[ord, c("V2", "V17")]
    gaf.ord$chr <- gsub(":.*", "", gaf.ord$V17)
    gaf.ord$start <- as.numeric(gsub("^.*:", "", gsub("-.*", "", gaf.ord$V17)))
    gaf.ord$end <- as.numeric(gsub(":.*", "", gsub("^.*?-", "", gaf.ord$V17)))
    
    ## Reorder all the matrices into genomic loci numerical order
    chr.ord <- paste0("chr", c(1:22, "X", "Y"))
    gaf.ord <- gaf.ord[order(gaf.ord$start),]
    gaf.ord <- gaf.ord[order(match(gaf.ord$chr, chr.ord)), ]

    ord <- match(gaf.ord$V2, rownames(df.ex))
    df.ex <- df.ex[ord,]
    z.ex <- z.ex[ord,]
    
    ## Order the list by chromosomes
    gaf.chr <- split(gaf.ord, f=gaf.ord$chr)
    z.chr <- split(z.ex, f=gaf.ord$chr)
    
    chrom.ord <- match(paste0("chr", c(1:22, "X", "Y")), names(gaf.chr))
    gaf.chr <- gaf.chr[chrom.ord]
    z.chr <- z.chr[chrom.ord]
}

print(paste0("Number of Affy6 samples: ", length(names(seg.ids))))
print(paste0("Number of Affy6 samples: ", length(names(z.ex))))
print(paste0("Number of Affy6 samples: ", length(names(seg.ids))))

length(intersect(names(seg.ids), gsub("\\.", "-", names(z.ex))))

library(gplots)
library(RColorBrewer)

rf <- colorRampPalette(c("white", "black"))
pf <- colorRampPalette(c("white", "white", "red", "darkred"))
cf <- colorRampPalette(c("blue4", "blue4", "blue", "white", "red", "darkred",  "darkred"))
r <- rf(32); p <- pf(1000); cn <- cf(100)
names(cn) <- round(seq(-4.9, 5.0, by=0.1),1)

#pdf(file.path(plotsdir, "legend.pdf"), width=6)
null <- split.screen(c(1,3))
screen(1); color.bar(p, min=0, max=-10, title="AB-Std.Res")
screen(2); color.bar(r, min=0, max=1, title="SNP density")
screen(3); color.bar(cn, min=-5, max=5, title="CN-LogR")
close.screen(all.screens=TRUE)
#dev.off()

bin.size <- 1000000
spacer.param <- 20

sample.stdres <- suppressWarnings(lapply(colnames(bs.chr[[1]]), mapAndPlotFeatures,
                        mapping.cov=mapping.cov, use.aff=use.affy,
                        plotsdir=plotsdir, snp6.chr=snp6.chr,
                        seg.chr=seg.chr, bs.chr=bs.chr, z.chr=z.chr,
                        r=r, p=p, cn=cn, gen.plot=TRUE, use.absolute=use.absolute))

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

head(all.stdres[['FISTS_p_TCGA_b107_121_SNP_N_GenomeWideSNP_6_D08_777884']],3)

head(sample.stdres[['FISTS_p_TCGA_b107_121_SNP_N_GenomeWideSNP_6_D08_777884']],3)

#save(all.stdres, sample.stdres, sample.stdres.bkup, file=file.path(tmpdir, paste0("tmp", use.affy, ".RData")))
#load(file.path(tmpdir, paste0("tmp", use.affy, ".RData")))

#save(all.stdres, sample.stdres, sample.stdres.bkup, file=file.path(tmpdir, paste0("tmpABSOLUTE", use.affy, ".RData")))
load(file.path(tmpdir, paste0("tmpABSOLUTE", use.affy, ".RData")))

if(use.absolute){
    seg.ids <- split(segmaf, f=segmaf$sample)
    lapply(seg.ids, plotSeg, outdir=outdir)
}

stdres.thresh <- -5
attributes <- lapply(seq_along(sample.stdres), generateIgvAttributes,
                     sample.stdres=sample.stdres, mapping.cov=mapping.cov,
                     stdres.thresh=-5)
attributes <- do.call("rbind", attributes)
head(attributes, 5)

write.table(attributes, file.path(outdir, "LOH_attributes.txt"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
save(sample.stdres, attributes, mapping.cov, file=file.path(outdir, "gene_stdRes.R"))


all.ctbl <- doTheChi(ctbl=NULL, attributes, mut.attr, 
                     tbl.idx=c(1,3), gene='ADAM10', mut='ADAM10_CNA')
ctbl <- all.ctbl[['ctbl']]

print(names(ctbl))

print(doTheChi(ctbl=ctbl, tbl.idx=c(1,3), gene='ADAM10', mut='ADAM10_CNA'))
print(doTheChi(ctbl=ctbl, tbl.idx=c(1,2), gene='NOTCH1', mut='NOTCH1_CNA'))
print(doTheChi(ctbl=ctbl, tbl.idx=c(1,3), gene='AJUBA', mut='AJUBA_CNA'))

#pdf(file.path(outdir, "stdres.pdf"), width=20)
m.stdres <- do.call("rbind", all.stdres)
m.stdres$ID <- gsub("\\..*", "", rownames(m.stdres))
m.stdres <- melt(m.stdres, id.vars = "ID", measure.vars = "stdres")
m.stdres$IDidx <- match(m.stdres$ID, unique(m.stdres$ID))
h2 <- hist2d(m.stdres[,c("IDidx", "value")],
             nbins=c(length(unique(m.stdres$ID)), 30),
             col=r, ylim=c(-20,10),
             xlab = "Genomic Loci (Mb)")
#dev.off()

ideo_hg19 <- getIdeo("hg19")

for(gene in rownames(sample.stdres[[1]])){
  pdf(file.path(outdir, paste0("ideo_", gene, ".pdf")), height=20)
  plotIdeoGene(ideo_hg19, sample.stdres, gene, thresh=-0.1)
  dev.off()
}

plotIdeoGene(ideo_hg19, sample.stdres, 'RIPK4', thresh=0)

plotIdeoGene(ideo_hg19, sample.stdres, 'ADAM10', thresh=0)

plotIdeoGene(ideo_hg19, sample.stdres, 'AJUBA', thresh=0)

c.tbl <- function(g1, g2){
    x <- lapply(names(sample.stdres), function(x.id){
        x <- sample.stdres[[x.id]]
        if(x[g1, 'seg.mean'] < 0) ad <- paste(g1, '_LOSS') else ad <- 'NA'
        if(x[g2, 'seg.mean'] < 0) aj <- paste(g2, '_LOSS') else aj <- 'NA'
        if(grepl('LOSS', ad) & grepl('LOSS', aj)) {
            id.df <- data.frame('affy'=x.id, 
                                'tcga'=gsub("01[AB]-.*", "01", 
                                            mapIds(x.id, mapping.cov, 
                                                   in.type='affy', out.type='tcga')))
        } else { id.df <- NULL }
        list("cnt"=c(ad, aj), "ids"=id.df)
    })
    
    z <- sapply(x, function(y) y[['cnt']])
    id.df <- do.call("rbind", lapply(x, function(y) y[['ids']]))
    list("tbl"=table(t(z)[,1], t(z)[,2]),
         "ids"=id.df)
}


x <- c.adam.aj[['tbl']]
x

log2((x[1,1] * x[2,2]) / (x[1,2] * x[2,1]))

#if(use.absolute){
    c.adam.aj <- c.tbl('ADAM10', 'AJUBA'); print(c.adam.aj)
    c.adam.rip <- c.tbl('ADAM10', 'RIPK4'); print(c.adam.rip[['tbl']])
    c.aj.rip <- c.tbl('AJUBA', 'RIPK4'); print(c.aj.rip[['tbl']])
#}


.reportCases <- function(id){
    id <- gsub("01[AB]-.*", "01", mapIds(id, mapping.cov, in.type='affy', out.type='tcga'))
    chr14 <- seg.chr[[id]][['chr14']]
    chr15 <- seg.chr[[id]][['chr15']]


    ajuba.idx <- which(with(chr14, loc.start < 23438410 & loc.end > 23453848))
    adam.idx <- which(with(chr15, loc.start < 58886510 & loc.end > 59044177))
    flank <- 1

    list("purity"=purity.ids[[id]][,c(1, 4,5,10 )],
         'ajuba'=chr14[c((ajuba.idx-flank):(ajuba.idx+flank)), c(17, 2:4, 6, 18)],
         'adam10'=chr15[c((adam.idx-flank):(adam.idx+flank)), c(17, 2:4, 6, 18)])
}

#lapply(c.adam.aj[['ids']][['affy']], .reportCases)

gene <- 'AJUBA'
gene.alt <- 'JUB' # Alternate name of AJUBA

overlapping.cases <- as.character(c.adam.aj[['ids']][['tcga']])
overlapping.purity <- sapply(purity.ids[overlapping.cases], function(x) x['purity'])
purity.range <- c(min(as.numeric(overlapping.purity)), 
                  max(as.numeric(overlapping.purity)))
print(purity.range)

valid.purity <- sapply(purity.ids, function(x) x['purity'] > min(purity.range) & x['purity'] < max(purity.range))
valid.samples <- names(purity.ids)[which(unlist(valid.purity))]
print(head(valid.samples))

ajuba.exp <- getGeneExp('JUB', df.ex)
ref.exp <- ajuba.exp[gsub("-", ".", valid.samples), 1, drop=FALSE]
ref.exp.num <- ref.exp[,1]
ref.z <- (ref.exp.num - mean(ref.exp.num, na.rm=TRUE)) / (sd(ref.exp.num, na.rm=TRUE))


plot(x=rep(1, length(ref.z)), ref.z)
points(x=rep(1, length(overlapping.cases)), 
       ref.z[match(gsub("-", ".", overlapping.cases), rownames(ref.exp))], 
       col='red', pch=16)

ex.by.mut <- parseIdsByMutation('AJUBA', getGeneExp('JUB', z.ex), 
                               seg.ids=seg.ids, lo.q=0.1, hi.q=0.9,
                               cbio=cbio.abs, keep.only.abs=TRUE)
null <- .mutBoxplot(ex.by.mut, add.purity=FALSE)
lapply(ex.by.mut, head, 3)

t.test(ex.by.mut[['cna']][,4], ex.by.mut[['NA']][,4])
t.test(ex.by.mut[['mut']][,4], ex.by.mut[['NA']][,4])

pdf(file.path(outdir, "ajuba_mut_expr.pdf"), width=5)
null <- .mutBoxplot(ex.by.mut, add.purity=FALSE, jitter=0.1)
dev.off()

ord <- order(ex.by.mut[['mut']][,4], decreasing=TRUE)
head(ex.by.mut[['mut']][ord,])

ex.by.mut <- parseIdsByMutation('ADAM10', getGeneExp('ADAM10', z.ex), 
                               seg.ids=seg.ids, lo.q=0.1, hi.q=0.9,
                               cbio=cbio.abs, keep.only.abs=TRUE)
null <- .mutBoxplot(ex.by.mut, add.purity=FALSE)
lapply(ex.by.mut, head, 3)

t.test(ex.by.mut[['cna']][,4], ex.by.mut[['NA']][,4])

pdf(file.path(outdir, "adam10_mut_expr.pdf"), width=5)
null <- .mutBoxplot(ex.by.mut, add.purity=FALSE, jitter=0.1)
dev.off()

genes <- unique(sort(gsub("_.*", "", colnames(mut.attr))))
genes <- genes[-grep("TRACK", genes)]

cbio.abs <- apply(mut.attr, 1, SchramekLOH::cbioConv, genes=genes, 
                  use.absolute=use.absolute, sample.stdres=sample.stdres)
cbio.abs <- do.call("rbind", cbio.abs)

cbio.meta <- apply(mut.attr, 1, SchramekLOH::cbioConv, genes=genes, 
                  use.absolute=FALSE)
cbio.meta <- do.call("rbind", cbio.meta)

head(cbio.abs[grep('AJUBA', cbio.abs$Gene),],)

head(cbio.meta[grep('AJUBA', cbio.meta$Gene),])

list(c('MIRES_p_TCGA_151_SNP_N_GenomeWideSNP_6_F06_831494', 'cnloh', 'AJUBA'),
           c('CLUBS_p_TCGA_186_188_SNP_N_GenomeWideSNP_6_E01_914086', 'purity', 'no_model'),
           c('PALPS_p_TCGA_265_266_267_N_GenomeWideSNP_6_A10_1306810', 'hetloss', 'no_model'),
           c('BUCKS_p_TCGA_272_273_N_GenomeWideSNP_6_C04_1320614', 'hetloss', 'no_model'),
           c('CONGA_p_TCGA_b_317_318_319_NSP_GenomeWideSNP_6_F06_1365262', 'NEAT', 'AJUBA'),
           c('CONGA_p_TCGA_b_317_318_319_NSP_GenomeWideSNP_6_G11_1365318', 'hetloss'),
          c("MAULS_p_TCGA_189_190_SNP_N_GenomeWideSNP_6_C11_932832 ", "NEAT", "AJUBA"))

hetloss.abs <- as.character(cbio.abs$Sample[which(with(cbio.abs, Gene == 'RIPK4' & mut == 'HETLOSS'))])
hetloss.meta <- as.character(cbio.meta$Sample[which(with(cbio.meta, Gene == 'RIPK4' & mut == 'HETLOSS'))])
if(0){
    data.frame("affy"=mapIds(setdiff(hetloss.meta, hetloss.abs), mapping.cov, in.type='tcga', out.type='affy'),
               "tcga"=setdiff(hetloss.meta, hetloss.abs))
}

## Isolate the samples containing a mutation on a particular gene
write.table(cbio.abs[which(cbio.abs$Gene == 'ADAM10'),],
           file=file.path(outdir, "adam10_absolute.tsv"),
           sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)

write.table(cbio.abs[which(cbio.abs$Gene == 'AJUBA'),],
           file=file.path(outdir, "ajuba_absolute.tsv"),
           sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)

write.table(cbio.abs[which(cbio.abs$Gene == 'RIPK4'),],
           file=file.path(outdir, "ripk4_absolute.tsv"),
           sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)

## cBioportal TSV to generate the oncoprint
write.table(cbio.abs, file=file.path(outdir, "cbioportal.tsv"),
            sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)

.genBlankDf <- function(cbio){
    missing.samples <- mut.attr$TRACK_ID[which(! mut.attr$TRACK_ID %in% cbio$Sample)]
    blank <- vector(mode = "character", length=length(missing.samples))
    blank <- rep("", length(missing.samples))
    blank.df <- data.frame("Sample"=missing.samples)
    blank.df
}

cbio.final <- cbio.abs
rm.idx <- lapply(c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "JAG1", "JAG2", "RIPK4"), function(gene){
    with(cbio.final, which(Gene == gene & mut != 'HOMDEL' & mut.type == 'CNA'))
})
rm.idx <- as.integer(as.character(unlist(rm.idx)))
cbio.final <- cbio.final[-rm.idx,]

#cbio.final <- rbind(cbio.abs, .genBlankDf(cbio.abs))
write.table(cbio.final, file=file.path(outdir, "cbioportal_final.tsv"),
            sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table( .genBlankDf(cbio.final), file=file.path(outdir, "cbioportal_final.tsv"),
            sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE, append=TRUE)

cnt <- sapply(c("ADAM10", "AJUBA", "RIPK4"), function(gene) 
    nrow(cbio.final[with(cbio.final, which(Gene == gene & mut != 'HETLOSS')),]))
data.frame("A"=cnt, "x"=(cnt / 504), "y"=(cnt / 504) * 100)

cnt <- sapply(c("ADAM10", "AJUBA"), function(gene) 
    nrow(cbio.final[with(cbio.final, which(Gene == gene & mut == 'HETLOSS')),]))
data.frame("A"=cnt, "x"=(cnt / 504), "y"=(cnt / 504) * 100)


head(sdf)

sdf <- lapply(names(all.stdres), combineSegLohPurity, seg.ids=seg.ids, all.stdres=all.stdres)

sdf <- do.call("rbind", sdf)
sdf$HETLOSS <- FALSE
sdf[which(sdf$seg.mean == -1),]$HETLOSS <- TRUE
sdf$mean <- sdf$stdres


train.sdf <- sdf[which(sdf$seg.mean == -1 | sdf$seg.mean == 0), c('HETLOSS', 'copy_ratio', 'mean',
                                                                  'purity', 'subclonal')]

glm.fit <- glm(HETLOSS ~ mean + purity, data = train.sdf, family = binomial)
summary(glm.fit)
glm.probs <- predict(glm.fit,type = "response")


boxplot(lapply(split(sdf, f=sdf$seg.mean), 
               function(x) predict(glm.fit, x, type='response')),
                   ylab="Probability", xlab="CN - 2", las=1,
                   ylim=c(0,1))

## Create a LOH/Allelic imbalance probability matrix
sdf.spl <- split(sdf[,-c(1:5)], f=sdf$ID)
loh.list <- sapply(sdf.spl, getLohProb, mod=glm.fit)
loh.mat <- round(do.call("cbind", loh.list), 2)

## Format the position of each segment
loh.pos <- sdf.spl[[1]][,c("seqnames", "start", "end")]
loh.pos <- loh.pos[-which(duplicated(loh.pos)), ]

#pdf(file.path(outdir, "loh_plots.pdf"))
split.screen(c(1,3))
screen(1); par(mar=c(5.1, 1, 2, 1)); all.loess <- plotLohProb(main="All")
screen(2); par(mar=c(5.1, 1, 2, 1)); ajuba.loess <- plotLohProb("AJUBA")
screen(3); par(mar=c(5.1, 1, 2, 1)); adam.loess <- plotLohProb("ADAM10")
close.screen(all.screens=TRUE)
#dev.off()

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
    abline(v=quantile(smoothed10, 0.01), lty=2, col="black")
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

split.screen(c(1,4))
screen(1); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(i=ajuba.loess, j=all.loess, diff=TRUE, main="Ajuba-All", xlab='Delta-Prob')
screen(2); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(i=ajuba.loess, j=adam.loess, diff=TRUE, main="Ajuba-Adam10", xlab='Delta-Prob')
screen(3); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(i=adam.loess, j=all.loess, diff=TRUE, main="Adam10-All", xlab='Delta-Prob')
screen(4); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(i=adam.loess, j=ajuba.loess, diff=TRUE, main="Adam10-Ajuba", xlab='Delta-Prob')
close.screen(all.screens=TRUE)


pdf(file.path(outdir, "gistic_plots", "delta_loess.pdf"))
split.screen(c(1,4))
screen(1); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(i=ajuba.loess, j=all.loess, diff=TRUE, main="Ajuba-All", xlab='Delta-Prob')
screen(2); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(i=ajuba.loess, j=adam.loess, diff=TRUE, main="Ajuba-Adam10", xlab='Delta-Prob')
screen(3); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(i=adam.loess, j=all.loess, diff=TRUE, main="Adam10-All", xlab='Delta-Prob')
screen(4); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(i=adam.loess, j=ajuba.loess, diff=TRUE, main="Adam10-Ajuba", xlab='Delta-Prob')
close.screen(all.screens=TRUE)
dev.off()

adam10.idx <- which(mut.attr$AJUBA_CNA != 'no_alteration')
attr.soi <- mut.attr[adam10.idx, ]$TRACK_ID
cbio.soi <- as.character(cbio.abs[which(cbio.abs$Gene == 'AJUBA'),]$Sample)

#pdf(file.path(outdir, "expr_plots.pdf"))
split.screen(c(1,3))
screen(1); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(soi=attr.soi, is.expr=TRUE, 
                                                         main='Attr AJUBA', minx=c(-1,1))
screen(2); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(soi=cbio.soi, is.expr=TRUE, 
                                                         main='cBio AJUBA', minx=c(-1,1))
screen(3); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(soi=setdiff(attr.soi, cbio.soi), is.expr=TRUE, 
                                                         main='Attr-only AJUBA', minx=c(-1,1))
close.screen(all.screens=TRUE)
#dev.off()

adam10.idx <- which(mut.attr$ADAM10_CNA != 'no_alteration')
attr.soi <- mut.attr[adam10.idx, ]$TRACK_ID
cbio.soi <- as.character(cbio.abs[which(cbio.abs$Gene == 'ADAM10'),]$Sample)

#pdf(file.path(outdir, "expr_plots.pdf"))
split.screen(c(1,3))
screen(1); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(soi=attr.soi, is.expr=TRUE, 
                                                         main='Attr ADAM10', minx=c(-1,1))
screen(2); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(soi=cbio.soi, is.expr=TRUE, 
                                                         main='cBio ADAM10', minx=c(-1,1))
screen(3); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(soi=setdiff(attr.soi, cbio.soi), is.expr=TRUE, 
                                                         main='Attr-only ADAM10', minx=c(-1,1))
close.screen(all.screens=TRUE)
#dev.off()

#pdf(file.path(outdir, "expr_plots.pdf"))
split.screen(c(1,3))
screen(1); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(goi='ADAM10', is.expr=TRUE, xlab='mean z-score')
screen(2); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(goi='AJUBA', is.expr=TRUE, xlab='mean z-score')
screen(3); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(goi='RIPK4', is.expr=TRUE, xlab='mean z-score')
close.screen(all.screens=TRUE)
#dev.off()

split.screen(c(1,4))
screen(1); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(goi='NOTCH1', is.expr=TRUE, xlab='mean z-score')
screen(2); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(goi='NOTCH2', is.expr=TRUE, xlab='mean z-score')
screen(3); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(goi='NOTCH3', is.expr=TRUE, xlab='mean z-score')
screen(4); par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(goi='NOTCH4', is.expr=TRUE, xlab='mean z-score')
close.screen(all.screens=TRUE)

dir.create(file.path(outdir, "gistic_plots"), recursive = TRUE, showWarnings = FALSE)

.outExpr <- function(gene){
    png(file.path(outdir, "gistic_plots", paste0(gene, "_expr_plots.png")), width=2.5, height=7.5, units = "in", res=300)
    par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(goi=gene, is.expr=TRUE, xlab='Mean z-score')
    dev.off()
}
.outLoh <- function(gene){
    png(file.path(outdir, "gistic_plots", paste0(gene, "_loh_plots.png")), width=2.5, height=7.5, units = "in", res=300)
    par(mar=c(5.1, 1, 2, 1)); null <- plotLohProb(goi=gene, xlab='Mean probability')
    dev.off()
}

.outExpr("ADAM10")
.outLoh("ADAM10")
.outExpr("AJUBA")
.outLoh("AJUBA")
.outExpr("RIPK4")
.outLoh("RIPK4")

train.gene <- FALSE

gene <- 'ADAM10'
mut <- paste0(gene, "_CNA")
rge <- 5
split.screen(c(1,2)); 
screen(1); range.stdres <- visOneGene(gene, rge, mut, use.absolute=TRUE, 
                                      purity=purity, plot.cr=TRUE, use.seg=TRUE)
screen(2); boxplot(lapply(split(range.stdres, f=range.stdres$cn), 
               function(x) predict(glm.fit, x, type='response')),
                   ylab="Probability", xlab="CN - 2", las=1,
                   ylim=c(0,1))
close.screen(all.screens=TRUE)

if(train.gene){
    train <- range.stdres[which(range.stdres$cn == -1 | range.stdres$cn == 0), ]
    glm.fit <- glm(HETLOSS ~ mean  + purity, data = train, family = binomial)
    summary(glm.fit)
    glm.probs <- predict(glm.fit,type = "response")
}

gene <- 'RIPK4'
mut <- paste0(gene, "_CNA")
rge <- 5
split.screen(c(1,2)); 
screen(1); range.stdres <- visOneGene(gene, rge, mut, use.absolute=TRUE, 
                                      purity=purity, plot.cr=TRUE, use.seg=TRUE)
screen(2); boxplot(lapply(split(range.stdres, f=range.stdres$cn), 
               function(x) predict(glm.fit, x, type='response')),
                   ylab="Probability", xlab="CN - 2", las=1,
                   ylim=c(0,1))
close.screen(all.screens=TRUE)

if(train.gene){
    glm.fit <- glm(HETLOSS ~ mean  + purity, data = range.stdres, family = binomial)
    summary(glm.fit)
    glm.probs <- predict(glm.fit,type = "response")
}

gene <- 'AJUBA'
mut <- paste0(gene, "_CNA")
rge <- 5

split.screen(c(1,2)); 
screen(1); range.stdres <- visOneGene(gene, rge, mut, use.absolute=TRUE, 
                                      purity=purity, plot.cr=TRUE, use.seg=TRUE)
screen(2); boxplot(lapply(split(range.stdres, f=range.stdres$cn), 
               function(x) predict(glm.fit, x, type='response')),
                   ylab="Probability", xlab="CN - 2", las=1,
                   ylim=c(0,1))
close.screen(all.screens=TRUE)

if(train.gene){
    glm.fit <- glm(HETLOSS ~ mean  + purity, data = range.stdres, family = binomial)
    summary(glm.fit)
    glm.probs <- predict(glm.fit,type = "response")
}

col.df <- data.frame("UID"=c("HOMDEL_HETLOSS", "HETLOSS_HOMDEL", "HETLOSS_HETLOSS", "HETLOSS-SNV_HETLOSS", "HETLOSS-SNV_SNV",
                             "HOMDEL_NA", "HOMDEL-SNV_NA", "HETLOSS_NA", "HETLOSS-SNV_NA", "SNV_NA", 
                             "NA_HOMDEL", "NA_HOMDEL-SNV", "NA_HETLOSS", "NA_SNV", "NA_NA"),
                     
                     "col"=c("#00441b", "#00441b", "#00441b", "#238b45", "#74c476",
                            "#a50f15", "#a50f15", "#a50f15", "#a50f15", "#fb6a4a",
                            "#08519c", "#08519c", "#08519c", "#6baed6", "grey"))

gene.y <- 'ADAM10'
gene.x <- 'AJUBA'
rge <- 5

range.stdres.xyz <- compTwoGenes(gene.x, gene.y, col.df, cex.type='xy', use.absolute=use.absolute)
range.stdres.xyz[which(range.stdres.xyz$UID == 'HETLOSS_HETLOSS'),]
plotSegSizes()

range.stdres.xyz <- compTwoGenes('ADAM10', 'RIPK4', col.df, cex.type = 'xy', 
                                 plot.legend=FALSE, use.absolute=use.absolute)
head(range.stdres.xyz[which(range.stdres.xyz$seg.y < -0.3),], 10)
#range.stdres.xyz[which(range.stdres.xyz$UID == 'HETLOSS_HETLOSS'),]

range.stdres.xyz <- compTwoGenes('AJUBA', 'RIPK4', col.df, cex.type = 'xy', 
                                 use.absolute=use.absolute, plot.legend=FALSE)
head(range.stdres.xyz[which(range.stdres.xyz$seg.y < -0.3),], 10)


