
library(Homo.sapiens)
library(taRifx) ## Removes factors
library(scales)
library(SchramekLOH)
library(gplots)
library(IdeoViz)
library(reshape)

detach("package:SchramekLOH", unload=TRUE)
require(SchramekLOH)

data("birdseed")  # df.bs, mapping.cov, mapping, ref.probe.ord
data("expr.mapping")  # expr.mapping
data("gaf")       # gaf
data("geneExpr")  # df.ex
data("mapping")   # mapping
data("snp6")      # snp6
data("Affyseg")   # affyseg
data("TCGAseg")   # seg


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
         "JAG1", "JAG2", "ADAM10", "AJUBA")

#### Read in Birdseed + Segs ####
use.affy <- FALSE
if(use.affy) seg <- affyseg

head(seg)

#### Chromosome order datasets ####
seg$chrom <- gsub("(chr).*\\1", "\\1", paste0("chr", seg$chrom))
seg.ids <- split(seg, f=seg$ID)
seg.chr <- lapply(seg.ids, function(seg.i){
  seg.tmp <- split(seg.i, f=seg.i$chrom)
  chrom.ord <- match(paste0("chr", c(1:22)), names(seg.tmp))
  seg.tmp[chrom.ord]
})

head(seg.chr[['BALMS_p_TCGAb54and67_SNP_N_GenomeWideSNP_6_A03_730402']][['chr3']])

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
                        r=r, p=p, cn=cn, gen.plot=TRUE))

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
load(file.path(tmpdir, paste0("tmp", use.affy, ".RData")))

stdres.thresh <- -5
attributes <- lapply(seq_along(sample.stdres), generateIgvAttributes,
                     sample.stdres=sample.stdres, mapping.cov=mapping.cov,
                     stdres.thresh=-5)
attributes <- do.call("rbind", attributes)
head(attributes, 5)

write.table(attributes, file.path(outdir, "LOH_attributes.txt"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
save(sample.stdres, attributes, mapping.cov, file=file.path(outdir, "gene_stdRes.R"))


segdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc'

mut.attr <- read.table(file.path(segdir, "hnsc_tcga_attributes.txt"),
                       sep="\t", header=TRUE, stringsAsFactors = FALSE,
                       check.names = FALSE)

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

plotIdeoGene(ideo_hg19, sample.stdres, 'ADAM10', thresh=-0.1)

plotIdeoGene(ideo_hg19, sample.stdres, 'AJUBA', thresh=-0.1)

gene <- 'AJUBA'
gene.alt <- 'JUB' # Alternate name of AJUBA

.mutBoxplot <- function(ex.by.mut){
    ex.id <- colnames(ex.by.mut[[1]])[3]

    boxplot(lapply(ex.by.mut, function(x) x[,ex.id ]), outline=FALSE)

    x <- sapply(seq_along(ex.by.mut), function(x){
        ex <- ex.by.mut[[x]]
        points(x=rep(x, nrow(ex)), y=ex[, ex.id],
              pch=16, col=alpha("black", rescale(as.numeric(ex[,4]), to=c(0.2,0.8))),
              cex=rescale(as.numeric(ex[,4]), to=c(1,2)))
    })
    NULL
}

ex.by.mut <- parseIdsByMutation('AJUBA', getGeneExp('JUB'), 
                               seg.ids=seg.ids, lo.q=0.1, hi.q=0.9)
null <- .mutBoxplot(ex.by.mut)
lapply(ex.by.mut, head, 3)

ord <- order(ex.by.mut[['mut']][,3], decreasing=TRUE)
ex.by.mut[['mut']][ord,]

#gene <- 'ADAM10'
gene <- 'AJUBA'
mut <- paste0(gene, "_CNA")
rge <- 5

range.stdres <- sapply(names(sample.stdres), aggregateStdRes,
                       range=rge, gene=gene, sample.stdres=sample.stdres, all.stdres=all.stdres)
range.stdres <- data.frame(t(range.stdres), stringsAsFactors=FALSE)

hetloss.idx <- which(mut.attr[,mut]  == 'HETLOSS')
homloss.idx <- which(mut.attr[,mut]  == 'HOMDEL')
mut.ids <- mut.attr[c(hetloss.idx, homloss.idx), 'TRACK_ID', drop=TRUE]
mapped.ids <- mapIds(mut.ids, mapping.cov, in.type='tcga', out.type='affy')

range.stdres$HETLOSS <- FALSE
range.stdres[which(rownames(range.stdres) %in% mapped.ids),]$HETLOSS <- TRUE

head(range.stdres)

with(range.stdres, plot(mean~seg, cex=sd, xlim=c(-1, 0.5),
                        col=alpha('grey', 1-rescale(sd, to=c(0,1))), 
                        pch=16, ylab="StdRes", xlab="Seg"))
axis(side=1, at = seq(-5, 4, by=0.1), 
     labels = rep("", length(seq(-5, 4, by=0.1))), tick = TRUE)

with(range.stdres[which(range.stdres$HETLOSS),], 
     points(mean~seg, cex=sd, 
          col=alpha('red', 1-rescale(sd, to=c(0,1))), 
          pch=16))

# Sample with a Homozygous Deletion
print("Homozygous Deletion case")
print(paste(c(gsub("01[AB]-.*", "01", mut.attr[homloss.idx,'TRACK_ID']),
              mapIds(mut.attr[homloss.idx,'TRACK_ID'], mapping.cov)), collapse=", "))
range.stdres[mapIds(mut.attr[homloss.idx,'TRACK_ID'], mapping.cov),]



# Sample with a -3.0 segment
print("Low LRR case")
id.x <- mapIds(rownames(range.stdres[which(range.stdres$seg < -1),]), mapping.cov,
               in.type='affy', out.type='tcga')
print(paste(c(gsub("01[AB]-.*", "01", id.x),
              rownames(range.stdres[which(range.stdres$seg < -1),])), collapse=", "))

mut.attr[mut.attr$TRACK_ID == gsub("01[AB]-.*", "01", id.x), , drop=FALSE]
range.stdres[mapIds(id.x , mapping.cov, in.type="tcga", out.type="affy"),]

gene.y <- 'ADAM10'
gene.x <- 'AJUBA'
rge <- 5

col.df <- data.frame("UID"=c("HOMDEL_HETLOSS", "HETLOSS_HOMDEL", "HETLOSS_HETLOSS", "HETLOSS-SNV_HETLOSS", "HETLOSS-SNV_SNV",
                             "HOMDEL_NA", "HOMDEL-SNV_NA", "HETLOSS_NA", "HETLOSS-SNV_NA", "SNV_NA", 
                             "NA_HOMDEL", "NA_HOMDEL-SNV", "NA_HETLOSS", "NA_SNV", "NA_NA"),
                     
                     "col"=c("#00441b", "#00441b", "#00441b", "#238b45", "#74c476",
                            "#a50f15", "#a50f15", "#a50f15", "#a50f15", "#fb6a4a",
                            "#08519c", "#08519c", "#08519c", "#6baed6", "grey"))

range.stdres.xyz <- compTwoGenes(gene.x, gene.y, col.df)
range.stdres.xyz[which(range.stdres.xyz$UID == 'HETLOSS_HETLOSS'),]

range.stdres.xyz <- compTwoGenes('ADAM10', 'NOTCH1', col.df, cex.type = 'xy', plot.legend=FALSE)
head(range.stdres.xyz[which(range.stdres.xyz$seg.y < -0.3),], 10)
#range.stdres.xyz[which(range.stdres.xyz$UID == 'HETLOSS_HETLOSS'),]

range.stdres.xyz <- compTwoGenes('AJUBA', 'NOTCH1', col.df, cex.type = 'xy', plot.legend=FALSE)
head(range.stdres.xyz[which(range.stdres.xyz$seg.y < -0.3),], 10)


