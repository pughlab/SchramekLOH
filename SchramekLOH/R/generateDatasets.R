
#' generateDatasets
#'
#' @param git Path to the git/SchramekLOH/data directory
#'
#' @return
#' @export
#'
#' @examples
generateDatasets <- function(git="~/git/SchramekLOH/SchramekLOH/data"){
  # Bed file obtained form http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/GenomeWideSNP_6.hg19.bed.zip
  pdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2'
  ref <- file.path(pdir, "ref")
  snp6 <- read.table(file.path(ref, "GenomeWideSNP_6.hg19.headerless.bed"),
                     header=FALSE, sep="\t", stringsAsFactors = FALSE,
                     check.names = FALSE)
  save(snp6, file=file.path(git, "snp6.Rdata"))

  ## Mapping of Birdseed genotype files
  mapping <- read.table(file.path(ref, "mappingIds.txt"),
                        header=FALSE, sep="\t", stringsAsFactors = FALSE,
                        check.names = FALSE)
  save(mapping, file=file.path(git, "mapping.Rdata"))

  ## Mapping of expression files
  expr.mapping <- read.table(file.path(pdir, "gene_expression", "ref", "mappingIds.txt"),
                             header=FALSE, sep="\t", stringsAsFactors = FALSE,
                             check.names = FALSE)
  save(expr.mapping, file=file.path(git, "expr.mapping.Rdata"))

  ## TCGA hg19 GAF file of knownGenes and locations
  gaf <- read.table(file.path(pdir, "gene_expression", "ref", "TCGA.hg19.June2011.gaf"),
                    header=FALSE, sep="\t", stringsAsFactors = FALSE,
                    check.names = FALSE)
  save(gaf, file=file.path(git, "gaf.Rdata"))


  ## Reading in the seg files
  outdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/summary_af'
  vcfdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/symlink'
  segdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc'
  exprdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/gene_expression/symlink'
  affysegdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/cnv_seg/symlink'

  seg <- read.table(file.path(segdir, "hnsc_tcga_segments.seg"),
                    header=TRUE, stringsAsFactors = FALSE,
                    check.names=FALSE)
  save(seg, file=file.path(git, "TCGAseg.Rdata"))

  affyseg <- read.table(file.path(affysegdir, "TCGA_HNSC_affy6.seg"),
                        header=TRUE, stringsAsFactors = FALSE,
                        check.names = FALSE)
  save(affyseg, file=file.path(git, "Affyseg.Rdata"))
  NULL
}

#' loadBirdseed
#' @description This is controlled access, so the resulting Rdata file will
#'   not be committed and psuehd in the git repo
#' @return
#'
#' @examples
loadBirdseed <- function(git="~/git/SchramekLOH/SchramekLOH/data"){
  vcfdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/symlink'

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

  save(df.bs, mapping, ref.probe.ord,
       mapping.cov, file=file.path(git, "birdseed.Rdata"))
}

#' loadGeneExpr
#' @description This is open access, so the resulting Rdata file is pushed
#' @return
#'
#' @examples
loadGeneExpr <- function(git="~/git/SchramekLOH/SchramekLOH/data"){
  exprdir <- '/Users/rquevedo/Onedrive/PughLab/Schramek_notch/IGV_segs/TCGA_hnsc_vcf2/gene_expression/symlink'

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

  save(df.ex, file=file.path(git, "geneExpr.Rdata"))
  NULL
}
