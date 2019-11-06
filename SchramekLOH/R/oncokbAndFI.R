#### OncoKB Functions ####
#' mapFIandTruncatingMutations
#' @description Maps the Functional Impact (FI) and the "Truncating Mutations"
#' to a mutation MAF data table
#' 
#' @param muts 
#' @param verbose 
#' @param FI 
#' @param truncating.ids 
#'
#' @return
#' @export
#'
#' @examples
mapFIandTruncatingMutations <- function(muts, verbose=T, FI,
                                        truncating.ids=c('Frame_Shift_Del', 'Frame_Shift_Ins', 
                                                         'Nonsense_Mutation', 'Splice_Region', 'Splice_Site')){
  if(any(grepl('Mutation Effect', colnames(muts)))) muts <- muts[,-grep('Mutation Effect', colnames(muts))]
  
  ## Assign mutation types to "Truncating Mutations" classification based on Variant_Classification column
  muts$HGVSp_Short2 <- 'NA'
  trunc.idx <- which(muts$Variant_Classification %in% truncating.ids)
  muts[trunc.idx,]$HGVSp_Short2 <- 'Truncating Mutations' 
  if(verbose) print(paste0("Number of truncating mutations: ", length(trunc.idx)))
  
  ## Assign the FI predicted functional status
  ## FI.verbose was computed in the FI Analysis section
  FI <- annotateFI(FI)
  nookb.fi <- merge(muts, FI, 
                    by.x=c('Hugo_Symbol', 'HGVSp_Short'), 
                    by.y=c('gene', 'mut'), all.x=T)
  
  for(del.ids in c("^Deletion$", "^Het. Deletion$", "^Hom. Deletion$")){
    if(any(grepl(del.ids, nookb.fi$HGVSp_Short))){
      nookb.fi[grep(del.ids, nookb.fi$HGVSp_Short),]$p.lof <- TRUE
    }
  }
  nookb.fi[is.na(nookb.fi$p.lof),]$p.lof <- FALSE
  nookb.fi[is.na(nookb.fi$c.lof),]$c.lof <- 'unknown'
  nookb.fi
}

.getMutDf <- function(merge.type, not=FALSE){
  if(not){
    merge.type <- paste0("not.", merge.type)
  }
  switch(merge.type,
         'annotated'=annotated.mut.ext,
         'actionable'=actionable.mut.ext,
         'not.annotated'=not.annotated.mut.ext,
         'not.actionable'=not.actionable.mut.ext)
}

.mergeMutToOncokb <- function(mut.to.merge, mk.uniq=TRUE){
  ## Merge genes with ACTIONABLE or ANNOTATED mutations with VARIANT STATUS oncokb
  actionable.mut.status <- merge(mut.to.merge, oncokb.var[,oncokb.var.cols], 
                                 all.x=TRUE, by.x=mut.merge.cols, by.y=oncokb.merge.cols)
  actionable.mut.status$'Mutation Effect'[is.na(actionable.mut.status$'Mutation Effect')] <- 'Not-in-OncoKB'
  
  if(mk.uniq){
    actionable.mut.status <- unique(actionable.mut.status[,c('Hugo_Symbol', 'HGVSp_Short', 
                                                             'Variant_Classification', 'Mutation Effect')])
  }
  actionable.mut.status
}

#### Functional Impact Functions ####
#' annotateFI
#' @description  Annotates the functional impact dataset based on number of tools supporting LOF
#' @param FI 
#'
#' @return
#' @export
#'
#' @examples
annotateFI <- function(FI){
  FI$p.lof <- apply(FI[,3:5], 1, function(x){
    ##VERY IMPORTANT: DEFINITION OF PREDICTED FUNCTIONAL IMPACT!!!
    sum(x %in% c('benign', 'tolerated', 'undefined', 'unknown')) <= 2
  })
  
  ## Compare against all counts of benigng/LOF
  FI.cnts <- sapply(c(0:3), function(cnt){
    apply(FI[,3:5], 1, function(x){
      sum(x %in% c('benign', 'tolerated', 'undefined', 'unknown')) == cnt
    })
  })
  colnames(FI.cnts) <- c('all.loss', 'two.loss', 'one.loss', 'all.benign')
  FI$c.lof <- colnames(FI.cnts)[apply(FI.cnts, 1, which)]
  FI
}

#' genFIdataframe
#' @description Wrapper ot get all the functional impact data and assemble into one data structure
#' 
#' @param gene.dat 
#'
#' @return
#' @export
#'
#' @examples
genFIdataframe <- function(gene.dat, ...){
  ## Parse cBio data and assemble into FI matrices
  gene.fi <- lapply(gene.dat, SchramekLOH:::parseFunctionalImpact)
  gene.fi.mats  <- lapply(gene.fi, function(gene, ...) SchramekLOH:::mkFIMatrix(gene, ...), ...)
  
  ## Assemble all FI matrices together into one large matrix
  FI.verbose <- as.data.frame(do.call(rbind, lapply(gene.fi.mats, function(i) i[[1]])))
  FI.verbose$gene <- gsub("\\.[0-9]*$", "", rownames(FI.verbose))
  FI.verbose <- FI.verbose[,c(5,1,2,3,4)]
  
  FI.score <- as.data.frame(do.call(rbind, lapply(gene.fi.mats, function(i) i[[2]])))
  FI.score$gene <- gsub("\\.[0-9]*$", "", rownames(FI.score))
  FI.score <- FI.score[,c(5,1,2,3,4)]
  
  ## Return all data structures
  list("FI"=gene.fi, "FI.mat"=gene.fi.mats, "FI.verbose"=FI.verbose, "FI.score"=FI.score)
}

## Fills 'undefined' for all errors
.genTrunc <- function(){
  undefined <- paste(c("MutationAssessor:", "SIFT:", "Polyphen-2:"), 
                     "impact: truncating, score: -1")
  paste(undefined, collapse=";")
}

## Given a gene-specific table of mutational data extracted from cBioportal, parse the 'Functional Impact' column for each gene
parseFunctionalImpact <- function(g){
  non.missense.idx <- which(!g$'Mutation Type' %in% c('Missense_Mutation', 'Splice_Region', 'Splice_Site', 'Translation_Start_Site'))
  if(length(non.missense.idx) > 0) g[non.missense.idx,]$'Functional Impact' <- .genTrunc()
  
  idx <- which(g$'Functional Impact' == '')
  if(length(idx) > 0) g <- g[-idx,]
  
  #if(length(idx) > 0) g <- g[-idx,]
  
  if(nrow(g) > 0){
    na.idx <- is.na(g$`Functional Impact`)
    if(any(na.idx)) {
      g$`Functional Impact`[which(na.idx)] <- SchramekLOH:::.genUndefined()
    }
    impacts <- strsplit(g$'Functional Impact', split=";")
    impacts <- lapply(impacts, function(i) gsub("Error", "impact: undefined, score: undefined", i))
    impact <- sapply(impacts, function(i) gsub("^.*impact: ", "", i) %>% gsub(",.*", "", .))
    score <- suppressWarnings(sapply(impacts, function(i) as.numeric(gsub("^.*score: ", "", i))))
    ids <- sapply(impacts, function(i) gsub(":.*", "", i))
    
    data.frame('mut'=as.character(matrix(rep(g$'Protein Change',3), nrow=3, byrow=TRUE)),
               "tool"=as.character(ids),
               "impact"=as.character(impact),
               "score"=as.numeric(score))
  } else {
    NULL
  }
  
}

## Fills 'undefined' for all errors
.genUndefined <- function(){
  undefined <- paste(c("MutationAssessor:", "SIFT:", "Polyphen-2:"), 
                     "impact: undefined, score: undefined")
  paste(undefined, collapse=";")
}

## Assemble the Functional Impact scores for each gene into a matrix for each gene and mutation
mkFIMatrix <- function(gene, mk.uniq=TRUE){
  require(reshape2)
  if(length(gene) > 0){
    lapply(c('impact', 'score'), function(v){
      if(any(gene$impact == '')){
        gene[which(gene$impact == ''),'impact'] <- 'undefined'
      }
      mut.chk <- table(unique(gene)$mut)>3
      if(any(mut.chk)){
        gene.idx <- which(gene$mut == names(which(mut.chk)))
        ma.idx <- grep("MutationAssessor", gene[gene.idx,]$tool)
        ma.mixup <- gene[gene.idx[ma.idx],]
        rm.idx <- which(ma.mixup$impact == 'undefined')
        gene <- gene[-gene.idx[ma.idx[rm.idx]],]
      }
      
      if(mk.uniq) {
        unmelt <- dcast(unique(gene), mut~tool, value.var=v)
      } else {
        tmp <- ddply(gene, .(mut, tool), transform, newid = paste(mut, seq_along(tool)))
        unmelt <- dcast(tmp, mut + newid ~ tool, value.var = v)
        unmelt <- unmelt[,-which(colnames(unmelt) == "newid")]
      }
      unmelt
    })
  }
}

#' addCnvToFi
#'
#' @param g 
#'
#' @return
#' @export
#'
#' @examples
addCnvToFi <- function(g, gene.fi.spl, all.lt.genes,
                       cnvs=c('Het. Deletion', 'Hom. Deletion', 'Amplification')) {
  i <- all.lt.genes[[g]]
  cnv.ids <- rev(sort(i$HGVSp_Short[i$HGVSp_Short %in% cnvs]))
  gene.fi.spl[[g]][,1] <- as.character(gene.fi.spl[[g]][,1])
  
  cnv.mat <- as.data.frame(matrix(rep(as.character(cnv.ids), 4), ncol=4))
  colnames(cnv.mat) <- colnames(gene.fi.spl[[g]][,2:5])
  
  rbind(gene.fi.spl[[g]][,-1], cnv.mat)
}

## Visualization of the FI matrices by first inistializing a blank plot, and then adding rectangles representing the different FI states
## Makes a blank plot to plot FI 
mkBlankPlot <- function(gene.fi.spl, xmax=50, fic){
  par(mar=c(10, 10, 4.1, 2.1))
  plot(0, type='n', xlim=c(1,xmax), ylim=c(0.5, (length(gene.fi.spl) + 0.5)), 
       yaxt='n', xaxt='n', ylab="", xlab="Mutations")
  axis(side=2, at=c(1:length(gene.fi.spl)), labels=names(gene.fi.spl), line = 5, 
       tick = FALSE, cex.axis=0.9, las=1)
  fi.idx <- c(0.7, 1, 1.3)
  axis(side=2, line = 0, las=2, cex.axis=0.7,
       at=as.numeric(sapply((seq_along(gene.fi.spl)-1), function(f) f + fi.idx)), 
       labels=rep(colnames(gene.fi.spl[[1]][,2:4]), length(gene.fi.spl)))
  legend("bottomright", fill=fic, legend=names(fic), box.lwd = 0)
}

## Visualization of the FI matrices by first inistializing a blank plot, and then adding rectangles representing the different FI states
## Adds rectangles for FI data
addRect <- function(fim, idx, y.cent=1, fi.col){
  if(any(fim=='')) fim[fim==''] <- 'undefined'
  y.bot <- switch(idx,
                  '1'=0.55,
                  '2'=0.85,
                  '3'=1.15)
  y.bot <- y.bot + (y.cent - 1)
  y.top <- y.bot + 0.3
  rect(xleft=0:(nrow(fim)-1), ybottom=rep(y.bot, nrow(fim)), 
       xright = 1:nrow(fim), ytop = rep(y.top, nrow(fim)),  
       col=as.character(fi.col[as.character(fim[,idx])]), border='white')
}