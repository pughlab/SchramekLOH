#' filterTcgaMuts
#'
#' @param muts 
#'
#' @return
#' @export
#'
#' @examples
#' muts <- read.table("~/git/SchramekLOH/SchramekLOH/data-raw/data_mutations_extended.trimmed.2019-06-11.txt",
#'                    quote="", fill=FALSE, stringsAsFactors = FALSE, check.names = FALSE,
#'                    sep="\t", header=TRUE)
#' x <- filterTcgaMuts(muts)
#' save(x[['muts']], x[['simple.mut']], 
#'      file="~/git/SchramekLOH/SchramekLOH/data/muts.extended.RData")
filterTcgaMuts <- function(muts){
  require(reshape2)
  idx <- grep("[*][0-9]*[*]", muts$HGVSp_Short, perl = TRUE) # Non-coding mutations
  idx <- c(idx, which(muts$HGVSp_Short == '')) # Non-coding mutations
  idx <- c(idx, grep("p.M1[?]", muts$HGVSp_Short)) # TSS: Arguable whether to remove this
  idx <- c(idx, grep("=$", muts$HGVSp_Short)) # Silent mutations
  muts <- muts[-idx,]
  
  simple.mut <- muts[,c('Tumor_Sample_Barcode',
                        'Hugo_Symbol',
                        'HGVSp_Short')]
  list('muts'=muts, 'simple.mut'=simple.mut)
}

library(reshape2)




