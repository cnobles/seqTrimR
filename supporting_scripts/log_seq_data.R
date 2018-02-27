#' Return general information about sequence data objects during processing
#' 
#' @param seqs ShortReadQ object of reads or unique sequences
#' @author Christopher Nobles, Ph.D.
log_seq_data <- function(seqs){
  require("BiocGenerics")
  require("ShortRead")
  require("pander")
  
  if(class(seqs) != "list"){
    cnt <- length(seqs)
    nts <- sum(width(seqs))
    ave <- mean(width(seqs))
    tbl <- data.frame(
      "Reads" = c("counts", format(cnt, big.mark = ",")),
      "Bases" = c("counts", format(nts, big.mark = ",")),
      "Ave.Length" = c("nts", round(ave, digits = 1)),
      row.names = NULL)
    return(tbl)
  }else{
    cnt <- sapply(seqs, length)
    nts <- sapply(seqs, function(x) sum(width(x)))
    ave <- sapply(seqs, function(x) mean(width(x)))
    tbl <- data.frame(
      "Reads" = c("counts", format(cnt, big.mark = ",")),
      "Bases" = c("counts", format(nts, big.mark = ",")),
      "Ave.Length" = c("nts", round(ave, digits = 1)),
      row.names = NULL)
    return(tbl)
  }
}