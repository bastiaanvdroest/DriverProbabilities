#' Calculate mutations and their probabilities
#' 
#' Calculate the number of mutations and the probabilities on this number of mutations per cells. 
#' 
#' @param cell_type Character indicating the cell type
#' @param gene Character indicating which gene to analyse
#' @param stat_type Character indicating the type of statistic to calculate, must be "mean",
#' "mean-std" or "mean+std"
#'
#' @return Data.frame with the number of cells, the number of mutations, the probabilities of 
#' the mutations and the input variables
#' 
#' @export
#'
#' @examples
#' probDriverMut(cell_type = "iPS", gene = "all", stat_type = "mean")
#' 

getCDSLength <- function(table1, table2){
  table1 <- subset(table1, region=='CDS')
  types <- unique(table2$Type)
  type_lengths <- foreach(t = types) %do% {
    group <- which(table2$Type == t)
    lengths <- foreach(i = group, .combine='c') %do% {
      s <- table2$Sample[i]
      if (s %in% table1$sample){
        length <- table1$surveyed_region_length[which(table1$sample==s)]
      }
    }
  }
  mean_types <- sapply(type_lengths,mean)
  names(mean_types) <- cell_types
  return(mean_types)
}
