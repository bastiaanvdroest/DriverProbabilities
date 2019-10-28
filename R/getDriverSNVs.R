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

getDriverSNVs <- function(table1, table2){
  drivers <- data.frame(gene=table1$gene, gdna=table1$gdna)
  drivers <- foreach(i=1:nrow(drivers), .combine='rbind') %do% {
    r <- drivers[i,]
    if (r$gene %in% table2$gene){
      rg <- which(table2$gene == r$gene)
      if (table2$gene_MoA[rg] == 'Act'){
        if (!grepl('del',r$gdna) & !grepl('ins',r$gdna) & !grepl('dup',r$gdna)){
          return(r)
        }
      }
    }
  }
}
