#' Counting driver activating mutations
#' 
#' Counts the number of driver activating mutations in a gene or in total
#' 
#' @param drivers Dataframe with 2 columns. First column must be gene ID and second the mutations
#' @param per_gene Indicating if driver mutations are counted per gene, default is FALSE
#'
#' @return Table of counts of the 6 mutation types 'C>A', 'C>G, 'C>T, 'T>A', 'T>C', 'T>G'
#' @export
#' @import foreach
#' @importFrom Biostrings complement
#' @importFrom Biostrings DNAStringSet
#'
#' @examples
#' drivers = data.frame(gene=c('ABL1','ABL1'),gdna=c('G>A','A>G'), stringsAsFactors=FALSE)
#' countDriverSNVs(drivers)
countDriverSNVs <- function(drivers, per_gene = F, tri = F){
  if (per_gene == T){
    SNVS <- foreach(g = unique(drivers[,1])) %:%
      foreach(x = drivers[drivers[,1]==g,2], .combine='c') %do% {
        if (grepl('__', x)){
          x <- unlist(strsplit(x, "__"))
        }
        if (all(endsWith(x, 'A') | endsWith(x, "C") | endsWith(x, "G") | endsWith(x, "T"))){
          y <- substr(x, nchar(x)-2, nchar(x))
          for (i in 1:length(y)){
            y_split <- unlist(strsplit(y[i], '>'))
            if (y_split[1] == "G" | y_split[1] == "A"){
              y_split[1] <- as.character(complement(DNAStringSet(y_split[1])))
              y_split[2] <- as.character(complement(DNAStringSet(y_split[2])))
            }
            y[i] <- paste(y_split, collapse='>')
          }
          return(y)
        }
      }
    tables <- foreach(g = 1:length(SNVS)) %do% table(SNVS[[g]])
    names(tables) <- unique(drivers[,1])
    return(tables)
  } else if (per_gene == F & tri == T){
    SNVS <- foreach(x = drivers[,2], .combine = "c") %do% {
      if (grepl('__', x)){
        x <- unlist(strsplit(x, "__"))
      }
      if (all(endsWith(x, 'A') | endsWith(x, "C") | endsWith(x, "G") | endsWith(x, "T"))){
        chromosome <- do.call(rbind, strsplit(x, ":g."))[,1]
        mut <- do.call(rbind, strsplit(x, ":g."))[,2]
        position <- as.numeric(substr(mut, 1, nchar(mut)-3))
        y <- substr(mut, nchar(mut)-2, nchar(mut))
        
        for (i in 1:length(y)){
          y_split <- unlist(strsplit(y[i], '>'))
          tn <- as.character(getSeq(Hsapiens, chromosome[i], position[i]-1, position[i]+1))
          if (y_split[1] == "G" | y_split[1] == "A"){
            y_split[1] <- as.character(complement(DNAStringSet(y_split[1])))
            y_split[2] <- as.character(complement(DNAStringSet(y_split[2])))
            tn <- as.character(reverseComplement(DNAStringSet(tn)))
          }
          y[i] <- sprintf("%s[%s>%s]%s", substr(tn,1,1), y_split[1], y_split[2], substr(tn,3,3))
        }
        return(y)
      }
    }
    return(list(all=table(SNVS)))
  }  else if (per_gene == F & tri == F){
    SNVS <- foreach(x = drivers[,2], .combine='c') %do% {
      if (grepl('__', x)){
        x <- unlist(strsplit(x, "__"))
      }
      if (all(endsWith(x, 'A') | endsWith(x, "C") | endsWith(x, "G") | endsWith(x, "T"))){
        y <- substr(x, nchar(x)-2, nchar(x))
        for (i in 1:length(y)){
          y_split <- unlist(strsplit(y[i], '>'))
          if (y_split[1] == "G" | y_split[1] == "A"){
            y_split[1] <- as.character(complement(DNAStringSet(y_split[1])))
            y_split[2] <- as.character(complement(DNAStringSet(y_split[2])))
          }
          y[i] <- paste(y_split, collapse='>')
        }
        return(y)
      }
    }
    return(list(all=table(SNVS)))
  } else {
    stop(sprintf('Cannot count the driver activating mutations. per_gene must be TRUE or FALSE'))
  }
}
