#' Calculate probability of obtaining driver gene mutation
#'
#' @param N Integer giving the generation number
#' @param N0 Integer giving the initial number of cells
#' @param mu Integer giving the mutation rate of the cell type
#' @param d Integer giving the doubling time of the cell type
#' @param mt Integer giving the mitosis rate of the cell type
#' @param dep Integer giving the depletion of mutations in the Coding Sequence
#' @param npos Integer giving the length of the Coding Sequence
#' @param perc Integer giving the percentage of coding sequence in whole genome
#' @param psnv Table with the mutation spectrum from the cell type
#' @param count Table with the number of driver gene mutations specified per mutation type
#'
#' @return Integer probability of obtaining at least one driver gene mutation
#' @export
#'
calcMutations <- function(N, N0, mu, d, mt, dep, npos, perc, psnv, count){
                            
  if (N == 0){
    return(0)
  } else {
    m <- perc * dep *  mu * N0*exp(log(2)*1/d*mt*N)  #Number of mutations in CDS
    
    nmut <- foreach(t = c('C>A','C>G','C>T','T>A','T>C','T>G'),.combine='c') %do% {
      if (t %in% names(count)){
        return(m*psnv[t]*count[t]/npos)
      } else {
        return(0)
      }
    }
    
    return(sum(nmut))
  }
}
