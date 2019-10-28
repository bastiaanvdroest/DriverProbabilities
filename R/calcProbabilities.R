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
#' @return Number of mutations accumulated after N generations
#' @export
#'
 
calcProbabilities <- function(N, N0 = 1, mu = NULL, d = NULL, mt = NULL, dep = NULL, npos = NULL, perc = 0.015, psnv = NULL, count = NULL, invivo = F)
{
  if (N == 0){
    return(0)
  } else {
    if (invivo){
      B = 1e6
      m <- N0*mu*perc*dep*N
      ms <- m*psnv
      p <- ((npos-1)/npos + (1-psnv)*1/npos)^ms
      return(1-p)
    
    } else {
      B = 1e8
      # ps=0
      # for (i in names(count)){
      #   ps = ps + psnv[[i]]
      # }
      
      m <- perc * dep *  mu * N0*exp(log(2)*1/d*mt*N)
      # print(N)
    
    
      p <- foreach(t = c('C>A','C>G','C>T','T>A','T>C','T>G'),.combine='c') %do% {
        if (t %in% names(count)){
          ms = m*psnv[t]
          # if (ms < 100){
          pm <- ((npos-count[t])/npos)^ms
          # } else if (ms > B){
          #   pr=1
          #   n=0
          #   while(n <= ms){
          #     pr <- prod(c(pr,rep((npos-count[t] + 1-psnv[t]),B)/rep(npos, B)))
          #     n = n + B
          #   }
          #   pr <- prod(c(pr,rep((npos-count[t] + 1-psnv[t]), ms-n+B)/rep(npos, ms-n+B)))
          #   pm <- prod(pr)
          # } else {
          #   pr <- rep((npos-count[t] + 1-psnv[t]), ms)/rep(npos, ms)
          #   pr_rest <- ((npos-count[t])/npos + (1-psnv[t])*1/npos)^(ms%%1)
          #   pm <- prod(c(pr, pr_rest))
          # }
          
          return(pm)
        } else {
          return(NULL)
        }
      }
    }
    
    return(1-prod(p))
  }
}
