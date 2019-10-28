#' Plot in vivo results
#'
#' Function to plot \(a subset of\) the number of in vitro grown SI ASC, with the
#' probability on oncogenic mutation, compared to number of years of in vivo 
#' mutation accumulation of SI ASC.
#'
#' @param r Dataframe with results from \code{\link{probDriverMut}}
#' @param gene Character string to indicate which gene must be analyzed
#' @param dep_invitro List with enrichment/depletion test results of invitro experiments
#' @param CDS_length Numeric with length of Coding Sequence of 'in vivo' cell type
#'
#' @return Extended version of results from \code{\link{probDriverMut}}, with
#' in column year the amount of years in vivo
#' 
#' @export
#'
#' @examples
#' r = probDriverMut(cell_type = cell_types, gene = "all")
#' 
#' invivo_plot = plotInvivoResults(r, gene = "all", dep_invitro = CDS_dep, 
#'                                 CDS_length = CDS_length[["SI"]])
#' 

invivoYears <- function(r,
                        gene = "all",
                        dep_invitro,
                        CDS_length)
{
  
  dep_invitro <- CDS_dep
  
  total <- vector()
  for (i in 1:6){
    #print(probs_snv$SI[i]*driver_counts$all[i])
    total[i] <- ((CDS_length - driver_counts$all[i])/CDS_length)^probs_snv$SI[i]
  }
  total_prod_vitro <- prod(total)
  
  total <- vector()
  for (i in 1:6){
    #print(probs_snv$SI[i]*driver_counts$all[i])
    total[i] <- ((CDS_length - driver_counts$all[i])/CDS_length)^probs_snv$In_vivo[i]
  }
  total_prod_vivo <- prod(total)
  
  year_prob <- vector()
  for (i in 1:nrow(r)){
    if (r$gene[i] == 'all'){
      year_prob[i] <- log(1-r$probs[i], base=total_prod_vitro) / (1e8*40*0.015*dep_invivo)
    } else if (r$gene[i] == 'BRAF.V600E'){
      year_prob[i] <- log(1-r$probs[i], base=(CDS_length-1)/CDS_length) / (1e8*40*0.015*dep_invivo*probs_snv$In_vivo['T>A'])
    } else {
      year_prob[i] <- 0
    }
  }
  
  # year_mut <- vector()
  # for (i in 1:nrow(r)){
  #   year_mut[i] <- (r$muts[i] * CDS_length[['SI']])/(4e9*probs_snv$In_vivo["T>A"]*0.015*dep_invivo)
  # }
  # 
  # f_mut <- Vectorize(function(x){
  #   if (x <= 0) return(0)
  #   m <- (x * dep_invitro$SI * mutation_rate$SI[1] * probs_snv$SI['T>A']) / (4e9 * probs_snv$In_vivo['T>A'] * dep_invivo * CDS_length[['SI']]/CDS_length[['SI']])
  # })
  
  f_prob <- Vectorize(function(x){
    if (x <= 0) return(0)
    t <- (dep_invitro$SI*mutation_rate$SI[1]*probs_snv$SI['T>A']*x)/(4e9*probs_snv$In_vivo['T>A']*dep_invivo)
  })
  
  f_prob_all <- Vectorize(function(x){
    if (x <= 0) return(0)
    t <- (log(total_prod_vivo, base=total_prod_vitro)*dep_invitro$SI*mutation_rate$SI[1]*x)/(4e9*dep_invivo)
  })
  
  r$year <- year_prob
  
  f_prob <<- f_prob
  f_prob_all <<- f_prob_all
  
  return(r)
}