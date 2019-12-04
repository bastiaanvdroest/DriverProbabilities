#' Calculate mutations and their probabilities
#' 
#' Calculate the number of mutations and the probabilities on this number of mutations per cells. 
#' 
#' @param r Data.frame with the results from \code{\link{probDriverMut}}
#' @param hline Numeric indicating the value on the y-axis to intersect with
#'
#' @return Numeric vector of intersections with \code{hline} for different cell types and genes
#' 
#' @export
#'
#' @examples
#' invitro_res = probDriverMut(cell_type = cell_types, gene = names(driver_counts), stat_type = "mean")
#' hline = 1
#' calcIntersect(r = invitro_res, hline = hline)
#' 

calcIntersects <- function(r, hline)
{
  res <- foreach(i = unique(r$cell_type), .combine = "c") %:%
    foreach(j = unique(r$gene), .combine = "c") %do% {
      scale = subset(r, cell_type == i & gene == j)$muts / subset(r, cell_type == i & gene == j)$cells
      return(hline/scale[length(scale)])
    }
}