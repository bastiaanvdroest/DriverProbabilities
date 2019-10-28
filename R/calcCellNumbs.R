#' Calculate cell numbers
#' 
#' Calculates the number of cells after \code{t} hours, by using the 
#' doubling time as a measure of mitosis and death of cells in the culture
#' 
#' The function calculates the number of cells by using the equation 
#' \eqn{N = N0 * 2^(t/d)}
#'
#' @param N0 Initial number of cells. Default = 1
#' @param d Doubling time of the cell culture. Default = 24 hours
#' @param t Time in hours
#'
#' @return The number of cells after \code{t} hours
#' @export
#'
#' @examples
#' calcCellNumber(1,24,24) = 2
calcCellNumb <- function(N0 = 1, d = 24, t){
  d = log(2)*(1/d)   # Doubling time adjusted to calculate with the power of 2
  N = N0*exp(d*t)
  return(N)
}