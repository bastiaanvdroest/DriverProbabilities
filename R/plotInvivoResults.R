#' Plot in vivo results
#'
#' Function to plot \(a subset of\) the number of in vitro grown SI ASC, with the
#' probability on oncogenic mutation, compared to number of years of in vivo 
#' mutation accumulation of SI ASC.
#'
#' @param r Dataframe with results from \code{\link{probDriverMut}}
#' @param gene Character string to indicate which gene must be plotted
#' @param color Character string for color of line in figure
#' @param xlab Character string with label for x-axis. Default is\cr
#' "Number of cells"
#' @param ylab Character string with label for y-axis. Default is\cr
#' "Number of SNVs in oncogenes"
#' @param title Character string with title of figure. Default is no title
#' @param pdf_out Character with file name where to save the figure. If not given,
#' figure will not be saved, only returned from function
#'
#' @return Plot of probablity on oncogenic mutations against cell numbers and years
#' of in vivo mutation accumulation
#' @export
#' @import ggplot2
#' @import gridExtra
#'
#' @examples
#' r = data.frame(cells = calcCellNumb, muts = calcMutEvents, probs = probDriverMut)
#' 
#' invivo_plot = plotInvivoResults(r, gene = "all", color = "red")
#' 

plotInvivoResults <- function(r,
                              cell_type = "SI",
                              gene ="all",
                              color = "black",
                              xlab = "Probability",
                              ylab = "Number of cells",
                              title = "",
                              pdf_out)
{

  if (!(gene %in% r$gene))
    stop("Given gene name not found in 'r'")
  
  if (gene == "all"){
    invivo_plot <- ggplot(r[r$cell_type == 'SI' & r$gene=='all' & r$year > 1e-5 & r$year < 1,], aes(x=probs, y=cells))+
      geom_point(color = color) +
      geom_line(size=2, color = color) +
      scale_y_log10(sec.axis = sec_axis(~f_prob_all(.), name="Years of adult life", breaks=c(c(1,5) %o% 10^(-4:0))), breaks = c(1 %o% 10^(4:7))) +
      labs(title=title, x = xlab, y=ylab)
  } else if (gene == "BRAF.V600E"){
    invivo_plot <- ggplot(r[r$cell_type == 'SI' & r$gene=='BRAF.V600E' &  r$year < 100,], aes(x=probs, y=cells))+
      geom_point(color = color) +
      geom_line(size=2, color = color) +
      scale_y_continuous(expand = c(0,0), sec.axis = sec_axis(~f_prob(.), name="Years of adult life" )) + #, breaks=c(seq(1,9,2) %o% 10^(0:7))), breaks=c(seq(1,9,2) %o% 10^(7:11))) +
      labs(title=title, x = xlab, y=ylab)
  }
  
  invivo_plot <- invivo_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               text = element_text(size=20))
  
  if(missing(pdf_out)){
    return(invivo_plot)
  } else {
    pdf(pdf_out, width=15, height=10)
    print(invivo_plot)
    dev.off()
   return(invivo_plot)
  }
}