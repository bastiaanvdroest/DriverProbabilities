#' Plot results per cell
#'
#' Function to plot \(a subset of\) the results from \code{\link{calcCellNumb}}, \code{\link{calcMutEvents}} and 
#' \code{\link{probDriverMut}}. Unit is cells
#'
#' @param r Dataframe with results from \code{\link{probDriverMut}}
#' @param xmax Numeric with x-axis maximum
#' @param ymax Numeric with y-axis maximum
#' @param logscale (Optional) Character vector indicating which axes must be set to 
#' logscale. Default is \code{c("x","y")} 
#' @param vline Numeric vector with vertical lines
#' @param hline Numeric vector with horizontal lines
#' @param xlab Character string with label for x-axis. Default is\cr
#' "Number of cells"
#' @param ylab Character string with label for y-axis. Default is\cr
#' "Number of SNVs in oncogenes"
#' @param title Character string with title of figure. Default is\cr
#' "Mutation accumulation of driver gene mutations"
#' @param pdf_out Character with file name where to save the figure. If not given,
#' figure will not be saved, only returned from function
#'
#' @return Plot of mutational events against cell numbers.
#' @export
#' @import ggplot2
#' @import gridExtra
#'
#' @examples
#' r = data.frame(cells = calcCellNumb, muts = calcMutEvents, probs = probDriverMut)
#' 
#' mutation_plot = plotMutationsPerCell(r, xmax = 1e13, ymax = 1e4, hline = 1, vline = 1e6)
plotMutationsPerCell <- function(r, 
                                 xmax,
                                 ymax,
                                 logscale = c("x","y"),
                                 vline = F,
                                 hline = F,
                                 xlab = "Number of cells",
                                 ylab = "Number of SNVs in oncogenes",
                                 title = "",
                                 pdf_out)
{

  mut_plot <- ggplot(data=r, aes(x=cells, y=muts, colour=gene))+
    geom_line(aes(linetype=cell_type), size=1) +
    scale_linetype_manual(values=c("dotted", "dashed","solid")) +
    labs(x=xlab, y=ylab, title = title, 
         color='Type of mutations',
         linetype='Cell types') 
  
  if (length(xmax) > 0 & length(ymax) > 0){
    mut_plot <- mut_plot + coord_cartesian(xlim = c(1, xmax), ylim = c(1e-10,ymax))
  } else if (length(xmax) > 0 & missing(ymax)){
    mut_plot <- mut_plot + coord_cartesian(xlim = c(1, xmax))
  } else if (length(ymax) > 0 & missing(xmax)){
    mut_plot <- mut_plot + coord_cartesian(xlim = c(1e-10, ymax))
  }
  
  if (any(vline!=F)){
    for (i in vline){
      mut_plot <- mut_plot + geom_vline(xintercept=as.numeric(i), linetype='dotted')
    }
  }
  if (any(hline!=F)){
    for (i in hline){
      mut_plot <- mut_plot +  geom_hline(yintercept=as.numeric(i), linetype='dotted') 
    }
  }
  
  breaks_x = sort(c(1 %o% 10^(c(2,5,12)), vline))
  breaks_x = unique(breaks_x)
  labels_x = foreach(b = breaks_x, .combine='c') %do% {
    if (b == 0){
      return("0")
    } else if (log(b,10)%%1==0){
      return(formatC(b,format="e",digits=0))
    } else {
      return(formatC(b,format="e",digits=3))
    }
  }
  
  breaks_y = sort(c(1 %o% 10^(c(-8,-4,0,4)), hline))
  breaks_y = unique(breaks_y)
  labels_y = foreach(b = breaks_y, .combine='c') %do% {
    if (b == 0){
      return("0")
    } else if (log(b,10)%%1==0){
      return(formatC(b,format="e",digits=0))
    } else {
      return(formatC(b,format="e",digits=3))
    }
  }
  
  if ("x" %in% logscale)
    mut_plot <- mut_plot + scale_x_log10(breaks = breaks_x, labels = labels_x)
  if ("y" %in% logscale)
    mut_plot <- mut_plot + scale_y_log10(breaks = breaks_y, labels = labels_y)

  mut_plot <- mut_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"),
                               text = element_text(size=20))
  if (any(vline != F)){
    mut_plot <- mut_plot + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  }
  
  if(missing(pdf_out)){
    return(mut_plot)
  } else {
    pdf(file=pdf_out, width=15, height=10)
    print(mut_plot)
    dev.off()
    
    return(mut_plot)
  }
}
