#' Calculate mutations and their probabilities
#' 
#' Calculate the number of mutations and the probabilities on this number of mutations per cells. 
#' 
#' @param cell_types Character vector indicating the cell types
#' @param genes Character vector indicating which genes to analyse
#' @param stat_types Character vector indicating the types of statistic to calculate:
#' "mean", "mean-std" and/or "mean+std"
#' @param verbose Boolean indicating verbose mode
#'
#' @return Data.frame with the number of cells, the number of mutations, the probabilities of 
#' the mutations and the input variables
#' 
#' @export
#'
#' @examples
#' probDriverMut(cell_type = "iPS", gene = "all", stat_type = "mean")
#' 

probDriverMut <- function(cell_types, genes, stat_type = "mean", verbose = F)
{
  results1 <- foreach(cell_type = cell_types, .combine='rbind') %:% 
    foreach(gene = genes, .combine='rbind') %:%
      foreach(stat_type = stat_type, .combine='rbind') %do% { 
        string <- sprintf('%s_%s_%s', cell_type, gene, stat_type)
        if(verbose)
          print(string)
    
        gene = driver_counts[gene]
        
        d_time <- doubling_time[[cell_type]]
        if (stat_type == 'mean'){
          m_rate <- unlist(mutation_rate[[cell_type]])[1]
        } else if (stat_type == 'mean-std'){
          m_rate <- unlist(mutation_rate[[cell_type]])[1] - unlist(mutation_rate[[cell_type]])[2]
        } else if (stat_type == 'mean+std'){
          m_rate <- unlist(mutation_rate[[cell_type]])[1] + unlist(mutation_rate[[cell_type]])[2]
        }
        
        cell_numbers <- foreach(i=seq(0,t_culture,mitosis_rate[[cell_type]]), .combine='c') %do% {
          calcCellNumb(N0 = init_numb_cells, d = d_time, t = i)
        }
        
        if(verbose)
          print("Cell numbers: DONE")
        
        mut = 0
        mutations <- foreach(n=seq(0,t_culture/mitosis_rate[[cell_type]],1), .combine='c') %do% {
          m <- calcMutations(N=n, N0=init_numb_cells, mu=m_rate, d=d_time, mt=mitosis_rate[[cell_type]],
                              dep=CDS_dep[[cell_type]], npos=6*CDS_length[[cell_type]], 
                              perc=CDS_perc, psnv=probs_snv[[cell_type]], count=gene[[1]])
          mut <- mut+m
          # print(m)
          return(mut)
        }
        
        if(verbose)
          print("Mutations: DONE")
        
        p=0
        probabilities <- foreach(n=seq(0,t_culture/mitosis_rate[[cell_type]],1), .combine='c') %do% {
          # print(c(p, 1-p))
          if (1-p < 1e-9){
            # print(T)
            return(1)
          } else {
            p <- calcProbabilities(N=n, N0=init_numb_cells, mu=m_rate, d=d_time, mt=mitosis_rate[[cell_type]],
                                   dep=CDS_dep[[cell_type]], npos=CDS_length[[cell_type]], 
                                   perc=CDS_perc, psnv=probs_snv[[cell_type]], count=gene[[1]])
            
            return(p)
          }
            
        }
        
        if(verbose)
          print("Probabilities: DONE")
        
        results <- data.frame(cells=cell_numbers, muts=mutations, probs=probabilities)
        results$cell_type <- rep(cell_type, length(cell_numbers))
        results$gene <- rep(names(gene), length(cell_numbers))
        results$stat <- rep(stat_type, length(cell_numbers))
        
        #f <- strsplit(args[2], split="/")[[1]]
        #f <- c(f[1:length(f)-1], paste0("results_", args[1], ".txt"))
        #write.table(results, file = paste(f, collapse = "/"))
        return(results)
      }
  return(results1)
}