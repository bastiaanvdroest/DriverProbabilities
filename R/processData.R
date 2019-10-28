#' Process incoming data
#' 
#' Read the data in the given files and process it to data.frames or vector needed by 
#' downstream functions. Processing of drivers is optional, if 'process_drivers' is set to False
#' a list of driver mutations must be given later on. 
#' 
#' @param mutation_rate_file Character indicating the file containing column 'Cell_type' and
#' 'Number_of_mutations_per_cell_cycle'
#' @param oncogenic_mutations_file Character indicating the file containing a oncogene driver
#' catalog with columns 'gene' and 'gdna' containing the variants
#' @param gene_MoA_file Character indicating the file containing columns 'gene' and 'gene_MoA'.
#' The column 'gene_MoA' contains either 'Act', 'ambiguous' or 'LoF' as mode of action of gene
#' @param prob_per_gene Boolean indicating if probability must be calculated per gene or for 
#' all genes together
#' @param trinucleotide Boolean indicating if trinucleotide context of drivers must be computed
#' @param CDS_dep_invitro_file Character indicating the file with the enrichment/depletion test
#' results from the mutationalPatterns package for invitro data
#' @param CDS_dep_invivo_file Character indicating the file with the enrichment/depletion test
#' results from the mutationalPatterns package for invivo data
#' @param sample_cell_type_file Character indicating the file with the 2 colums 'Sample' and 'Type'
#' @param invivo_tissue Character indicating the invivo tissue
#' @param mut_spec_invitro_file Character indicating the file with the invitro mutational spectrum
#' @param mut_spec_invitro_file Character indicating the file with the invivo mutational spectrum
#'
#' @return All necessary data for processing downstream will be imported as global variables
#'

processData <- function(mutation_rate_file,
                        process_drivers = T,
                        oncogenic_mutations_file,
                        gene_MoA_file,
                        prob_per_gene = F,
                        trinucleotide = F,
                        CDS_dep_invitro_file,
                        CDS_dep_invivo_file,
                        sample_celltypes_file,
                        invivo_tissue,
                        mut_spec_invitro_file,
                        mut_spec_invivo_file){
  
  mutations <- read.table(file=mutation_rate_file, header=T)
  mutation_rate <- foreach(t = unique(mutations$Cell_type)) %do% {
    m <- mean(subset(mutations, Cell_type==t)[,3])
    s <- sd(subset(mutations, Cell_type==t)[,3])
    return(c(m,s))
  }
  names(mutation_rate) = unique(mutations$Cell_type)
  
  if (process_drivers){
    driver_table1 <- read.table(file=oncogenic_mutations_file, header=T, sep='\t')
    driver_table2 <- read.table(file=gene_MoA_file, header=T, sep='\t')
    
    drivers <- getDriverSNVs(driver_table1, driver_table2)
    driver_counts <- countDriverSNVs(drivers, per_gene = prob_per_gene, tri = trinucleotide)
  }
  
  CDS_table1 <- read.table(file=CDS_dep_invitro_file, header=T, sep='\t')
  CDS_table2 <- read.table(file=sample_celltypes_file, header=T, sep='\t')
  
  CDS_length <- getCDSLength(CDS_table1, CDS_table2)
  
  dep_table_invitro <- CDS_table1
  dep_table_invitro$tissue <- rep(CDS_table2$Type, 3)
  
  dep_table_invivo <- read.table(CDS_dep_file_invivo, header=T, sep=' ')
  dep_table_invivo$tissue <- rep(invivo_tissue, nrow(dep_table_invivo))
  
  dep_table <- rbind(dep_table_invitro, dep_table_invivo)
  
  dep_table_CDS <- subset(dep_table, region=='CDS')
  dep_invitro <- foreach(t = unique(dep_table_invitro$tissue)) %do% {
    dep <- dep_table_CDS[dep_table_CDS$tissue==t, 'observed']/dep_table_CDS[dep_table_CDS$tissue==t, 'expected']
    return(mean(dep))
  }
  names(dep_invitro) <- unique(dep_table_invitro$tissue)
  
  dep_invivo <- mean(dep_table_invivo[dep_table_invivo$region=="CDS",'observed']/dep_table_invivo[dep_table_invivo$region=="CDS",'expected'])
  
  probs_snv_table <- read.table(file = mut_spec_invitro_file, header=T)
  probs_snv <- foreach(type = unique(probs_snv_table$by)) %do% {
    tb <- subset(probs_snv_table, by == type)
    p <- tb$mean[order(tb$variable)]
    p[3] <- p[3]+p[4]
    p <- p[-4]
    names(p) <- unique(probs_snv_table$sub_type)
    return(p)
  }
  names(probs_snv) <- unique(probs_snv_table$by)
  
  probs_snv$In_vivo <- read.table(file = mut_spec_in_vivo_file, header=T)$mean
  names(probs_snv$In_vivo) <- unique(read.table(file = mut_spec_in_vivo_file, header=T)$sub_type)
  
  if (!prob_per_gene){
    wanted_genes <- 'all'
  }
  
  mutation_rate <<- mutation_rate
  driver_counts <<- driver_counts
  CDS_length <<- CDS_length
  dep_invivo <<- dep_invivo
  probs_snv <<- probs_snv
  CDS_dep <<- dep_invitro
  wanted_genes <<- wanted_genes
}
