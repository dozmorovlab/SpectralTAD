#' Contact matrix from Rao 2014, chromosome 20, 25kb resolution
#'
#' A sparse 3-column contact matrix
#'
#' @format A data.frame with 3 columns and 2125980 rows:
#' \describe{
#' \item{V1}{The genomic loci corresponding to a given row of the
#'  contact matrix}
#' \item{V2}{The genomic loci corresponding to a given column of the
#'  contact matrix}
#' \item{V3}{Number of contacts between Loci1 and Loci2}
#' }
#' 
#' @source Data from Rao SS, Huntley MH, Durand NC, Stamenova EK et al. 
#' A 3D map of the human genome at kilobase resolution reveals principles
#' of chromatin looping. Cell 2014 Dec 18;159(7):1665-80. PMID: 25497547. 
#' Available at \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525}
#' @return A data.frame
#' @usage data(rao_chr20_25_rep)
#'
"rao_chr20_25_rep"
