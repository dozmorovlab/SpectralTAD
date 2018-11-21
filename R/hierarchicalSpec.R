#' Hierarchical Spectral Clustering of TADs
#'
#' Given a sparse 3 column or n x n contact matrix, hierarchicalSpec returns a list of TAD coordinates in bed form.
#' @param cont_mat Contact matrix in either sparse 3 column or n x n form. If a non-sparse matrix is used, the column names must correspond to the start point of the corresponding bin. Required.
#' @param chr The chromosome of the contact matrix being analyzed. Required.
#' @param levels The number of levels of the TAD hierarchy to be calculated, The default setting is 1.
#' @param qual_filter Option to turn on quality filtering which removes TADs with negative silhouette scores. The default setting is TRUE.
#' @param z_clust Option to filter sub-TADs based on the z-score of their eigenvector gaps. The default setting is FALSE.
#' @param eigenvalues The number of eigenvectors to be calculated. The default and suggested setting is 2.
#' @param min_size The minimum allowable TAD size measured in bins. Defaults to 5.
#' @param resolution The resolution of the contact matrix. If none selected the resolution is estimated by taking the difference in endpoints between the first and second bin.
#' @keywords spectral clustering
#' @export
#' @details
#' hierarchicalSpec()

hierarchicalSpec = function(cont_mat, chr, levels = 1, qual_filter = TRUE, z_clust = FALSE, eigenvalues = 2, min_size =5, resolution = "auto") {
  source("windowedSpec.R")

  row_test = dim(cont_mat)[1]
  col_test = dim(cont_mat)[2]

  if (col_test == 3) {
    cont_mat = HiCcompare::sparse2full(cont_mat)
    message("Converting to full matrix")
  } else if (col_test!=3 & row_test != col_test) {
    stop("Contact matrix must be sparse or square!")
    break
  }

  if (resolution == "auto") {
    resolution = as.numeric(colnames(cont_mat)[2])-as.numeric(colnames(cont_mat)[1])
  } else {
    resolution = resolution
  }

  bed = windowedSpec(cont_mat, chr = chr, resolution = resolution, z_clust = z_clust, eigenvalues = eigenvalues, min_size = min_size, qual_filter = qual_filter) %>% mutate(Level = 1)

  coords = cbind(match(bed$start, as.numeric(colnames(cont_mat))), match(bed$end-resolution, as.numeric(colnames(cont_mat))))

  tads = apply(coords, 1, function(x) cont_mat[x[1]:x[2], x[1]:x[2]])


  called_tads = list(bed)

  curr_lev = 2

  while (curr_lev != (levels +1) ) {

  coords = cbind(match(called_tads[[curr_lev-1]]$start, as.numeric(colnames(cont_mat))), match(called_tads[[curr_lev-1]]$end-resolution, as.numeric(colnames(cont_mat))))

  #Get rows that are less than 10 in length and thus not seperable

  less_5 = which( (coords[,2]-coords[,1])<10  )

  if (length(less_5)>0) {

  #Pull them out

  pres_tads = called_tads[[curr_lev-1]][less_5,]

  coords = coords[-less_5, ]

  } else {
    pres_tads = c()
  }

  tads = apply(coords, 1, function(x) cont_mat[x[1]:x[2], x[1]:x[2]])

  #Remove sub-tads with too many zeros

  zeros = which(unlist(lapply(tads, function(x) nrow(x)-sum(rowSums(x)==0)))<10)

  if (length(zeros)>0) {
  pres_tads = rbind(pres_tads, called_tads[[curr_lev-1]][zeros,])
  tads[zeros] = NULL
  }


  sub_tads = lapply(tads, function(x) {
    windowedSpec(x, chr =chr, resolution = resolution, qual_filter = qual_filter, z_clust = FALSE, min_size = min_size)
    })

  called_tads[[curr_lev]] = bind_rows(sub_tads, pres_tads) %>% mutate(Level = curr_lev) %>% arrange(start)

  curr_lev = curr_lev+1

  }

  names(called_tads) = paste0("Level_", 1:levels)

  return(called_tads)
}
