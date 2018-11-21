#' Hierarchical Spectral Clustering of TADs
#'
#' Given a sparse 3 column, an n x n contact matrix or n x (n+3) contact matrix, hierarchicalSpec returns a list of TAD coordinates in bed form.
#' @param cont_mat Contact matrix in either sparse 3 column, n x n or n x (n+3) form where the first 3 columns are a bed file of coordinates. If a n x n matrix is used, the column names must correspond to the start point of the corresponding bin. Required.
#' @param chr The chromosome of the contact matrix being analyzed. Required.
#' @param levels The number of levels of the TAD hierarchy to be calculated, The default setting is 1.
#' @param qual_filter Option to turn on quality filtering which removes TADs with negative silhouette scores. The default setting is TRUE.
#' @param z_clust Option to filter sub-TADs based on the z-score of their eigenvector gaps. The default setting is FALSE.
#' @param eigenvalues The number of eigenvectors to be calculated. The default and suggested setting is 2.
#' @param min_size The minimum allowable TAD size measured in bins. Defaults to 5.
#' @param resolution The resolution of the contact matrix. If none selected the resolution is estimated by taking the difference in start points between the first and second bin. For n x (n+3) contact matrices this value is automatically calculated from the first 3 columns.
#' @keywords spectral clustering
#' @export
#' @details
#' hierarchicalSpec()

hierarchicalSpec = function(cont_mat, chr, levels = 1, qual_filter = TRUE, z_clust = FALSE, eigenvalues = 2, min_size =5, resolution = "auto") {
  source("windowedSpec.R")

  #Calculate the number of rows and columns of the contact matrix

  row_test = dim(cont_mat)[1]
  col_test = dim(cont_mat)[2]

  if (col_test == 3) {

    #Convert sparse matrix to n x n matrix

    cont_mat = HiCcompare::sparse2full(cont_mat)
    message("Converting to n x n matrix")
  } else if (col_test-row_test == 3) {

    message("Converting to n x n matrix")

    #Find the start coordinates based on the second column of the bed file portion of matrix

    start_coords = cont_mat[,2]

    #Calculate resolution based on given bin size in bed file

    resolution = as.numeric(cont_mat[1,3])-as.numeric(cont_mat[1,2])

    #Remove bed file portion

    cont_mat = cont_mat[,-c(1:3)]

    #Make column names correspond to bin start

    colnames(cont_mat) = start_coords

  } else if (col_test!=3 & row_test != col_test & col_test-row_test != 3) {

    #Throw error if matrix does not correspond to known matrix type

    stop("Contact matrix must be sparse or n x n or n x (n+3)!")
    break
  }

  #Automatically estimate resolution based on distance between first two bin start points

  if (resolution == "auto") {
    resolution = as.numeric(colnames(cont_mat)[2])-as.numeric(colnames(cont_mat)[1])
  } else {
    resolution = resolution
  }

  #Performed window spectral clustering

  bed = windowedSpec(cont_mat, chr = chr, resolution = resolution, z_clust = z_clust, eigenvalues = eigenvalues, min_size = min_size, qual_filter = qual_filter) %>% mutate(Level = 1)

  #Calculate the end point of TADs based on bin instead of genomic coordinate

  coords = cbind(match(bed$start, as.numeric(colnames(cont_mat))), match(bed$end-resolution, as.numeric(colnames(cont_mat))))

  #Create a list of tad start and end points

  tads = apply(coords, 1, function(x) cont_mat[x[1]:x[2], x[1]:x[2]])


  called_tads = list(bed)

  #Initialize TAD level

  curr_lev = 2

  while (curr_lev != (levels + 1) ) {

  #Get a list of TAD coordinates at the previous level

  coords = cbind(match(called_tads[[curr_lev-1]]$start, as.numeric(colnames(cont_mat))), match(called_tads[[curr_lev-1]]$end-resolution, as.numeric(colnames(cont_mat))))

  #Get tads that are less than the twice the minmium length and thus not seperable

  less_5 = which( (coords[,2]-coords[,1])<min_size*2  )

  if (length(less_5)>0) {

  #Remove TADs which cannot be seperate

  pres_tads = called_tads[[curr_lev-1]][less_5,]

  coords = coords[-less_5, ]

  } else {
    pres_tads = c()
  }

  tads = apply(coords, 1, function(x) cont_mat[x[1]:x[2], x[1]:x[2]])

  #Remove sub-tads with too many zeros

  zeros = which(unlist(lapply(tads, function(x) nrow(x)-sum(rowSums(x)==0)))<min_size*2)

  if (length(zeros)>0) {
  pres_tads = rbind(pres_tads, called_tads[[curr_lev-1]][zeros,])
  tads[zeros] = NULL
  }

  #Calculate sub-TADs for each seperable TAD

  sub_tads = lapply(tads, function(x) {
    windowedSpec(x, chr =chr, resolution = resolution, qual_filter = qual_filter, z_clust = FALSE, min_size = min_size)
    })

  #Create a bed file for sub-TADs

  called_tads[[curr_lev]] = bind_rows(sub_tads, pres_tads) %>% mutate(Level = curr_lev) %>% arrange(start)

  curr_lev = curr_lev+1

  }

  #Assign names based on levels

  names(called_tads) = paste0("Level_", 1:levels)

  return(called_tads)
}
