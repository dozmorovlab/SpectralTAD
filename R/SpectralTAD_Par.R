#' Parallelized Hierarchical Spectral Clustering of TADs
#'
#' @importFrom BiocParallel bplapply MulticoreParam register SnowParam
#' @param cont_list List of contact matrices where each is in either
#' sparse 3 column, n x n or n x (n+3) form, where the first 3 columns
#' are chromosome, start and end coordinates of the regions.
#' If an x n matrix is used, the column names must correspond to
#' the start point of the corresponding bin. Required.
#' @param chr Vector of chromosomes in the same order as their
#' corresponding contact matrices. Must be same length as cont_list. Required.
#' @param levels The number of levels of the TAD hierarchy to be calculated.
#' The default setting is 1.
#' @param qual_filter Option to turn on quality filtering which removes TADs
#' with negative silhouette scores (poorly organized TADs). Default is FALSE.
#' @param z_clust Option to filter sub-TADs based on the z-score of
#' their eigenvector gaps. Default is TRUE.
#' @param eigenvalues The number of eigenvectors to be calculated.
#' The default and suggested setting is 2.
#' @param min_size The minimum allowable TAD size measured in bins. Default is 5.
#' @param window_size The size of the sliding window for calculating TADs.
#' Smaller window sizes correspond to less noise from long-range contacts
#' but limit the possible size of TADs
#' @param resolution The resolution of the contact matrix. If none selected,
#' the resolution is estimated by taking the most common distance between bins.
#' For n x (n+3) contact matrices, this value is automatically calculated
#' from the first 3 columns.
#' @param gap_threshold Corresponds to the percentage of zeros allowed before
#' a column/row is removed from analysis. 1=100\%, .7 = 70\%, etc. Default is 1.
#' @param grange Parameter to determine whether the result should be a 
#' GRangeList object. Defaults to FALSE
#' @param cores Number of cores to use. Defaults to total available
#' cores minus one.
#' @param labels Vector of labels used to name each contact matrix. Must be
#' same length as cont_list. Default is NULL.
#' @export
#' @return List of lists where each entry is a list of data frames or GRanges
#'  in BED format corresponding to TADs seperated by hierarchies
#' @details This is the parallelized version of the SpectralTAD() function.
#' Given a sparse 3 column, an n x n contact matrix,
#' or n x (n+3) contact matrix, SpectralTAD returns a list of TAD coordinates
#' in BED format. SpectralTAD works by using a sliding window that moves along
#' the diagonal of the contact matrix. By default we use the biologically
#' relevant maximum TAD size of 2Mb and minimum size of 5 bins to determine
#' the size of this window. Within each window, we calculate a Laplacian matrix
#' and determine the location of TAD boundaries based on gaps between
#' eigenvectors calculated from this matrix. The number of TADs in a given
#' window is calculated by finding the number that maximize the silhouette score.
#' A hierarchy of TADs is created by iteratively applying the function to
#' sub-TADs. The number of levels in each hierarchy is determined by the user.
#' @examples
#' #Read in data
#' data("rao_chr20_25_rep")
#' #Make a list of matrices
#' mat_list = list(rao_chr20_25_rep, rao_chr20_25_rep)
#' #Make a vector of chromosomes
#' chr = c("chr20", "chr20")
#' #Make a vector of labels
#' labels = c("run1", "run2")
#' spec_table <- SpectralTAD_Par(mat_list, chr= chr, labels = labels)

SpectralTAD_Par = function(cont_list, chr, levels = 1,
                           qual_filter = FALSE, 
                           z_clust = FALSE, 
                           eigenvalues = 2,
                           min_size =5,
                           window_size = 25,
                           resolution = "auto", grange = FALSE, 
                           gap_threshold = 1, cores = "auto", 
                           labels = NULL) {

  if ( (is.data.frame(cont_list)) | (is.matrix(cont_list))) {
    stop("Input must be list of contact matrices")
  }
  
  if (length(chr) != length(cont_list)) {
    stop("Length of chromosome vector must match number of contact matrices")
  }
  

  if (cores == "auto") {
  # Check how many cores you have
  numCores <- parallel::detectCores()
  } else {
    numCores = cores
  }

  # Set the number of cores at least one less than the total number
  if(Sys.info()['sysname'] == "Windows") {
    # Windows settings
    BiocParallel::register(BiocParallel::SnowParam(workers = numCores - 1), default = TRUE)
  } else {
    # Unix settings
    BiocParallel::register(BiocParallel::MulticoreParam(workers = numCores - 1), default = TRUE)
  }

  #Run SpectralTAD simultaneously on each contact matrix

  bed = BiocParallel::bplapply(seq_len(length(cont_list)), function(x, spec_fun, cont_list,
                                                             eigenvalues, z_clust, qual_filter, levels, min_size,
                                                             chr, gap_threshold, grange) spec_fun(cont_mat = cont_list[[x]],
                                                             chr[[x]], levels = levels, qual_filter = qual_filter,
                                                             z_clust = z_clust, eigenvalues = eigenvalues, min_size = min_size,
                                                             window_size = window_size),
                                                             spec_fun = SpectralTAD, cont_list = cont_list, eigenvalues = eigenvalues,
                                                             min_size = min_size, levels = levels,z_clust = z_clust, 
                                                             window_size = window_size, qual_filter = qual_filter,
                                                             chr = chr, gap_threshold = gap_threshold, grange = grange)

  #Assign labels to identify each contact matrix

  if (!is.null(labels)) {
    if (length(labels) != length(cont_list)) {
      warning("Labels ignored: Labels not same length as list of contact matrices")
      labels = NULL
    }
    names(bed) = labels
  }


return(bed)
}





