#' Hierarchical Spectral Clustering of TADs
#'
#' @import dplyr
#' @param cont_mat Contact matrix in either sparse 3 column, n x n or n x (n+3)
#' form where the first three columns are a bed file of coordinates.
#' If an x n matrix is used, the column names must correspond to the start
#' point of the corresponding bin. Required.
#' @param chr The chromosome of the contact matrix being analyzed. Required.
#' @param levels The number of levels of the TAD hierarchy to be calculated.
#' The default setting is 1.
#' @param qual_filter Option to turn on quality filtering which removes TADs
#' with negative silhouette scores (poorly organized TADs). Default is TRUE.
#' @param z_clust Option to filter sub-TADs based on the z-score of
#' their eigenvector gaps. Default is FALSE.
#' @param eigenvalues The number of eigenvectors to be calculated.
#' The default and suggested setting is 2.
#' @param min_size The minimum allowable TAD size measured in bins. Default is 5.
#' @param resolution The resolution of the contact matrix. If none selected,
#' the resolution is estimated by taking the most common distance between bins.
#' For n x (n+3) contact matrices, this value is automatically calculated
#' from the first three columns.
#' @param gap_threshold Corresponds to the percentage of zeros allowed before
#' a column/row is removed from the analysis. 1=100\%, .7 = 70\%, etc. Default is 1.
#' @return A list where each entry is a bed file corresponding to the level of the hierarchy.
#' @export
#' @details Given a sparse 3 column, an n x n contact matrix,
#' or n x (n+3) contact matrix, SpectralTAD returns a list of TAD coordinates
#' in BED format. SpectralTAD works by using a sliding window that moves along
#' the diagonal of the contact matrix. By default, we use the biologically
#' relevant maximum TAD size of 2Mb and minimum size of 5 bins to determine
#' the size of this window. Within each window, we calculate a Laplacian matrix
#' and determine the location of TAD boundaries based on gaps between
#' eigenvectors calculated from this matrix. The number of TADs in a given
#' window is calculated by finding the number that maximizes the silhouette score.
#' A hierarchy of TADs is created by iteratively applying the function to
#' sub-TADs. The number of levels in each hierarchy is determined by the user.
#' @examples
#' #Read in data
#' data("rao_chr20_25_rep")
#' #Find TADs
#' spec_table <- SpectralTAD(rao_chr20_25_rep, chr= 'chr20')


SpectralTAD = function(cont_mat, chr, levels = 1, qual_filter = TRUE,
                       z_clust = FALSE, eigenvalues = 2, min_size = 5, resolution = "auto", gap_threshold = 1) {

  #Calculate the number of rows and columns of the contact matrix

  if (missing("chr")) {
    stop("Must specify chromosome")
  }

  row_test = dim(cont_mat)[1]
  col_test = dim(cont_mat)[2]

  if (row_test == col_test) {
    if (all(is.finite(cont_mat)) == FALSE) {
      stop("Contact matrix must only contain real numbers")
    }
  }

  if (col_test == 3) {

    #Convert sparse matrix to n x n matrix

    message("Converting to n x n matrix")

    cont_mat = HiCcompare::sparse2full(cont_mat)

    if (all(is.finite(cont_mat)) == FALSE) {
      stop("Contact matrix must only contain real numbers")
    }

    if (resolution == "auto") {
      message("Estimating resolution")
      resolution = as.numeric(names(table(as.numeric(colnames(cont_mat))-lag(as.numeric(colnames(cont_mat)))))[1])
    }

  } else if (col_test-row_test == 3) {

    message("Converting to n x n matrix")

    #Find the start coordinates based on the second column of the bed file portion of matrix

    start_coords = cont_mat[,2]

    #Calculate resolution based on given bin size in bed file

    resolution = as.numeric(cont_mat[1,3])-as.numeric(cont_mat[1,2])

    #Remove bed file portion

    cont_mat = as.matrix(cont_mat[,-c(1:3)])

    if (all(is.finite(cont_mat)) == FALSE) {
      stop("Contact matrix must only contain real numbers")
    }

    #Make column names correspond to bin start

    colnames(cont_mat) = start_coords

  } else if (col_test!=3 & (row_test != col_test) & (col_test-row_test != 3)) {

    #Throw error if matrix does not correspond to known matrix type

    stop("Contact matrix must be sparse or n x n or n x (n+3)!")
    break
  } else if ( (resolution == "auto") & (col_test-row_test == 0) ) {
      message("Estimating resolution")

      #Estimating resolution based on most common distance between loci

      resolution = as.numeric(names(table(as.numeric(colnames(cont_mat))-lag(as.numeric(colnames(cont_mat)))))[1])
  }

  #Performed window spectral clustering

  bed = .windowedSpec(cont_mat, chr = chr, resolution = resolution, z_clust = z_clust, eigenvalues = eigenvalues, min_size = min_size, qual_filter = qual_filter, gap_threshold = gap_threshold) %>% mutate(Level = 1)

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

  #Account for the situation where there is only 1 potential sub-tad

  if (is.null(nrow(coords))) {
    coords = t(as.matrix(coords))
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
    .windowedSpec(x, chr =chr, resolution = resolution, qual_filter = qual_filter, z_clust = TRUE, min_size = min_size)
    })

  #Create a bed file for sub-TADs

  called_tads[[curr_lev]] = bind_rows(sub_tads, pres_tads) %>% mutate(Level = curr_lev) %>% arrange(start)

  curr_lev = curr_lev+1

  }

  #Assign names based on levels

  names(called_tads) = paste0("Level_", 1:levels)

  return(called_tads)
}




#Function to perform the actual sliding window Spectral clustering
#Used within SpectralTAD

.windowedSpec = function(cont_mat, resolution, chr,
                        gap_filter = TRUE,z_clust = FALSE,  qual_filter = TRUE, eigenvalues = 2, min_size = 5, gap_threshold = 1) {

  #Set window sized based on biologically maximum TAD size of 2000000

  window_size = ceiling(2000000/resolution)

  #Find all regions which aren't completely zero and remove those that are


  #Get end point of the first window

  Group_over = dplyr::bind_rows()

  #Initialize first window

  start = 1
  end = window_size

  #Set parameter for determining end of loop

  end_loop = 0

  #Test if start+window is larger than the contact matrix and correct end point

  if (end+window_size>nrow(cont_mat)) {
    end = nrow(cont_mat)
  }

  #Begin sliding window clustering

  while (end_loop == 0) {

    #Subset matrix based on window size

    sub_filt = cont_mat[start:end, start:end]

    zero_thresh = round(nrow(sub_filt)*(gap_threshold))
    non_gaps_within = which((colSums(sub_filt == 0))<zero_thresh)

    sub_filt = sub_filt[non_gaps_within, non_gaps_within]

    if (nrow(sub_filt) < min_size*2) {
      start = end
      end = start+window_size
      next
    }

    # sub_gaps = colSums(sub_filt)>0
    # sub_filt = sub_filt[sub_gaps, sub_gaps]

    #Calculate distance matrix for silhouette score

    dist_sub = 1/(1+sub_filt)

    #Get degree matrix

    dr = rowSums(abs(sub_filt))

    #Creating the normalized laplacian

    Dinvsqrt = diag((1/sqrt(dr)))

    P_Part1 = Matrix::crossprod(as.matrix(sub_filt), Dinvsqrt)
    sub_mat = Matrix::crossprod(Dinvsqrt, P_Part1)

    colnames(sub_mat) = colnames(cont_mat)[non_gaps_within]

    sub_mat[is.nan(sub_mat)] = 0

    #Get first k eigenvectors

    Eigen = PRIMME::eigs_sym(sub_mat, NEig = eigenvalues)

    eig_vals = Eigen$values
    eig_vecs = Eigen$vectors

    #Get order of eigenvalues from largest to smallest

    large_small = order(-eig_vals)

    eig_vals = eig_vals[large_small]
    eig_vecs = eig_vecs[,large_small]

    index = 1
    Group_mem = list()

    #Calculate the range of possible clusters

    clusters = 1:ceiling( (end-start+1)/min_size)

    #Normalize the eigenvectors from 0-1

    norm_ones = sqrt(dim(sub_mat)[2])

    for (i in 1:dim(eig_vecs)[2]) {
      eig_vecs[,i] = (eig_vecs[,i]/sqrt(sum(eig_vecs[,i]^2)))  * norm_ones
      if (eig_vecs[1,i] !=0) {
        eig_vecs[,i] = -1*eig_vecs[,i] * sign(eig_vecs[1,i])
      }
    }

    n = dim(eig_vecs)[1]
    k = dim(eig_vecs)[2]

    #Project eigenvectors onto a unit circle

    eig_vecs = crossprod(diag(diag(tcrossprod(eig_vecs))^(-1/2)), eig_vecs)

    #Get distance between points on circle

    point_dist = sqrt(rowSums( (eig_vecs-rbind(NA,eig_vecs[-nrow(eig_vecs),]))^2  ))

    #Use z-score to select significant gaps

    if (z_clust) {

      #Get statisticaly significant boundaries

      sig_bounds = which(scale(point_dist[-length(point_dist)])>2)

      #Remove boundaries within the minimum size

      sig_bounds = subset(sig_bounds, sig_bounds>min_size)

      #2*min_size is to offset and remove the second occurence

      dist_bounds = which(c(min_size*2,diff(sig_bounds))<min_size)

      #Remove bounds within the mininum size if they exist

      if (length(dist_bounds) > 0) {
        sig_bounds = sig_bounds[-dist_bounds]
      }

      #Create TADs using significant boundaries

      TAD_start = c(1, sig_bounds+1)

      TAD_end = c(sig_bounds, nrow(sub_filt))

      widths = (TAD_end-TAD_start)+1

      memberships = unlist(lapply(1:length(TAD_start), function(x) rep(x,widths[x])))

      #Create groups

      if (length(sig_bounds) == 0) {

        #Create empty set if non-significant

        end_group = dplyr::bind_rows()
      } else {

        sig_bounds = which(scale(point_dist[-length(point_dist)])>2)

        #Remove boundaries within the minimum size

        sig_bounds = subset(sig_bounds, sig_bounds>min_size)

        #2*min_size is to offset and remove the second occurence

        dist_bounds = which(c(min_size*2,diff(sig_bounds))<min_size)

        #Assign IDs based on coordinate and groups based on significant boundaries

        end_group = data.frame(ID = as.numeric(colnames(sub_filt)), Group = memberships)

        #Compile into bed file

        end_group = end_group %>% dplyr::mutate(group_place = Group) %>% dplyr::group_by(group_place) %>% dplyr::mutate(Group = last(ID)) %>% dplyr::ungroup() %>% dplyr::select(ID, Group)

      }


    } else {


      #Find largest gaps

      gap_order = order(-point_dist)

      #Remove boundaries occuring before minimum size at the very beginning of window

      #gap_order = gap_order[-which(gap_order<min_size)]

      #Initialize silhouette score

      sil_score = c()


      for (cluster in clusters) {

        #Loop through first k gaps and remove repeating boundaries

        #Set intial cutpoints to the number of clusters

        k = 1
        partition_found = 0
        first_run = TRUE
        cutpoints = c()

        #Loop through cluster numbers by iteratively adding new candidate boundaries and testing

        while(partition_found == 0) {

          #Get candidate gaps

          new_gap = gap_order[k]

          cutpoints = c(cutpoints, new_gap)

          #Identify gaps which are closer together than the minimum TAD size

          diff_points = which( abs(new_gap-cutpoints[-length(cutpoints)]) <= min_size)

          #If a point exists that is too close to another, remove it

          if (length(diff_points)>0) {
            cutpoints = cutpoints[-length(cutpoints)]
          }

          #If not these are final clusters

          if (length(cutpoints) == cluster) {
            partition_found = 1
          } else {
            k = k+1
          }
        }

        #If the new candidate cluster is an NA value, ignore

        if (any(is.na(cutpoints))) {
          next
        }

        #Order

        cutpoints = cutpoints[order(cutpoints)]

        #Combine cutpoints with start and end of window

        cutpoints = c(1, cutpoints, length(non_gaps_within)+1)

        #Find size of each cluster (TAD)

        group_size = diff(cutpoints)

        #Assign locations of the window memberships based on cutpoints

        memberships = c()
        for (i in 1:length(group_size)) {
          memberships = c(memberships, rep(i,times = group_size[i]))
        }

        #Get silhouette score for current number of clusters (TADs)

        sil = summary(cluster::silhouette(memberships,dist_sub))

        #Save silhouette scores for each configuration in vector

        sil_score = c(sil_score, sil$si.summary[4])

        #Save memberships in list

        Group_mem[[cluster]] = memberships

      }



      #Pull out the cutpoints which maximize silhouette score

      end_group = Group_mem[[which(diff(sil_score)<0)[1]]]

      #Put coordinates and group IDs into data frame

      if (length(end_group) == 0) {
        end_group = dplyr::bind_rows()
      } else {

      end_group = data.frame(ID = as.numeric(colnames(sub_filt)), Group = end_group)

      #Convert IDs to coordinates of endpoint to avoid overlaps

      end_group = end_group %>%dplyr::mutate(group_place = Group) %>%dplyr::group_by(group_place) %>%dplyr::mutate(Group = max(ID)) %>% ungroup() %>% dplyr::select(ID, Group)
      }
    }

    #End while loop if window reaches end of contact matrix

    if (end == nrow(cont_mat)) {
      Group_over = dplyr::bind_rows(Group_over, end_group)
      end_loop = 1
    } else {

      #Remove the last group (To account for overlap across windows) and set new start to start of group

      end_IDs = which(end_group$Group == last(end_group$Group))

      start = end-length(end_IDs)+1

      #Account for cases when final TAD can't be removed

      if (length(start) == 0 ) {
        start = end
      }

      #Set new window end

      end = start+window_size

      #Remove final group to avoid repeating

      end_group = end_group[-end_IDs, ]

      #Combine TAD coordinates into single bed file

      Group_over = dplyr::bind_rows(Group_over, end_group)

      #Set end point to end of contact matrix if window is larger than end of matrix

      if ( (end + (2000000/resolution)) > nrow(cont_mat) ) {
        end = nrow(cont_mat)

      }
    }
  }


  #Organize final results based on options selected

  if (z_clust) {

    if (nrow(Group_over) > 0) {
      bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>%dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%
        dplyr::filter((end-start)/resolution >= min_size) %>%dplyr::arrange(start)
    } else {
      bed = Group_over
    }
  } else {

    if (qual_filter) {

      #Calculate an overall distance matrix for calculating silhouette score for filtering

      #Get range of values in the contact matrix

      fin_range = match(Group_over$ID,colnames(cont_mat))

      over_dist_mat = 1/(1+cont_mat[fin_range, fin_range])


      #Calculate group-wise silhouette

      sil = cluster::silhouette(Group_over$Group, over_dist_mat)

      ave_sil = summary(sil)$clus.avg.widths

      #Subset results based on silhouette score depending on qual_filter option

      bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>% dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%
        dplyr::mutate(Sil_Score = ave_sil) %>% dplyr::filter( ((end-start)/resolution >= min_size) & Sil_Score > .15)  %>%dplyr::arrange(start)
    } else {
      bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>% dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%dplyr::filter((end-start)/resolution >= min_size) %>% dplyr::arrange(start)
    }
  }

  return(bed)
}
