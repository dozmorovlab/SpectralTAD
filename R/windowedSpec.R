windowedSpec = function(cont_mat, resolution, chr,
                        gap_filter = TRUE,z_clust = FALSE,  qual_filter = TRUE, eigenvalues = 2, min_size = 5) {

  #Set window sized based on biologically maximum TAD size of 2000000

  window_size = ceiling(2000000/resolution)

  #Find all regions which aren't completely zero and remove those that are

  non_gaps = which(colSums(cont_mat) !=0)

  cont_mat = cont_mat[non_gaps,non_gaps]

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

    #Calculate distance matrix for silhouette score

    dist_sub = 1/(1+sub_filt)

    #Get degree matrix

    dr = rowSums(abs(sub_filt))

    #Creating the normalized laplacian

    Dinvsqrt = diag((1/sqrt(dr+2e-16)))

    P_Part1 = Matrix::crossprod(as.matrix(sub_filt), Dinvsqrt)
    sub_mat = Matrix::crossprod(Dinvsqrt, P_Part1)

    colnames(sub_mat) = colnames(cont_mat)[start:end]

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

    vm = matrix(kronecker(rep(1,k), as.matrix(sqrt(rowSums(eig_vecs^2)))),n,k)
    eig_vecs = eig_vecs/vm

    #Get distance between points on circle

    point_dist = sqrt(rowSums( (eig_vecs-rbind(NA,eig_vecs[-nrow(eig_vecs),]))^2  ))

    #Use z-score to select significant gaps

    if (z_clust) {

      #Get statisticaly significant boundaries

      sig_bounds = which(scale(point_dist[-length(point_dist)])>2)

      #Remove boundaries within the minimum size

      sig_bounds = subset(sig_bounds, sig_bounds>min_Size)

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

        #Assign IDs based on coordinate and groups based on significant boundaries

        end_group = data.frame(ID = as.numeric(colnames(sub_filt)), Group = memberships)

        #Compile into bed file

        end_group = end_group %>% dplyr::mutate(group_place = Group) %>% dplyr::group_by(group_place) %>% dplyr::mutate(Group = last(ID)) %>% dplyr::ungroup() %>% dplyr::select(ID, Group)

      }


    } else {

    #Find largest gaps

    gap_order = order(-point_dist)

    #Remove boundaries occuring before minimum size at the very beginning of window

    gap_order = gap_order[-which(gap_order<min_size)]

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

      cutpoints = c(1, cutpoints, (end-start)+2)

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

    end_group = Group_mem[[which.max(sil_score)]]

    #Put coordinates and group IDs into data frame

    end_group = data.frame(ID = as.numeric(colnames(sub_filt)), Group = end_group)

    #Convert IDs to coordinates of endpoint to avoid overlaps

    end_group = end_group %>%dplyr::mutate(group_place = Group) %>%dplyr::group_by(group_place) %>%dplyr::mutate(Group = max(ID)) %>% ungroup() %>% dplyr::select(ID, Group)

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

    over_dist_mat = 1/(1+cont_mat)

    sil = cluster::silhouette(Group_over$Group, over_dist_mat)

    ave_sil = summary(sil)$clus.avg.widths

    #Subset results based on silhouette score depending on qual_filter option

    bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>% dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%
     dplyr::mutate(Sil_Score = ave_sil) %>%dplyr::filter( ((end-start)/resolution >= min_size) & Sil_Score > .15)  %>%dplyr::arrange(start)
  } else {
    bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>% dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%dplyr::filter((end-start)/resolution >= min_size) %>% dplyr::arrange(start)
  }
  }

  return(bed)
}
