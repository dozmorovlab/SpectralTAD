windowedSpec = function(cont_mat, resolution, chr,
                        gap_filter = TRUE,z_clust = FALSE,  qual_filter = TRUE, eigenvalues = 2, min_size = 5) {

  require(PRIMME)
  require(Matrix)
  require(dplyr)
  require(cluster)
  source("GroupsFix2.R")

  window_size = ceiling(2000000/resolution)

  non_gaps = which(colSums(cont_mat) !=0)

  cont_mat = cont_mat[non_gaps,non_gaps]

  max_end = ceiling(2000000/resolution)
  max_size = window_size/ceiling(200000/resolution)
  min_size = ceiling(200000/resolution)
  Group_over = bind_rows()

  start = 1
  end = max_end

  end_loop = 0

  if (end+window_size>nrow(cont_mat)) {
    end = nrow(cont_mat)
  }

    while (end_loop == 0) {
    sub_filt = cont_mat[start:end, start:end]

    dist_sub = 1/(1+sub_filt)

    dr = rowSums(abs(sub_filt))

    #Creating the normalized laplacian

    Dinvsqrt = diag((1/sqrt(dr+2e-16)))

    P_Part1 = crossprod(as.matrix(sub_filt), Dinvsqrt)
    sub_mat = crossprod(Dinvsqrt, P_Part1)

    #sub_mat = crossprod(diag(dr^-(1/2)), as.matrix(sub_filt))

    colnames(sub_mat) = colnames(cont_mat)[start:end]

    #Find gaps at 2mb and remove


    #Remove gaps if true

    # if (length(sub_gaps) != 0) {
    #   sub_mat = sub_mat[!sub_gaps, !sub_gaps]
    # }


    #Get first two eigenvectors

    Eigen = eigs_sym(sub_mat, NEig = eigenvalues)

    eig_vals = Eigen$values
    eig_vecs = Eigen$vectors

    #Get order of eigenvalues from largest to smallest

    large_small = order(-eig_vals)

    eig_vals = eig_vals[large_small]
    eig_vecs = eig_vecs[,large_small]

    index = 1
    Group_mem = list()

    #When wide search is on we have an extra step that does a preliminary search for the best starting point

    clusters = 1:ceiling( (end-start+1)/min_size)

    #Normalize the eigenvectors

    norm_ones = sqrt(dim(sub_mat)[2])

    for (i in 1:dim(eig_vecs)[2]) {
      eig_vecs[,i] = (eig_vecs[,i]/sqrt(sum(eig_vecs[,i]^2)))  * norm_ones
      if (eig_vecs[1,i] !=0) {
        eig_vecs[,i] = -1*eig_vecs[,i] * sign(eig_vecs[1,i])
      }
    }


    eps = 2.2204e-16

    n = dim(eig_vecs)[1]
    k = dim(eig_vecs)[2]

    #Project eigenvectors onto a unit circle

    vm = matrix(kronecker(rep(1,k), as.matrix(sqrt(rowSums(eig_vecs^2)))),n,k)
    eig_vecs = eig_vecs/vm

    #Get distance between points on circle


    point_dist = sqrt(rowSums( (eig_vecs-rbind(NA,eig_vecs[-nrow(eig_vecs),]))^2  ))

    #Only consider points with z-scores higher than 2

    # if (strict ) {
    #
    #   if (sum(scale(point_dist)>2, na.rm = TRUE) == 0 ) {
    #   Group_over = bind_rows()
    #   end_loop = 1
    #   break
    #   }
    #
    # }

    if (z_clust) {

      #Get statisticaly significant boundaries

      sig_bounds = which(scale(point_dist[-length(point_dist)])>2)

      #Remove boundaries before 5 or within 5 of others

      sig_bounds = subset(sig_bounds, sig_bounds>5)

      #10 is to offset and remove the second occurence

      dist_bounds = which(c(10,diff(sig_bounds))<5)


      if (length(dist_bounds) > 0) {
      sig_bounds = sig_bounds[-dist_bounds]
      }

      #Create TADs

      TAD_start = c(1, sig_bounds+1)

      TAD_end = c(sig_bounds, nrow(sub_filt))

      widths = (TAD_end-TAD_start)+1

      memberships = unlist(lapply(1:length(TAD_start), function(x) rep(x,widths[x])))

      #Create groups

      if (length(sig_bounds) == 0) {

        end_group = bind_rows()
      } else {
      end_group = data.frame(ID = as.numeric(colnames(sub_filt)), Group = memberships)


      end_group = end_group %>% dplyr::mutate(group_place = Group) %>% dplyr::group_by(group_place) %>% dplyr::mutate(Group = last(ID)) %>% dplyr::ungroup() %>% dplyr::select(ID, Group)

      }


    } else {

    #Find largest gaps

    gap_order = order(-point_dist)

    #Remove boundaries occuring before 200kb

    gap_order = gap_order[-which(gap_order<5)]

    sil_score = c()


    for (cluster in clusters) {

      #Loop through first k gaps and remove repeating boundaries

      #Set intial cutpoints to the number of clusters

      k = 1
      partition_found = 0
      first_run = TRUE
      cutpoints = c()

      while(partition_found == 0) {

        #Get candidate gaps

        new_gap = gap_order[k]

        cutpoints = c(cutpoints, new_gap)

        diff_points = which( abs(new_gap-cutpoints[-length(cutpoints)]) <= min_size)

        #Remove the point

        if (length(diff_points)>0) {
          cutpoints = cutpoints[-length(cutpoints)]
        }

        if (length(cutpoints) == cluster) {
          partition_found = 1
        } else {
          k = k+1
        }
      }

      if (any(is.na(cutpoints))) {
        next
      }

      #Order


      cutpoints = cutpoints[order(cutpoints)]
      cutpoints = c(1, cutpoints, (end-start)+2)

      #Add 1 to final cutpoint to account for offset

      #Find size of groups

      group_size = diff(cutpoints)

      #Add start and end

      memberships = c()
      for (i in 1:length(group_size)) {
        memberships = c(memberships, rep(i,times = group_size[i]))
      }

      #Get silhouette score
      sil = summary(silhouette(memberships,dist_sub))

      sil_score = c(sil_score, sil$si.summary[4])

      Group_mem[[cluster]] = memberships

    }

    #Get all IDs in the last group



    end_group = Group_mem[[which.max(sil_score)]]

    end_group = data.frame(ID = as.numeric(colnames(sub_filt)), Group = end_group)


    end_group = end_group %>%dplyr::mutate(group_place = Group) %>%dplyr::group_by(group_place) %>%dplyr::mutate(Group = max(ID)) %>% ungroup() %>% dplyr::select(ID, Group)

    }

    if (end == nrow(cont_mat)) {
      Group_over = bind_rows(Group_over, end_group)
      end_loop = 1
    } else {

      #Remove the last group and set new start to start of group

      end_IDs = which(end_group$Group == last(end_group$Group))


      start = end-length(end_IDs)+1

      if (length(start) == 0 ) {
        start = end
      }

      end = start+max_end

      end_group = end_group[-end_IDs, ]
      Group_over = bind_rows(Group_over, end_group)

      if ( (end + (2000000/resolution)) > nrow(cont_mat) ) {
        end = nrow(cont_mat)

      }
    }
    }


  if (z_clust) {

    if (nrow(Group_over) > 0) {
    bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>%dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%
   dplyr::filter((end-start)/resolution >= min_size) %>%dplyr::arrange(start)
    } else {
      bed = Group_over
    }
  } else {

  if (qual_filter) {

    over_dist_mat = 1/(1+cont_mat)

    sil = silhouette(Group_over$Group, over_dist_mat)

    ave_sil = summary(sil)$clus.avg.widths

    bed = Group_over %>%dplyr::group_by(Group) %>%dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>%dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%
     dplyr::mutate(Sil_Score = ave_sil) %>%dplyr::filter( ((end-start)/resolution >= min_size) & Sil_Score > .15)  %>%dplyr::arrange(start)
  } else {
    bed = Group_over %>%dplyr::group_by(Group) %>%dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>%dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%dplyr::filter((end-start)/resolution >= min_size) %>%dplyr::arrange(start)
  }
  }

  return(bed)
}
