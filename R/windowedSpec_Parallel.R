par_spec = function(cont_list, chr_over, labels = NA, levels = 1,  qual_filter = FALSE, z_clust = FALSE, eigenvalues = 2, min_size =5, cores = "auto") {

  require(BiocParallel)
  source("hierarchicalSpec.R")
  
  if (cores == "auto") {
  # Check how many cores you have
  numCores <- parallel::detectCores()
  } else {
    numCores = cores
  }
  
  # Set the number of cores at least one less than the total number
  if(Sys.info()['sysname'] == "Windows") {
    # Windows settings
    register(SnowParam(workers = numCores - 1), default = TRUE)
  } else {
    # Unix settings
    register(MulticoreParam(workers = numCores - 1), default = TRUE)
  }
  
bed = bplapply(1:length(cont_list), function(x, spec_fun, cont_list,eigenvalues, z_clust, qual_filter, levels, min_size, chr_over) spec_fun(cont_mat = cont_list[[x]],  chr_over[[x]], levels = levels, qual_filter = qual_filter, z_clust = z_clust, eigenvalues = eigenvalues, min_size = min_size), spec_fun = hierarchicalSpec, cont_list = cont_list, eigenvalues = eigenvalues, min_size = min_size, levels = levels,z_clust = z_clust, qual_filter = qual_filter, chr_over = chr_over)
#bed = bind_rows(bed)

if (!is.na(labels)) {
  names(bed) = labels
}

return(bed)
}




