library(glue)
library(Matrix)

############################################################
# Specificity scores
############################################################
get_average_expression = function(scaled_matrix, grouping_vector) {
  if (ncol(scaled_matrix) == 0 || nrow(scaled_matrix) == 0) { stop('CDS has either no cells or no features, must have at least one cell and feature.')}
  if (! length(grouping_vector) == ncol(scaled_matrix)) { stop('grouping vector is not the same length as the number of cells.') }

  groups = names(table(grouping_vector))[table(grouping_vector) > 30]

  # Now get average profile for all groups
  average_profiles = lapply(groups, function(group) {
    cells = colnames(scaled_matrix)[grouping_vector == group]
    mat = scaled_matrix[, cells]

    # Calculate an average for each group on size factor normalized data
    return(Matrix::rowMeans(mat, na.rm=TRUE))
  })

  names(average_profiles) = groups
  average_profiles = as.matrix(dplyr::bind_cols(average_profiles))
  rownames(average_profiles) = rownames(scaled_matrix)
  return(average_profiles)
}

get_median_depth = function(mat, grouping_vector) {
  if (ncol(mat) == 0 || nrow(mat) == 0) { stop('CDS has either no cells or no features, must have at least one cell and feature.')}
  if (! length(grouping_vector) == ncol(mat)) { stop('grouping vector is not the same length as the number of cells.') }

  groups = names(table(grouping_vector))[table(grouping_vector) > 30]

  # Now get average profile for all groups
  cell_totals = Matrix::colSums(mat)

  median_total_per_cluster = lapply(groups, function(group) {
    group_totals = cell_totals[grouping_vector == group]
    return(median(group_totals))
  })

  median_total_per_cluster = unlist(median_total_per_cluster)
  names(median_total_per_cluster) = groups
  return(median_total_per_cluster)
}

make_probability_vector = function(p){
  phat = p / sum(p)
  phat[is.na(phat)] = 0
  return(phat)
}

shannon_entropy = function(p) {
  if (min(p) < 0 || sum(p) <= 0)
    return(Inf)
  p_norm = p[p > 0] / sum(p)
  return(-sum(log2(p_norm) * p_norm))
}

js_distance = function(p, q){
  js_div = shannon_entropy((p + q) / 2) - (shannon_entropy(p) + 
                                           shannon_entropy(q)) * 0.5
  js_div[is.infinite(js_div)] = 1
  js_div[js_div < 0] = 0
  js_dist = sqrt(js_div)
  return(js_dist)
}

calculate_specificity_helper = function(average_profile_matrix, groups_to_test=NULL) {
  total_groups = ncol(average_profile_matrix)
  group_names = colnames(average_profile_matrix)

  if (! is.null(groups_to_test)) {
    column_indices_to_test = seq(1, total_groups)[group_names %in% groups_to_test]
  } else {
    column_indices_to_test = seq(1, total_groups)
  }

  cluster_progress_count = 0
  total_groups_to_test = length(column_indices_to_test)

  marker_gene_specificities = lapply(column_indices_to_test, function(cell_type_index){
    cluster_progress_count <<- cluster_progress_count + 1

    message(glue('Calculating specificity scores for group: {group_names[cell_type_index]} ({cluster_progress_count} of {total_groups_to_test})...'))
    
    # will compare to perfect specificity where only seen in one cluster
    perfect_specificity <- rep(0.0, total_groups)
    perfect_specificity[cell_type_index] <- 1.0

    # Set up a progress bar
    pb = txtProgressBar(min = 0, max = nrow(average_profile_matrix), style = 3, file = stderr())
    i = 0

    specificity_scores_cluster_i = apply(average_profile_matrix, 1, function(x) {
      ## Calculate specificity for feature
      if (sum(x) > 0) {
        specificity = 1 - js_distance(make_probability_vector(x), perfect_specificity)
      } else {
        specificity = 0
      }

      # Increment progress bar and return
      i <<- i + 1
      setTxtProgressBar(pb = pb, value = i)
      return(specificity)
    })
    close(pb)
    return(specificity_scores_cluster_i)
  })
  
  marker_gene_specificities = do.call(cbind, marker_gene_specificities)
  colnames(marker_gene_specificities) = group_names[column_indices_to_test]
  return(marker_gene_specificities)
}

# Function to calculate specificity scores. Helpful in finding features that act as best markers of a group/cluster.
# 
# Args:
#  mat (dcGMatrix): sparse, binary matrix of features (rows) by cells (columns)
#  groups (vector): vector assigning each cell in mat to a group
#  groups_to_test (vector): Vector of group names (drawn from categories in group) to return specificity scores for
#    note that all groups are considered when calculating scores, but specificity score calculations only run for specified groups (faster runtime)
#
# Returns:
#  data.frame: 
get_specificity_scores = function(mat, groups, groups_to_test=NULL) {
  # Average expression on binarized matrix is just proportion of site in each cluster
  message('Getting scaled proportions...')
  proportions = get_average_expression(mat, groups)

  # Get median total counts for each cluster for normalization and calculate scaling factor
  cluster_median_depth = get_median_depth(mat, groups)
  scaling_factor = mean(log10(cluster_median_depth))/log10(cluster_median_depth)
  proportions_scaled = t(t(proportions)*scaling_factor)

  # Now use normalized proportions to calculate specificity
  feature_specificities = calculate_specificity_helper(proportions_scaled, groups_to_test)

  # Square specificities and scale by proportion to favor specific sites that still are seen a reasonable proportion of cells
  if(!is.null(groups_to_test)) {
    proportions_scaled = proportions_scaled[, groups_to_test]
  }

  final_specificity_scores = feature_specificities^2 * proportions_scaled
  final_specificity_scores = reshape2::melt(final_specificity_scores)
  colnames(final_specificity_scores) = c('feature', 'group', 'specificity_score')
  return(final_specificity_scores)
}

