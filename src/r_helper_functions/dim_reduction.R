library(Matrix)
library(stringr)
library(irlba)

########################################
# Functions for LSI
########################################
# Helper function to do fast version of row scaling of sparse TF matrix by IDF vector.
# Exploits the way that data is stored within sparseMatrix object. Seems to be much more memory efficient than tf * idf and faster than DelayedArray.
# Args:
#   tf (sparse matrix): term frequency matrix
#   idf (vector): IDF vector
# Returns:
#   sparse matrix: TF-IDF matrix
safe_tfidf_multiply = function(tf, idf) {
   tf = t(tf)
   tf@x <- tf@x * rep.int(idf, diff(tf@p))
   tf = t(tf)
   return(tf)
}

# Perform TF-IDF on binary matrix
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
#   frequencies (bool): divide bmat by colSums (if FALSE simply use bmat for TF matrix)
#   log_scale_tf (bool): log scale TF matrix if TRUE
#   scale_factor (float): multiply terms in TF matrix by scale_factor prior to log1p. Equivalent to adding small pseudocount but doesn't cast to dense matrix at any point.
# Returns:
#   sparse matrix: TF-IDF matrix
tfidf = function(bmat, frequencies=TRUE, log_scale_tf=TRUE, scale_factor=100000) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  
  # TF-IDF
  tf_idf_counts = safe_tfidf_multiply(tf, idf)
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  return(tf_idf_counts)
}

do_svd = function(tf_idf_counts, dims=50, center=FALSE, scale=FALSE) {
  pca.results = irlba(t(tf_idf_counts), nv=dims, center=center, scale=scale)
  final_result = pca.results$u %*% diag(pca.results$d)
  rownames(final_result) = colnames(tf_idf_counts)
  colnames(final_result) = paste0('PC_', 1:dims)
  return(final_result)
}

filter_features = function(count_matrix, cells=10) {
  count_matrix = count_matrix[Matrix::rowSums(count_matrix) >= cells, ]
  return(count_matrix)
}

filter_cells = function(count_matrix, features_above_threshold=100, threshold=0) {
  count_matrix = count_matrix[, Matrix::colSums(count_matrix > threshold) >= features_above_threshold]
  return(count_matrix)
}

regress_from_pca = function(pca_coords, vector_to_regress) {
  data_df = as.data.frame(covariate=vector_to_regress)

  model_matrix <- Matrix::sparse.model.matrix(
    ~covariate,
    data = metadata,
    drop.unused.levels = TRUE)

  fit <- limma::lmFit(Matrix::t(pca_coords), model_matrix)
  beta <- fit$coefficients[, -1, drop = FALSE]
  beta[is.na(beta)] <- 0
  new_pca_cooords <- Matrix::t(as.matrix(Matrix::t(pca_coords)) -
                                  beta %*% Matrix::t(model_matrix[, -1]))
  return(new_pca_cooords)
}
