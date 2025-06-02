#' Sparse-safe propagation of signals over a graph
#'
#' Given a Seurat object and a sparse adjacency graph, iteratively diffuse
#' raw module scores while re-injecting the original values to preserve local peaks.
#'
#' @param seurat_obj  A Seurat object with raw signal columns in `@meta.data`.
#' @param signal_cols Character vector of column names in `@meta.data` to propagate
#'                    (e.g. the columns created by `AddModuleScore()`).
#' @param graph       Sparse adjacency matrix (`dgCMatrix`) or `igraph` object
#'                    matching the spots in `seurat_obj`.
#' @param alpha       Numeric between 0 and 1; weight on propagated vs. raw signal.
#' @param n_iter      Integer; number of diffusion iterations.
#'
#' @importFrom Matrix Diagonal rowSums
#' @export
propagate_signal_sparse <- function(seurat_obj,
                                    signal_cols,
                                    graph,
                                    alpha  = 0.8,
                                    n_iter = 5) {
  
  # 1 extract the raw signal matrix (spots × modules)
  mat <- as.matrix(seurat_obj@meta.data[, signal_cols, drop = FALSE])
  
  # 2 row-normalise the adjacency -> transition probability matrix P
  inv_sum <- 1 / Matrix::rowSums(graph)                # 1/degree for each spot
  P       <- Matrix::Diagonal(x = inv_sum) %*% graph   # ← fully-qualified call
  
  # 3 iterative propagation with re-injection
  prop <- mat
  for (i in seq_len(n_iter)) {
    prop <- alpha * (P %*% prop) + (1 - alpha) * mat
  }
  
  # 4 write results back into the Seurat meta.data
  for (j in seq_along(signal_cols)) {
    seurat_obj[[paste0(signal_cols[j], "_prop")]] <- prop[, j]
  }
  
  return(seurat_obj)
}
