#' Propagate Signal Using a Sparse Connectivity Graph
#'
#' Given an initial signal vector (e.g. gene-module scores or PD1 expression) for each spot and a sparse
#' adjacency matrix (or igraph) representing spot-to-spot connections, perform iterative propagation of
#' the signal across the spatial graph. This is useful for smoothing a discrete signal over tissue
#' architecture in a memory-efficient way.
#'
#' @param seurat_obj A \code{\link[Seurat]{Seurat}} object containing spatial transcriptomics data.
#'        Must have spatial coordinates stored in \code{seurat_obj@meta.data}.
#' @param initial_signal A numeric vector of length equal to the number of spots in \code{seurat_obj}.
#'        This vector holds the raw "mossy" signal for each spot. Order must match \code{colnames(seurat_obj)}.
#' @param graph Either:
#'   - A sparse adjacency matrix (any \code{Matrix} package class, e.g. \code{dgCMatrix}, \code{ddiMatrix}, etc.),
#'   - An \code{igraph} object encoding spot-to-spot connectivity. If using an \code{igraph}, it will be
#'     converted internally to a sparse \code{dgCMatrix}.
#'   The resulting adjacency must be \code{n_spots x n_spots}, where \code{n_spots = ncol(seurat_obj)}.
#' @param n_iter Integer: the number of propagation iterations to perform. Default is \code{10}.
#' @param alpha Numeric between 0 and 1: damping factor controlling how much of the neighbor signal is
#'        retained at each iteration (default = 0.5).
#'
#' @return A numeric vector of length equal to the number of spots. Each element is the propagated
#'         signal score for the corresponding spot after \code{n_iter} iterations.
#'
#' @details
#' Internally, this function:
#' 1. Converts an \code{igraph} object to a \code{dgCMatrix} if necessary.
#' 2. If \code{graph} is any other sparse \code{Matrix} class (e.g. \code{ddiMatrix}), it coerces it to \code{dgCMatrix}.
#' 3. Normalizes the adjacency matrix row-wise so that each row sums to 1.
#' 4. Iteratively updates the signal via:
#'    \preformatted{
#'      new_signal <- alpha * (P %*% old_signal) + (1 - alpha) * old_signal
#'    }
#'    until \code{n_iter} steps have been performed.
#'
#' @importFrom Matrix rowSums Diagonal
#' @importFrom igraph as_adjacency_matrix
#' @importFrom methods as
#' @export
propagate_signal_sparse <- function(seurat_obj,
                                    initial_signal,
                                    graph,
                                    n_iter = 10,
                                    alpha  = 0.5) {
  # 1. Validate inputs
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  n_spots <- ncol(seurat_obj)
  if (!is.numeric(initial_signal) || length(initial_signal) != n_spots) {
    stop("`initial_signal` must be a numeric vector of length = ncol(seurat_obj).")
  }
  
  # 2. Convert igraph to sparse adjacency if needed
  if (inherits(graph, "igraph")) {
    adj <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
  } else if (inherits(graph, "Matrix")) {
    adj <- as(graph, "dgCMatrix")
  } else {
    stop("`graph` must be either an igraph or a Matrix package sparse matrix.")
  }
  
  # 3. Ensure adj is a dgCMatrix with correct dimensions
  if (!inherits(adj, "dgCMatrix")) {
    stop("After conversion, adjacency is not a dgCMatrix. Check your 'graph' input.")
  }
  if (nrow(adj) != n_spots || ncol(adj) != n_spots) {
    stop("Adjacency `graph` must be square with dimension = number of spots (ncol(seurat_obj)).")
  }
  
  # 4. Normalize each row so rowSums(adj_norm) = 1
  row_sums <- Matrix::rowSums(adj)
  if (any(row_sums == 0)) {
    stop("Some rows of adjacency have sum = 0; cannot normalize.")
  }
  inv_sums <- 1 / row_sums
  P <- Matrix::Diagonal(x = inv_sums) %*% adj
  
  # 5. Iterative propagation
  prop <- initial_signal  # length = n_spots
  for (i in seq_len(n_iter)) {
    prop <- alpha * (P %*% prop) + (1 - alpha) * initial_signal
  }
  
  # 6. Return as a plain numeric vector
  result <- as.numeric(prop)
  if (length(result) != n_spots) {
    stop("Internal error: propagated vector length mismatch.")
  }
  return(result)
}



