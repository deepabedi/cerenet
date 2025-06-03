#' Compute and Correlate Coordination Scores
#'
#' 1. Propagate a "mossy" signal across spatial spots (using a sparse graph).
#' 2. Store both raw and propagated values in \code{@meta.data}.
#' 3. Compute a single Spearman \eqn{\rho} across all spots between raw and propagated signals,
#'    and also write that same \eqn{\rho} back into metadata (so the user can inspect it).
#'
#' @param seurat_obj    A Seurat object with spatial data already loaded.
#' @param signal_vector A named numeric vector of raw signal values (names must match spot IDs).
#' @param graph         A sparse adjacency matrix (class \code{dgCMatrix}) or \code{igraph}
#'                      object encoding spot-to-spot connectivity. Must match \code{seurat_obj}.
#' @param alpha         Numeric between 0 and 1 for diffusion damping; default is 0.5.
#' @param n_iter        Integer for number of propagation iterations; default is 10.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item{\code{seurat_obj}}{The input Seurat object, but now with three new columns in \code{@meta.data}:}
#'       \itemize{
#'         \item{\code{signal_raw}}{the original \code{signal_vector} repeated per-spot.}
#'         \item{\code{signal_coord}}{the propagated (coordination) score per spot.}
#'         \item{\code{signal_rho}}{the single Spearman \eqn{\rho}, repeated for every spot (so you can inspect it easily).}
#'       }
#'     \item{\code{rho}}{A single numeric value: Spearman \eqn{\rho} between \code{signal_raw} and \code{signal_coord}.}
#'   }
#'
#' @details
#' Internally, this function calls \code{\link{propagate_signal_sparse}()} to do the graph-based diffusion.
#' Once you have the vector of propagated scores, it computes
#' \eqn{ \rho = \text{cor}( \text{raw}, \text{coord}, \text{method} = "spearman" ) }.
#' That one number is stored back in \code{seurat_obj@meta.data\$signal_rho} for convenience.
#'
#' @importFrom Matrix rowSums Diagonal
#' @importFrom igraph as_adjacency_matrix
#' @importFrom methods as
#' @export
#' @examples
#' \dontrun{
#'   library(CereNet)
#'   library(Seurat)
#'   library(Matrix)
#'   # Create a tiny Seurat object:
#'   mat <- matrix(runif(6), nrow = 2, ncol = 3)
#'   rownames(mat) <- c("Gene1", "Gene2")
#'   colnames(mat) <- c("SpotA", "SpotB", "SpotC")
#'   so <- CreateSeuratObject(counts = mat)
#'   # Build a simple signal vector and identity graph:
#'   sig <- runif(3); names(sig) <- colnames(so)
#'   graph <- Diagonal(3)
#'   # Run function:
#'   out <- compute_coord_and_rho(so, sig, graph)
#'   print(out$rho)
#' }
compute_coord_and_rho <- function(seurat_obj,
                                  signal_vector,
                                  graph,
                                  alpha  = 0.5,
                                  n_iter = 10) {
  # 1. Propagate
  propagated <- propagate_signal_sparse(
    seurat_obj     = seurat_obj,
    initial_signal = signal_vector,
    graph          = graph,
    n_iter         = n_iter,
    alpha          = alpha
  )
  
  # 2. Compute a single Spearman rho across all spots
  #    If signal_vector is named, match names to columns
  if (!is.null(names(signal_vector))) {
    common_spots <- intersect(names(signal_vector), colnames(seurat_obj))
    if (length(common_spots) != length(signal_vector)) {
      warning("Some names in 'signal_vector' do not match spots in 'seurat_obj'. Dropping mismatches.")
    }
    signal_vector <- signal_vector[common_spots]
    propagated    <- propagated[match(common_spots, colnames(seurat_obj))]
  }
  
  rho <- stats::cor(
    x      = as.numeric(signal_vector),
    y      = as.numeric(propagated),
    method = "spearman",
    use    = "complete.obs"
  )
  
  # 3. Write everything back into seurat_obj@meta.data
  raw_full   <- rep(NA_real_, ncol(seurat_obj))
  coord_full <- rep(NA_real_, ncol(seurat_obj))
  rho_full   <- rep(rho, ncol(seurat_obj))
  
  if (!is.null(names(signal_vector))) {
    idx <- match(names(signal_vector), colnames(seurat_obj))
    raw_full[idx]   <- as.numeric(signal_vector)
    coord_full[idx] <- as.numeric(propagated)
  } else {
    raw_full   <- as.numeric(signal_vector)
    coord_full <- as.numeric(propagated)
  }
  
  seurat_obj@meta.data$signal_raw   <- raw_full
  seurat_obj@meta.data$signal_coord <- coord_full
  seurat_obj@meta.data$signal_rho   <- rho_full
  
  # 4. Return list
  return(list(
    seurat_obj = seurat_obj,
    rho        = rho
  ))
}


