#' Run the Full CERE-Net Pipeline
#'
#' This function takes a Seurat object, a named raw-signal vector, and a graph,
#' then:
#' 1. Propagates the raw signal across the graph using \code{propagate_signal_sparse()}.
#' 2. Inserts raw and propagated scores into \code{@meta.data}.
#' 3. Computes Spearman coordination values with \code{compute_coordination()}.
#' 4. Optionally produces a ggplot of one specified feature.
#'
#' @param seurat_obj    A \code{\link[Seurat]{Seurat}} object containing spatial data loaded.
#' @param signal_vector A named numeric vector of length = \code{ncol(seurat_obj)}, where names match
#'                      spot IDs (i.e., \code{colnames(seurat_obj)}). These are the raw scores.
#' @param graph         A sparse adjacency matrix (class \code{dgCMatrix}) or \code{igraph} object
#'                      encoding spot-to-spot connectivity. Must match the spots in \code{seurat_obj}.
#' @param alpha         Numeric between 0 and 1; damping factor for propagation (default = 0.5).
#' @param n_iter        Integer; number of propagation iterations (default = 10).
#' @param plot_feat     Character (optional). Name of the feature in \code{@meta.data} to plot after
#'                      propagation (e.g. \code{"signal_coord"}). If \code{NULL}, no plot is generated.
#' @param pt_size       Numeric; point size for the spatial plot (default = 1.0).
#' @param colormap      Character; name of a viridis palette for plotting (default = "viridis").
#'
#' @return A list with components:
#'   \itemize{
#'     \item{\code{seurat_obj}}{The input Seurat object, now with these new \code{@meta.data} columns:}
#'       \itemize{
#'         \item{\code{signal_raw}}{the original \code{signal_vector} repeated per-spot.}
#'         \item{\code{signal_coord}}{the propagated (coordination) score per spot.}
#'         \item{\code{<raw>_coord}}{Spearman ρ per raw/propagated pair (e.g. \code{"signal_raw_coord"}).}
#'       }
#'     \item{\code{plot}}{A \code{\link[ggplot2]{ggplot}} object of \code{plot_feat}, or \code{NULL} if no plot requested.}
#'     \item{\code{rho_values}}{A named numeric vector of the Spearman ρ values for each \code{raw_col}.}
#'   }
#'
#' @details
#' Internally, this function calls:
#' - \code{\link{propagate_signal_sparse}()} to get a numeric \code{signal_coord} vector.
#' - Writes \code{signal_raw} and \code{signal_coord} into \code{seurat_obj@meta.data}.
#' - Calls \code{\link{compute_coordination}()} to compute and store Spearman ρ for each raw/propagated pair.
#' - If \code{plot_feat} is not \code{NULL}, calls \code{\link{plot_cerenet_feature}()} to build a spatial plot.
#'
#' @importFrom Seurat FetchData
#' @importFrom Matrix Diagonal
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c theme_minimal labs
#' @export
run_cerenet <- function(seurat_obj,
                        signal_vector,
                        graph,
                        alpha     = 0.5,
                        n_iter    = 10,
                        plot_feat = NULL,
                        pt_size   = 1.0,
                        colormap  = "viridis") {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  n_spots <- ncol(seurat_obj)
  if (!is.numeric(signal_vector) || length(signal_vector) != n_spots) {
    stop("`signal_vector` must be a numeric vector of length = ncol(seurat_obj).")
  }
  if (!is.null(names(signal_vector)) && !all(names(signal_vector) %in% colnames(seurat_obj))) {
    stop("Names of `signal_vector` must match spot IDs (colnames of seurat_obj).")
  }
  
  # 1. Propagate the raw signal
  propagated <- propagate_signal_sparse(
    seurat_obj     = seurat_obj,
    initial_signal = signal_vector,
    graph          = graph,
    n_iter         = n_iter,
    alpha          = alpha
  )
  
  # 2. Write raw and propagated into metadata
  raw_full   <- rep(NA_real_, n_spots)
  coord_full <- rep(NA_real_, n_spots)
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
  
  # 3. Compute Spearman coordination (ρ) for each raw/propagated pair
  #    In our simple pipeline, we only have one raw column ("signal_raw")
  #    paired with one prop column ("signal_coord"). But compute_coordination can handle multiples.
  raw_cols  <- "signal_raw"
  prop_cols <- "signal_coord"
  
  seurat_obj <- compute_coordination(seurat_obj, raw_cols, prop_cols)
  # Collect the numeric ρ value(s) into a named vector:
  rho_name   <- paste0(raw_cols, "_coord")
  rho_values <- seurat_obj@meta.data[[rho_name]][1]
  
  # 4. Optionally plot one feature
  plot_obj <- NULL
  if (!is.null(plot_feat)) {
    if (!plot_feat %in% colnames(seurat_obj@meta.data)) {
      stop(paste0("`plot_feat`='", plot_feat, "' not found in seurat_obj@meta.data"))
    }
    plot_obj <- plot_cerenet_feature(
      seurat_obj  = seurat_obj,
      feature_name = plot_feat,
      pt_size     = pt_size,
      colormap    = colormap
    )
  }
  
  # 5. Return a list
  return(list(
    seurat_obj = seurat_obj,
    plot       = plot_obj,
    rho_values = rho_values
  ))
}

