#' Compute Spearman Coordination Between Raw and Propagated Signals
#'
#' For each pair of raw/propagated columns in a Seurat object's metadata, computes the Spearman
#' correlation (ρ) and stores it in a new metadata column named `{raw}_coord`.
#'
#' @param seurat_obj A \code{\link[Seurat]{Seurat}} object. Its \code{@meta.data} must contain
#'        columns named in \code{raw_cols} and \code{prop_cols}.
#' @param raw_cols   Character vector of column names (in \code{seurat_obj@meta.data}) that hold the
#'        raw signal values (e.g. \code{"signal_raw"}).
#' @param prop_cols  Character vector of column names (in \code{seurat_obj@meta.data}) that hold the
#'        propagated signal values (e.g. \code{"signal_coord"}). Must be the same length as
#'        \code{raw_cols}; each \code{prop_cols[i]} corresponds to \code{raw_cols[i]}.
#'
#' @return The same \code{\link[Seurat]{Seurat}} object, but with additional metadata columns.
#'         For each \code{raw_cols[i]}, a new column named \code{paste0(raw_cols[i], "_coord")}
#'         is created, containing the Spearman correlation ρ between raw and propagated values
#'         (repeated for every spot to make it easy to inspect).
#'
#' @details
#' This function loops over each index \code{i} in \code{seq_along(raw_cols)}, then:
#' 1. Extracts the raw vector \code{seurat_obj@meta.data[[ raw_cols[i] ]]}.
#' 2. Extracts the propagated vector \code{seurat_obj@meta.data[[ prop_cols[i] ]]}.
#' 3. Computes \code{rho <- cor(raw_vector, prop_vector, method="spearman", use="complete.obs")}.
#' 4. Stores \code{rho} in a new metadata column \code{paste0(raw_cols[i], "_coord")} of length
#'    \code{nrow(seurat_obj@meta.data)}, so that each spot has the same \code{rho} value.
#'
#' @importFrom stats cor
#' @export
compute_coordination <- function(seurat_obj, raw_cols, prop_cols) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  if (length(raw_cols) != length(prop_cols)) {
    stop("`raw_cols` and `prop_cols` must be the same length.")
  }
  
  # Access the metadata data.frame
  df <- seurat_obj@meta.data
  n_spots <- nrow(df)
  
  for (i in seq_along(raw_cols)) {
    raw_name  <- raw_cols[i]
    prop_name <- prop_cols[i]
    
    if (!raw_name %in% colnames(df)) {
      stop(paste0("Raw column '", raw_name, "' not found in seurat_obj@meta.data"))
    }
    if (!prop_name %in% colnames(df)) {
      stop(paste0("Propagated column '", prop_name, "' not found in seurat_obj@meta.data"))
    }
    
    raw_vec  <- df[[raw_name]]
    prop_vec <- df[[prop_name]]
    
    # Compute Spearman correlation (single value)
    rho <- stats::cor(raw_vec, prop_vec, method = "spearman", use = "complete.obs")
    
    # Create a new column named paste0(raw_name, "_coord"), repeating the same rho
    new_col_name <- paste0(raw_name, "_coord")
    seurat_obj@meta.data[[new_col_name]] <- rep(rho, n_spots)
  }
  
  return(seurat_obj)
}
