#' Compute Spearman coordination between raw and propagated signals
#'
#' For each pair of raw/propagated columns, computes Spearman ρ and
#' stores it in `@meta.data` under `{raw}_coord`.
#'
#' @param seurat_obj A Seurat object with `raw_cols` and `prop_cols` in `@meta.data`.
#' @param raw_cols Character vector of raw‐signal column names.
#' @param prop_cols Character vector of propagated‐signal column names.
#' @return The input Seurat object, augmented with `{raw}_coord` columns.
#' @export
compute_coordination <- function(seurat_obj, raw_cols, prop_cols) {
  df <- seurat_obj@meta.data
  for (i in seq_along(raw_cols)) {
    rho <- cor(
      df[[raw_cols[i]]],
      df[[prop_cols[i]]],
      method = "spearman",
      use    = "complete.obs"
    )
    seurat_obj@meta.data[[paste0(raw_cols[i], "_coord")]] <- rho
  }
  seurat_obj
}

