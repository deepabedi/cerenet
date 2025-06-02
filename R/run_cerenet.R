#' One-liner to run CereNet on a Visium folder
#'
#' @param path    Folder with 10x Visium files (outs/).
#' @param modules Named list: each element is a gene vector.
#' @param n_iter  Diffusion iterations (default 5).
#' @export
run_cerenet <- function(path, modules, n_iter = 5) {
  stopifnot(dir.exists(path))
  obj <- Seurat::Load10X_Spatial(path,
                                 filename = "filtered_feature_bc_matrix.h5")
  Seurat::DefaultAssay(obj) <- "Spatial"
  obj <- Seurat::SCTransform(obj, assay = "Spatial", verbose = FALSE) |>
    Seurat::RunPCA(verbose = FALSE) |>
    Seurat::FindNeighbors(dims = 1:20)
  Seurat::DefaultAssay(obj) <- "SCT"
  
  obj <- Seurat::AddModuleScore(obj, features = modules,
                                name = names(modules), n.random = 0)
  score_cols <- grep("1$", tail(colnames(obj[[]]), length(modules)), value = TRUE)
  graph      <- as(obj@graphs$SCT_snn, "dgCMatrix")
  
  obj <- propagate_signal_sparse(obj, score_cols, graph, n_iter = n_iter)
  obj
}
