#' Bright, high-contrast spatial map for a CereNet score
#'
#' @param seurat  Seurat object returned by `run_cerenet()`
#' @param feature Name of the column in `meta.data`
#'                (e.g. "Immune_PD11_prop", "Hypoxia1_prop")
#' @param palette Viridis palette ("inferno", "plasma", "magma", ...)
#' @return        A ggplot2 object
#'
#' @export
plot_cerenet_feature <- function(seurat, feature, palette = "inferno") {
  p <- Seurat::SpatialFeaturePlot(
    seurat,
    features       = feature,
    image.alpha    = 0,      # hide histology
    pt.size.factor = 7,      # big solid dots
    stroke         = 0,
    min.cutoff     = "q05",
    max.cutoff     = "q95"
  ) +
    viridis::scale_fill_viridis(
      option   = palette,
      na.value = "white"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")
  return(p)
}
