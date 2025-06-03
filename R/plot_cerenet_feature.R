#' Plot a CERE-Net Feature
#'
#' Given a Seurat object with computed coordination (or raw feature) scores in \code{@meta.data}, this
#' function produces a spatial scatter plot of any one feature (either raw or propagated) across spots.
#'
#' @param seurat_obj A \code{\link[Seurat]{Seurat}} object containing spatial spots and metadata.
#' @param feature_name Character string naming the metadata column to plot (e.g. "signal_coord" or "signal_raw").
#' @param pt_size Numeric; size of the points on the plot. Default is \code{1.0}.
#' @param colormap Character vector of length \code{2} giving start and end colors (e.g. c("viridis::viridis(10)")). Default is \code{c("viridis")}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object showing the spatial distribution of \code{feature_name}.
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c theme_minimal labs
#' @export
plot_cerenet_feature <- function(seurat_obj, feature_name, pt_size = 1.0, colormap = "viridis") {
  # Extract metadata and spatial coordinates
  df <- seurat_obj@meta.data
  coords <- FetchData(seurat_obj, vars = c("imagerow", "imagecol"))
  df$imagerow <- coords$imagerow
  df$imagecol <- coords$imagecol
  
  if (!feature_name %in% colnames(df)) {
    stop(paste0("'", feature_name, "' is not found in seurat_obj@meta.data"))
  }
  
  ggplot2::ggplot(df, ggplot2::aes(x = imagecol, y = -imagerow, color = .data[[feature_name]])) +
    ggplot2::geom_point(size = pt_size) +
    ggplot2::scale_color_viridis_c(option = colormap) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0("Spatial plot of ", feature_name),
      x = "X Coordinate",
      y = "Y Coordinate",
      color = feature_name
    )
}

