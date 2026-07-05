#' Create a GIF Animation of the Clustering Process
#'
#' Generates a GIF animation visualizing the step-by-step merging of spatial units
#' from a `flowbca` analysis. It creates a frame for each round of the clustering
#' process and combines them into a single animated GIF.
#'
#' @param flowbca_gis_layer A named list of `sf` objects, typically the output from the
#'   `flowbca_gis_layer` function.
#' @param unit_set A data frame from the `flowbca` result, containing details
#'   about the merges, including `round`, `sourceunit`, and `clusterid`.
#' @param width The width of the output GIF in pixels. The height is calculated
#'   automatically to maintain the aspect ratio. Defaults to 1000.
#' @param file_nm The file name for the output GIF. "The default value is set to
#'   the name of unit_set, but the user can specify it (e.g., flowbca_animation.gif)."
#' @param keep_frames A logical value. If `TRUE`, the individual PNG frames used to
#'   create the GIF are kept in a sub-directory. Defaults to `FALSE`.
#' @return This function does not return a value but writes a GIF file to the specified
#'   location.
#' @note Requires the `sf` and `magick` packages. For faster and
#'   higher-quality PNG generation, it is recommended to also install the `ragg` package.
#' @importFrom grDevices colors png dev.off
#' @importFrom graphics par text points
#' @export
flowbca_ani <- function(flowbca_gis_layer, unit_set, width = 1000,
                        file_nm = NULL, keep_frames = FALSE) {

  # --- Package Checks ---
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
  if (!requireNamespace("magick", quietly = TRUE)) stop("Package 'magick' is required.")
  use_ragg <- requireNamespace("ragg", quietly = TRUE)
  
  if (is.null(file_nm) || nchar(as.character(file_nm)) == 0) {
    file_nm <- paste0(deparse(substitute(unit_set)),'_ani.gif')
  }

  if (grepl('\\.gif$', file_nm) == FALSE) {
    stop('The file name is required to have the .gif extension.')
  }


  # --- Initial Setup ---
  unit_set$round <- as.character(unit_set$round)
  round_names <- names(flowbca_gis_layer)
  if (length(round_names) < 2) {
    stop("'flowbca_gis_layer' must contain at least two rounds to animate.")
  }
  gis_base <- flowbca_gis_layer[[1]]
  sourceunit <- gis_base$sourceunit

  g_colour <- sample(colors(), length(sourceunit), replace = length(sourceunit) > length(colors()))
  col_db <- data.frame(sourceunit, g_colour)
  gis_png <- lapply(flowbca_gis_layer, \(x) merge(x, col_db, by = 'sourceunit', all.x = TRUE))

  gis_base_center <- sf::st_centroid(gis_base)
  sourceunit_filter <- lapply(flowbca_gis_layer, \(x) sf::st_drop_geometry(x)$sourceunit)
  gis_center <- lapply(sourceunit_filter, \(x) gis_base_center[gis_base_center$sourceunit %in% x, ])

  bbox <- sf::st_bbox(gis_base)
  width_units  <- bbox["xmax"] - bbox["xmin"]
  height_units <- bbox["ymax"] - bbox["ymin"]
  aspect_ratio <- height_units / width_units
  png_width <- width
  png_height <- round(png_width * aspect_ratio)

  # --- Dynamic Scaling Setup ---
  cex_scaler <- width / 1000
  lwd_scaler <- 0.0005 * width + 0.5

  # --- Frame Generation Setup ---
  frame_dir <- paste0(sub("\\.[^.]*$", "", basename(file_nm)), "_frames")
  dir.create(frame_dir, showWarnings = FALSE)

  if (keep_frames) {
    message(paste("Frames will be kept in:", file.path(getwd(), frame_dir)))
  } else {
    on.exit(unlink(frame_dir, recursive = TRUE), add = TRUE)
  }

  # --- Frame Generation Logic ---
  generate_frame <- function(i) {
    j <- round_names[i]
    g_title <- subset(unit_set, round==j)
    g_png <- gis_png[[j]]
    g_png_attention <- subset(g_png,sourceunit==g_title$destinationunit)
    g_center <- gis_center[[j]]
    frame_path <- file.path(frame_dir, sprintf(paste0('%04d_','frame','_r',j,'.png'), i))

    # Use ragg for faster, high-quality PNG generation if available
    if (use_ragg) {
      ragg::agg_png(frame_path, width = png_width, height = png_height)
    } else {
      png(frame_path, width = png_width, height = png_height)
    }
    
    par(mar = c(0, 0, 0, 0))
    plot(sf::st_geometry(g_png), col = g_png$g_colour, lwd = 1)
    text(x = bbox[1], y = bbox[4], labels = paste0('round: ', j), adj = 0, cex = 3 * cex_scaler, col = "blue")
    text(x = bbox[1], y = bbox[4] * 0.99, labels = bquote(.(g_title$sourceunit) %->% .(g_title$destinationunit)), adj = 0, cex = 3 * cex_scaler, col = "blue")
    plot(sf::st_geometry(gis_base), col = NA, border = 'grey50', lwd = 0.5 * lwd_scaler, add = TRUE)
    plot(sf::st_geometry(g_png), col = NA, lwd = 1.5 * lwd_scaler, add = TRUE)
    plot(sf::st_geometry(g_png_attention), col = NA, lwd = 2.5 * lwd_scaler, border = 'red', add = TRUE)
    points(sf::st_coordinates(g_center), type = 'p', pch = 17, col = 'blue', cex = 2 * cex_scaler)
    dev.off()
    
    return(NULL)
  }

  # --- Execute Frame Generation (Sequentially) ---
  message("Generating frames sequentially.")
  lapply(seq_along(round_names)[-1], generate_frame)
  
  # --- GIF Animation ---
  frame_files <- list.files(frame_dir, full.names = TRUE, pattern = "\\.png$")
  if (length(frame_files) == 0) {
      stop("No frames were generated. Check for errors in the plotting code.")
  }
  ani_db <- lapply(frame_files, magick::image_read)
  ani_gif <- magick::image_animate(do.call(c, ani_db), fps = 2, loop = 1)
  magick::image_write(ani_gif, file_nm)
  
  print(paste("Animation saved to:", file_nm))
}