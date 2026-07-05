#' Calculate Internal Flow Statistics (Optimized)
#'
#' Computes summary statistics for internal relative flows from a list of matrices.
#' This optimized version uses efficient base R functions and is well-documented.
#'
#' @param matrix_list A list of square numeric matrices. The list must be named
#'   with the round numbers, which correspond to the number of clusters.
#' @return A data frame with columns: `round`, `mean`, `min`, `median`, `max`.
#' @importFrom stats median
#' @export
flowbca_stat <- function(matrix_list) {
  # sapply is an efficient way to iterate over a list and collect results.
  stats_matrix <- sapply(matrix_list, function(mat) {
    row_sums <- rowSums(mat)
    # Handle cases where row sum is 0 to prevent division by zero (NaN/Inf).
    internal_sums <- ifelse(row_sums == 0, 0, diag(mat))
    internal_flows <- ifelse(row_sums == 0, 0, diag(mat) / row_sums)

    # Return a named vector of summary statistics for each matrix.
    c(
      mean = mean(internal_flows, na.rm = TRUE),
      min = min(internal_flows, na.rm = TRUE),
      median = median(internal_flows, na.rm = TRUE),
      max = max(internal_flows, na.rm = TRUE),
      intra_flow_ratio= sum(internal_sums) / sum((row_sums)),
      inter_flow_ratio= 1 - (sum(internal_sums) / sum((row_sums)))
    )
  })

  # The result from sapply is a matrix with stats in rows and rounds in columns.
  # Transpose it to have rounds in rows, then convert to a data frame.
  stats_df <- as.data.frame(t(stats_matrix))

  # Add round numbers from the list names as a proper column.
  stats_df$round <- as.numeric(rownames(stats_df))

  # Ensure the column order is logical.
  stats_df <- stats_df[, c("round", "mean", "min", "median", "max", "intra_flow_ratio", "inter_flow_ratio")]
  rownames(stats_df) <- NULL # Reset row names for a clean data frame.

  return(stats_df)
}

#' Plot Flow Statistics Over Clustering Rounds (Optimized)
#'
#' Visualizes flow statistics with a fixed y-axis (0-1) and dynamic plot styles.
#' This optimized version is more robust, readable, and handles edge cases.
#'
#' @title Plot Flow Statistics
#' @description
#' This function takes a data frame of flow statistics and creates a side-by-side plot
#' to visualize the results. The first plot shows the mean, min, median, and max
#' internal flow ratios. The second plot shows the intra- and inter-cluster flow ratios.
#'
#' @param stat_data A data frame with flow statistics, typically the output of `flowbca_stat`.
#' @param upper_bound An optional integer. If provided, only data points with a 'round'
#'   value less than or equal to this upper bound will be plotted.
#'
#' @return Invisibly returns `NULL`. This function is called for its side effect of
#'   creating a plot.
#'
#' @importFrom graphics axis grid legend matplot par text
#' @export
#'
#' @examples
#' # The function is designed to be used with `flowbca_run()`,
#' # which returns a list containing the `stats` data frame.
#' # For a standalone example, we first need to generate some data.
#'
#' # Create a list of matrices for demonstration
#' matrix_list <- list(
#'   `10` = matrix(runif(100), 10, 10),
#'   `9` = matrix(runif(81), 9, 9),
#'   `8` = matrix(runif(64), 8, 8),
#'   `7` = matrix(runif(49), 7, 7),
#'   `6` = matrix(runif(36), 6, 6),
#'   `5` = matrix(runif(25), 5, 5),
#'   `4` = matrix(runif(16), 4, 4),
#'   `3` = matrix(runif(9), 3, 3)
#' )
#'
#' # Generate statistics
#' stats <- flowbca_stat(matrix_list)
#'
#' # Plot the statistics
#' flowbca_plot(stats)
#'
#' # Plot with an upper bound
#' flowbca_plot(stats, upper_bound = 7)
flowbca_plot <- function(stat_data, upper_bound = NULL) {
  # --- 1. Input Validation ---
  required_cols <- c("round", "mean", "min", "median", "max", "intra_flow_ratio", "inter_flow_ratio")
  if (!all(required_cols %in% names(stat_data))) {
    stop("Input data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!is.numeric(stat_data$round)) {
    stat_data$round <- as.numeric(as.character(stat_data$round))
  }

  # --- 2. Data Filtering ---
  plot_df <- stat_data
  if (!is.null(upper_bound) && is.numeric(upper_bound)) {
    plot_df <- plot_df[plot_df$round <= upper_bound, ]
  }
  if (nrow(plot_df) == 0) {
    warning("No data to plot after applying the filter.")
    return(invisible(NULL))
  }
  plot_df <- plot_df[order(plot_df$round, decreasing = TRUE), ]

  # --- 3. Dynamic Plot Style & Edge Case Handling ---
  num_points <- nrow(plot_df)
  plot_type <- if (num_points <= 20) "b" else "l"
  plot_pch <- if (num_points <= 20) 19 else NA
  x_values <- plot_df$round
  xlim_range <- if (num_points > 1) {
    c(max(x_values, na.rm = TRUE), min(x_values, na.rm = TRUE))
  } else {
    x_values + c(0.5, -0.5) # Add padding for a single point plot
  }

  # --- 4. Plotting Setup ---
  y_values_1 <- plot_df[, c("mean", "min", "median", "max")]
  y_values_2 <- plot_df[, c("intra_flow_ratio", "inter_flow_ratio")]
  # Use shorter, consistent names for the second plot's legend and labels
  names(y_values_2) <- c("intra", "inter")
  
  plot_colors_1 <- c("blue", "red", "green", "purple")
  plot_colors_2 <- c("orange", "darkgreen")

  # Set up the plotting area for two plots and ensure it's reset on exit
  op <- par(mfrow = c(1, 2))
  on.exit(par(op))

  # --- Plot 1: Internal Relative Flow Statistics ---
  matplot(x = x_values, y = y_values_1, type = plot_type, pch = plot_pch,
          lty = 1, lwd = 2, col = plot_colors_1, xlim = xlim_range, ylim = c(0, 1),
          xaxt = "n",
          xlab = "Round (Number of units)", ylab = "Proportion",
          main = "Statistics of Internal Relative Flow\nby Cluster")
  grid()
  
  # --- 5. Dynamic X-axis Labeling ---
  axis_ticks <- if (diff(range(x_values, na.rm = TRUE)) > 20) {
    seq(from = floor(min(x_values, na.rm = TRUE)/10)*10, to = ceiling(max(x_values, na.rm = TRUE)/10)*10, by = 10)
  } else {
    pretty(x_values)
  }
  axis(1, at = axis_ticks[axis_ticks == floor(axis_ticks)], labels = axis_ticks[axis_ticks == floor(axis_ticks)])

  if (num_points <= 20) {
    mapply(function(y_col, color) {
      y_vals <- plot_df[[y_col]]
      text(x = x_values, y = y_vals, labels = round(y_vals, 3),
           col = color, pos = 3, cex = 0.75)
    }, names(y_values_1), plot_colors_1)
  }
  legend("topleft", legend = names(y_values_1), col = plot_colors_1, lwd = 2, bty = "n")

  # --- Plot 2: Intra- and Inter-Cluster Flow Proportions ---
  matplot(x = x_values, y = y_values_2, type = plot_type, pch = plot_pch,
          lty = 1, lwd = 2, col = plot_colors_2, xlim = xlim_range, ylim = c(0, 1),
          xaxt = "n",
          xlab = "Round (Number of units)", ylab = "Proportion",
          main = "Proportion of Intra- and Inter-Cluster\nFlows in Overall Flow")
  grid()
  axis(1, at = axis_ticks[axis_ticks == floor(axis_ticks)], labels = axis_ticks[axis_ticks == floor(axis_ticks)])

  # --- 6. Find and Plot Intersection ---
  intersection <- find_intersection(plot_df)
  if (!is.null(intersection)) {
    points(intersection$round, intersection$ratio, pch = 19, col = "red", cex = 1.5)
    text(intersection$round, intersection$ratio, 
         labels = paste0("Round: ", round(intersection$round, 2), "\nRatio: ", round(intersection$ratio, 2)),
         pos = 4, col = "red")
  }

  if (num_points <= 20) {
    mapply(function(y_col, color) {
      y_vals <- y_values_2[[y_col]]
      text(x = x_values, y = y_vals, labels = round(y_vals, 3),
           col = color, pos = 3, cex = 0.75)
    }, names(y_values_2), plot_colors_2)
  }
  legend("topleft", legend = names(y_values_2), col = plot_colors_2, lwd = 2, bty = "n")

  return(invisible(NULL))
}

#' Find Intersection of Intra- and Inter-flow Ratios
#'
#' This helper function calculates the intersection point of the 
#' `intra_flow_ratio` and `inter_flow_ratio` curves using linear interpolation.
#'
#' @param plot_df A data frame containing the columns `round`, `intra_flow_ratio`, and `inter_flow_ratio`.
#' @return A list containing the `round` and `ratio` of the intersection point, or `NULL` if no intersection is found.
#' @noRd
find_intersection <- function(plot_df) {
  if (nrow(plot_df) < 2) return(NULL)
  for (i in seq_len(nrow(plot_df) - 1)) {
    p1_intra <- plot_df$intra_flow_ratio[i]
    p2_intra <- plot_df$intra_flow_ratio[i+1]
    p1_inter <- plot_df$inter_flow_ratio[i]
    p2_inter <- plot_df$inter_flow_ratio[i+1]

    if ((p1_intra > p1_inter && p2_intra < p2_inter) || (p1_intra < p1_inter && p2_intra > p2_inter)) {
      # Linear interpolation to find the intersection point
      # (y1_intra - y1_inter) + (x - x1) * ((y2_intra - y1_intra)/(x2 - x1) - (y2_inter - y1_inter)/(x2 - x1)) = 0
      x1 <- plot_df$round[i]
      x2 <- plot_df$round[i+1]
      y1_intra <- p1_intra
      y2_intra <- p2_intra
      y1_inter <- p1_inter
      y2_inter <- p2_inter

      denominator <- (y2_intra - y1_intra) - (y2_inter - y1_inter)
      if (denominator == 0) next # Parallel lines

      t <- (y1_inter - y1_intra) / denominator
      round_val <- x1 + t * (x2 - x1)
      ratio_val <- y1_intra + t * (y2_intra - y1_intra)

      return(list(round = round_val, ratio = ratio_val))
    }
  }
  return(NULL)
}