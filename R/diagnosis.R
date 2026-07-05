#' Plot Diagnosis
#'
#' This function plots the diagnosis data.
#'
#' @param diagnosis_data A list containing relative and absolute data.
#' @param y_var Name of the y-axis variable to plot (one of: 'mean', 'min', 'median', 'max', 'intra_flow_ratio', 'inter_flow_ratio', 'g', 'relative_g').
#' @return A plot.
#' @importFrom graphics axis legend lines rect
#' @importFrom grDevices rgb
#' @export
plot_diagnosis <- function(diagnosis_data, y_var) {
  # Extract relative and absolute data
  relative_data <- diagnosis_data$relative
  absolute_data <- diagnosis_data$absolute

  # Determine the range for x and y axes
  x_range <- range(unique(as.integer(c(relative_data$round, absolute_data$round))), na.rm = TRUE)
  y_range <- range(c(relative_data[[y_var]], absolute_data[[y_var]]), na.rm = TRUE)

  # Plot the relative data
  plot(relative_data$round, relative_data[[y_var]],
       type = "l", col = "blue",
       xlim = rev(x_range), ylim = y_range,
       xlab = "Round", ylab = y_var,
       main = paste(y_var, "Diagnosis"))
    
  # Add x-axis with 10-unit intervals
  x_ticks <- seq(from = floor(min(x_range)/10)*10, to = ceiling(max(x_range)/10)*10, by = 10)
  axis(side = 1, at = x_ticks)

  # Add the absolute data to the plot
  lines(absolute_data$round, absolute_data[[y_var]], col = "red")
  
  if(y_var == 'modularity'){
        # Add background color for significant modularity zone (0.3 ~ 0.7)
  usr <- par("usr")  # Get current plot coordinates: c(xmin, xmax, ymin, ymax)
  rect(xleft = usr[1], xright = usr[2],
       ybottom = 0.3, ytop = 0.7,
       col = rgb(173, 216, 230, maxColorValue = 255, alpha = 80),
       border = NA)
    
  # Add a legend
  legend("topleft",
         legend = c("Relative", "Absolute", "Significant Modularity (0.3-0.7)"),
         col = c("blue", "red", rgb(173, 216, 230, maxColorValue = 255)),
         lty = c(1, 1, NA), pch = c(NA, NA, 15),
         pt.cex = 1.2, bty = "n", y.intersp=1.2, cex=0.9)
  } else {
  # Add a legend
  legend("topleft", legend = c("Relative", "Absolute"),
       col = c("blue", "red"), lty = 1, cex=0.9)
  }
}

#' Flow-based Community Analysis Diagnosis
#'
#' This function performs a diagnosis of the flow-based community analysis.
#'
#' @param flow_input The input flow data.
#' @param is_directed A logical value indicating whether the graph is directed
#' @param data_name The name of the data to display in the plot. If NULL, the name of the flow_input object is used.
#' @return A list containing the diagnosis statistics.
#' @importFrom graphics mtext par
#' @export
flowbca_diagnosis <- function(flow_input, is_directed=TRUE, data_name=NULL){
    if (is.null(data_name)) {
        data_name <- deparse(substitute(flow_input)) 
    }
    if(is_directed==TRUE){
        dg_data <- list(flowbca(flow_input,opt_f=1,save_k=TRUE),
                        flowbca(flow_input,opt_f=3,save_k=TRUE))
        names(dg_data) <- c('relative','absolute')
    } else {
        # Symmetrize the flow matrix (F + t(F)) so that the input itself
        # represents an undirected network, as documented. The symmetric input
        # also makes the downstream (directed) modularity formula equal to the
        # undirected one.
        flow_mat <- as.matrix(flow_input[, -1])
        flow_mat <- flow_mat + t(flow_mat)
        flow_input <- data.frame(flow_input[, 1, drop = FALSE], flow_mat,
                                 check.names = FALSE, stringsAsFactors = FALSE)
        dg_data <- list(flowbca(flow_input,opt_f=2,save_k=TRUE),
                        flowbca(flow_input,opt_f=4,save_k=TRUE))
        names(dg_data) <- c('relative','absolute')
    }

    dg_unit_set <- lapply(dg_data, \(x) x[['unit_set']][, c('round', 'g')]) |>
                    lapply(\(x) x[!is.na(x$g), ]) |>
                    lapply(\(x) { x$round <- as.numeric(x$round); x }) |>
                    lapply(\(x) { x$relative_g <- (x$g / max(x$g)); x})
    dg_stat <- lapply(dg_data, \(x) flowbca_stat(x[['F_matrix_history']]))
    dg_moduality <- lapply(dg_data, \(x) flowbca_modularity(x[['unit_set']],x[['F_matrix_history']]))
    
    diagnosis_stat <- Map(function(...) Reduce(function(x, y) merge(x, y, by = "round"), list(...)), dg_stat, dg_moduality, dg_unit_set)

    diag_par <- par(mfrow=c(2,2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0), cex=0.8)
    on.exit(par(diag_par))
    plot_diagnosis(diagnosis_stat, 'mean')
    plot_diagnosis(diagnosis_stat, 'relative_g')
    plot_diagnosis(diagnosis_stat, 'intra_flow_ratio')
    plot_diagnosis(diagnosis_stat, 'modularity')
    mtext(paste0("Diagnosis: ", data_name), outer=TRUE, cex=1.5, line=0)

    invisible(diagnosis_stat)
}
