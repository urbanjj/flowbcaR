#' Create Step-by-step Merged Spatial Layers
#'
#' Creates a series of merged spatial layers based on the clustering results
#' from the `flowbca` function. It iteratively joins polygons according to the
#' merge rules for each round of the algorithm.
#'
#' @param unit_set A data frame from the `flowbca` result, containing the merge
#'   rules (e.g., `round`, `sourceunit`, `destinationunit`).
#' @param unit_gis An initial `sf` object with polygons corresponding to the units.
#' @param join_col Specifies the join key columns. It can be a single string for a
#'   common column name (e.g., `"ID"`) or a named vector for different names
#'   (e.g., `c("unit_set_ID" = "gis_ID")`). Defaults to `"sourceunit"`.
#' @return A named list of `sf` objects. Each element represents the spatial state
#'   at a specific round of the clustering process, with the names corresponding
#'   to the round number.
#' @note This function requires the `sf` package.
#' @importFrom stats aggregate
#' @export
flowbca_gis_layer <- function(unit_set, unit_gis, join_col = "sourceunit") {

  # --- 1. Input Validation & Column Name Parsing ---
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("The 'sf' package is required for this function.")
  }

  if (is.character(join_col) && !is.null(names(join_col)) && length(join_col) == 1) {
    unit_set_id_col <- names(join_col)[1]
    unit_gis_id_col <- join_col[1]
  } else if (is.character(join_col) && length(join_col) == 1) {
    unit_set_id_col <- join_col
    unit_gis_id_col <- join_col
  } else {
    stop("'join_col' must be a single string or a single named vector.")
  }

  required_set_cols <- c("round", unit_set_id_col, "destinationunit")
  if (!all(required_set_cols %in% names(unit_set))) {
    stop("unit_set must contain columns: ", paste(required_set_cols, collapse = ", "))
  }
  if (!unit_gis_id_col %in% names(unit_gis)) {
    stop(paste("unit_gis must contain the join column:", unit_gis_id_col))
  }

  # --- 2. Initialization & Standardization ---
  all_units_in_set <- unique(unit_set[[unit_set_id_col]])
  gis_filtered <- unit_gis[unit_gis[[unit_gis_id_col]] %in% all_units_in_set, ]

  gis_simple <- gis_filtered[, c(unit_gis_id_col)]
  
  # Standardize the column name to ensure consistent output for downstream functions
  names(gis_simple)[names(gis_simple) == unit_gis_id_col] <- unit_set_id_col

  merge_rules <- unit_set[!is.na(unit_set$round), ]
  Z <- list()
  Z[[1]] <- gis_simple
  rounds_to_merge <- sort(unique(merge_rules$round), decreasing = TRUE)
  current_sf <- gis_simple

  # --- 3. Optimized Main Loop ---
  for (i in seq_along(rounds_to_merge)) {
    round_val <- rounds_to_merge[i]
    rule_for_round <- merge_rules[merge_rules$round == round_val, ]

    if (nrow(rule_for_round) == 0) {
      Z[[i + 1]] <- current_sf
      next
    }

    involved_ids <- unique(c(rule_for_round[[unit_set_id_col]], rule_for_round$destinationunit))
    is_involved <- current_sf[[unit_set_id_col]] %in% involved_ids

    sf_to_process <- current_sf[is_involved, ]
    sf_static <- current_sf[!is_involved, ]

    if (nrow(sf_to_process) == 0) {
        Z[[i + 1]] <- current_sf
        next
    }

    match_indices <- match(sf_to_process[[unit_set_id_col]], rule_for_round[[unit_set_id_col]])
    to_update <- !is.na(match_indices)
    sf_to_process[[unit_set_id_col]][to_update] <- rule_for_round$destinationunit[match_indices[to_update]]

    merged_geoms <- aggregate(
      sf::st_geometry(sf_to_process),
      by = list(ID = sf_to_process[[unit_set_id_col]]),
      FUN = sf::st_union
    )
    names(merged_geoms)[names(merged_geoms) == "ID"] <- unit_set_id_col
    processed_sf <- sf::st_as_sf(merged_geoms)

    geom_col_name <- attr(sf_static, "sf_column")
    if(is.null(geom_col_name) && nrow(sf_static) == 0) {
        geom_col_name <- attr(current_sf, "sf_column")
    }
    sf::st_geometry(processed_sf) <- geom_col_name

    processed_sf <- remove_holes_sf(processed_sf)

    current_sf <- rbind(sf_static, processed_sf)

    Z[[i + 1]] <- current_sf
    message(paste("Processed round:", round_val, "(Step:", i, ")"))
  }

  if (length(rounds_to_merge) > 0) {
    names(Z) <- c(max(rounds_to_merge, na.rm = TRUE) + 1, rounds_to_merge)
  } else {
    # No merges occurred: the single layer is the initial state.
    names(Z) <- as.character(nrow(gis_simple) + 1)
  }

  return(Z)
}

