#' A Flow-based Cluster Algorithm in R
#'
#' This function is an R translation of the flowbca.ado Stata package
#' by Jordy Meekes (2018). It implements a hierarchical clustering algorithm
#' for flow data, including the detailed tie-breaking rules from the original
#' Mata implementation.
#'
#' @param flow_input A data frame where the first column contains source unit identifiers
#'   and subsequent columns represent the flows to destination units. The destination
#'   columns must be in the same order as the source unit rows.
#' @param opt_f An integer specifying the optimization function:
#'   1: Directed relative flows (default)
#'   2: Undirected relative flows
#'   3: Directed absolute flows
#'   4: Undirected absolute flows
#' @param q A numeric flow threshold (can be relative or absolute) that will
#'   adjust the default stopping criterion. The algorithm stops if the maximum
#'   flow for merging is below this value. Defaults to NULL (not used).
#' @param k An integer specifying the desired number of distinct clusters. The
#'   algorithm stops when this number of clusters is reached. Defaults to NULL;
#'   if no stopping condition is specified at all, `k = 1` is used.
#' @param la A numeric value for the average of the internal relative flows of all
#'   clusters. The algorithm terminates if the calculated average (La) exceeds
#'   this value. Defaults to NULL (not used).
#' @param lw A numeric value for the weighted average of the internal relative flows
#'   of all clusters. The algorithm terminates if the calculated weighted average
#'   (Lw) exceeds this value. Defaults to NULL (not used).
#' @param lm A numeric value for the minimum internal relative flow. The algorithm
#'   terminates if the calculated minimum (Lm) exceeds this value. Defaults to
#'   NULL (not used).
#' @param smaller A numeric value between 0 and 1. If the destination unit's
#'   inter-flow is smaller than the source unit's inter-flow multiplied by this
#'   ratio, the merge is stopped. Defaults to NULL (off). 
#'   Regardless of the value of opt_f, the inter-flow size is evaluated based on
#'   absolute values.
#' @param non_zero A logical value. If TRUE, the algorithm stops when all potential
#'   source units have non-zero inter-flows. Defaults to FALSE.
#' @param save_k Specifies whether to return the F_matrix and C_matrix for all
#'   iterations. Defaults to FALSE.
#'
#' @details Only one stopping condition (`q`, `k`, `la`, `lw`, `lm`, `smaller`,
#'   or `non_zero`) can be specified at a time; specifying more than one is an
#'   error. When none is specified, the algorithm runs until a single cluster
#'   remains (`k = 1`).
#'
#' @return A list containing:
#'   - `unit_set`: Details of cluster assignment for each unit.
#'   - `cluster_set`: Statistics for the final clusters.
#'   - `F_matrix`: The final aggregated flow matrix.
#'   - `F_matrix_history`: (Optional) A list of matrices from each clustering round.
#'   - `C_matrix_history`: (Optional) A list of transformation matrices.
#' 
#' @importFrom stats weighted.mean
#' @export
flowbca <- function(flow_input, opt_f = 1, q = NULL, k = NULL, la = NULL, lw = NULL, lm = NULL, smaller = NULL, non_zero = NULL, save_k = FALSE) {

  # --- 0. Stop condition validation ---
  stop_conditions <- c(!is.null(q), !is.null(k), !is.null(la), !is.null(lw), !is.null(lm), !is.null(smaller), !is.null(non_zero))
  if (sum(stop_conditions) > 1) {
    stop("Only one stopping condition (q, k, la, lw, lm, smaller, or non_zero) can be specified at a time.")
  }

  # Set default values if no stop condition is provided
  if (sum(stop_conditions) == 0) {
    k <- 1
  }

  if (!is.null(smaller) && (smaller < 0 || smaller > 1)) {
    stop("'smaller' must be between 0 and 1.")
  }

  # --- 1. Initial Setup ---
  source_units <- flow_input[[1]]
  destination_units <- colnames(flow_input)[-1]
  
  if (length(source_units) != length(unique(source_units))) {
    stop("Source unit IDs must be unique names.")
  }

  if (length(destination_units) != length(unique(destination_units))) {
    stop("Destination unit IDs must be unique names.")
  }

  if (!identical(as.character(source_units), destination_units)) {
    stop("Source unit IDs must be identical to and in the same order as destination column names.")
  }

  F_matrix <- as.matrix(flow_input[, -1])
  rownames(F_matrix) <- colnames(F_matrix) <- source_units
  
  if (nrow(F_matrix) != ncol(F_matrix)) {
    stop("The number of source units must equal the number of destination units.")
  }
  if (!is.null(k) && (k < 1 || k > nrow(F_matrix))) {
    stop(paste("k must be between 1 and", nrow(F_matrix)))
  }

  merge_history <- list()
  F_matrix_history <- list()
  C_matrix_history <- list()
  
  # --- 2. Main Clustering Loop ---
  first_try <- TRUE
  while (TRUE) {
    
    if (!is.null(k) && nrow(F_matrix) <= k) {
      message("Stopping: The number of units in the F_matrix is less than k.")
      break
    }

    # --- 2z. Check custom stopping conditions ---
    if (isTRUE(non_zero) && all(diag(F_matrix) > 0)) {
      message("Stopping: All source units have non-zero inter-flows.")
      break
    }

    K <- nrow(F_matrix)
    current_ids <- rownames(F_matrix)

    # --- 2a. Prepare Search Matrix (G_matrix) and Flow Matrix (F_prime) ---
    search_matrix <- NULL
    F_prime <- F_matrix # Base matrix for tie-breaking absolute flows
    
    if (opt_f %in% c(1, 2)) { # Relative flows
      row_sums <- rowSums(F_matrix)
      row_sums[row_sums == 0] <- 1 
      G_matrix <- F_matrix / row_sums
      
      if (opt_f == 1) { # Directed relative
        search_matrix <- G_matrix
      } else { # Undirected relative
        search_matrix <- G_matrix + t(G_matrix)
      }
      F_prime <- G_matrix # Use relative flows for tie-breaking
    } else { # Absolute flows
      if (opt_f == 3) { # Directed absolute
        search_matrix <- F_matrix
      } else { # Undirected absolute
        search_matrix <- F_matrix + t(F_matrix)
      }
    }
    
    diag(search_matrix) <- -Inf 
    
    # --- 2b. Identify Units to Merge (r and s) ---
    if (max(search_matrix, na.rm = TRUE) == -Inf) {
      message("Stopping: All values in search_matrix are NA, suggesting that no further clustering is possible as all source units have already been merged.")
      break
    }

    max_flow <- max(search_matrix, na.rm = TRUE)
    
    if (!is.null(q) && max_flow < q) {
      message("Stopping: Maximum flow is below threshold q.")
      break
    }
    
    indices <- which(search_matrix == max_flow, arr.ind = TRUE)
    
    # --- 2c. Disambiguation for non-unique (r,s) pairs (Mata caveats) ---
    if (nrow(indices) > 1) {
      cand_r <- sort(unique(indices[, "row"]))
      cand_s <- sort(unique(indices[, "col"]))
      
      # Use directed flows for tie-breaking in undirected cases.
      # The diagonal (internal flows) is zeroed, not set to -Inf: the caveat
      # logic below sums candidate sub-matrices, and a -Inf diagonal would make
      # every sum -Inf so the tie-breaking rules could never fire.
      tie_break_matrix <- if (opt_f %in% c(2, 4)) F_matrix else F_prime
      diag(tie_break_matrix) <- 0

      # Caveat 1 & 2: One dimension is unique, the other is not.
      if (length(cand_r) > 1 && length(cand_s) == 1) {
        # Multiple r, one s. Find r with max flow from other r's.
        sub_matrix <- tie_break_matrix[cand_r, cand_r]
        if(sum(sub_matrix) > 0) {
          col_sums <- colSums(sub_matrix)
          r_idx <- cand_r[which.max(col_sums)]
        } else {
          r_idx <- cand_r[1] # Fallback
        }
        s_idx <- cand_s[1]
      } else if (length(cand_r) == 1 && length(cand_s) > 1) {
        # One r, multiple s. Find s with max flow from other s's.
        sub_matrix <- tie_break_matrix[cand_s, cand_s]
        if(sum(sub_matrix) > 0) {
          col_sums <- colSums(sub_matrix)
          s_idx <- cand_s[which.max(col_sums)]
        } else {
          s_idx <- cand_s[1] # Fallback
        }
        r_idx <- cand_r[1]
      } else { # Caveat 3 & 4: Multiple r and multiple s, or complex ties
        # First, try to select the best 's' based on flows between candidate s's
        s_sub_matrix <- tie_break_matrix[cand_s, cand_s]
        if (sum(s_sub_matrix) > 0) {
          s_col_sums <- colSums(s_sub_matrix)
          s_idx <- cand_s[which.max(s_col_sums)]
        } else {
          # Fallback: if no internal flows, just pick the first one
          s_idx <- cand_s[1]
        }
        
        # Now, given the chosen 's', find the best 'r' that flows to it
        possible_r_for_s <- indices[indices[, "col"] == s_idx, "row"]
        if (length(possible_r_for_s) == 1) {
          r_idx <- possible_r_for_s
        } else { # Still multiple r's for our chosen s, apply caveat 1 logic
          r_sub_matrix <- tie_break_matrix[possible_r_for_s, possible_r_for_s]
          if(sum(r_sub_matrix) > 0) {
            r_col_sums <- colSums(r_sub_matrix)
            r_idx <- possible_r_for_s[which.max(r_col_sums)]
          } else {
            r_idx <- possible_r_for_s[1] # Fallback
          }
        }
      }
    } else {
      r_idx <- indices[1, "row"]
      s_idx <- indices[1, "col"]
    }

    # --- 2z. Check smaller condition ---
    if (!is.null(smaller)) {
      internal_r <- F_matrix[r_idx, r_idx]
      internal_s <- F_matrix[s_idx, s_idx]
      if (internal_s < internal_r * smaller) {
        message("Stopping: Destination inter-flow is smaller than source inter-flow * ratio.")
        break
      }
    }

    # --- 2d. Check Other Stopping Conditions (la, lw, lm) ---
    if (!is.null(la) || !is.null(lw) || !is.null(lm)) {
      row_flows_stop <- rowSums(F_matrix)
      internal_flows <- diag(F_matrix) / row_flows_stop
      internal_flows[is.nan(internal_flows) | is.infinite(internal_flows)] <- 0

      La <- mean(internal_flows, na.rm = TRUE)
      Lw <- weighted.mean(internal_flows, row_flows_stop, na.rm = TRUE)
      Lm <- min(internal_flows, na.rm = TRUE)
      
      if ((!is.null(la) && la <= La) || (!is.null(lw) && lw <= Lw) || (!is.null(lm) && lm <= Lm)) {
        message("Stopping: Condition (la, lw, or lm) met.")
        break
      }
    }
    
    # --- 2e. Aggregate Clusters using Transformation Matrix C ---
    r_id <- current_ids[r_idx]
    s_id <- current_ids[s_idx]
    
    merge_history[[length(merge_history) + 1]] <- list(
      round = K,
      r = r_id, 
      s = s_id, 
      q = max_flow
    )

    C <- diag(K)
    C[r_idx, s_idx] <- 1
    C <- C[, -r_idx, drop = FALSE]

    if (save_k && first_try) {
      F_matrix_history[[1]] <- list(F_matrix)
      first_try <- FALSE
    }

    F_matrix <- t(C) %*% F_matrix %*% C
    new_ids <- current_ids[-r_idx]
    rownames(F_matrix) <- colnames(F_matrix) <- new_ids

    if (save_k) {
      C_matrix <- C
      rownames(C_matrix) <- current_ids
      colnames(C_matrix) <- new_ids

      F_matrix_history[[length(F_matrix_history) + 1]] <- list(
        F_matrix
      )
      C_matrix_history[[length(C_matrix_history) + 1]] <- list(
        C_matrix
      )
    }

  }
  
  # --- 3. Generate Final Result Sets ---
  
  # --- 3a. Create unit_set ---
  unit_set <- data.frame(
    sourceunit = source_units,
    destinationunit = NA,
    clusterid = as.character(source_units),
    g = NA,
    round = NA,
    stringsAsFactors = FALSE
  )
  
  for (merge in rev(merge_history)) {
    is_source_unit <- unit_set$sourceunit == merge$r
    unit_set$destinationunit[is_source_unit & is.na(unit_set$destinationunit)] <- merge$s
    unit_set$g[is_source_unit & is.na(unit_set$g)] <- merge$q
    unit_set$round[is_source_unit & is.na(unit_set$round)] <- (as.numeric(merge$round) - 1)
  }
  
  # Find the final cluster for each unit
  for (i in 1:nrow(unit_set)) {
    current_unit <- unit_set$sourceunit[i]
    dest <- unit_set$destinationunit[i]
    while(!is.na(dest)) {
      current_unit <- dest
      dest_row <- unit_set[unit_set$sourceunit == current_unit, ]
      dest <- dest_row$destinationunit
    }
    unit_set$clusterid[i] <- current_unit
  }

  unit_set$core <- ifelse(is.na(unit_set$g), 1, 0)
   
  unit_set_order <- unit_set[order(unit_set$round, decreasing = TRUE, na.last = TRUE), ]
  
  # --- 3b. Create cluster_set ---
  final_clusters <- rownames(F_matrix)
  
  row_flows <- rowSums(F_matrix)
  internal <- diag(F_matrix)
  
  internal_relative <- ifelse(row_flows == 0, 0, internal / row_flows)
  
  if (length(F_matrix_history) > 0) {
    F_matrix_history <- lapply(F_matrix_history, \(x) x[[1]])
    names(F_matrix_history) <- (max(unit_set$round,na.rm=TRUE)+1):min(unit_set$round,na.rm=TRUE)

    C_matrix_history <- lapply(C_matrix_history, \(x) x[[1]])
    names(C_matrix_history) <- max(unit_set$round,na.rm=TRUE):min(unit_set$round,na.rm=TRUE)
  }

  cluster_set <- data.frame(
    clusterid = final_clusters,
    internal = internal,
    rowflows = row_flows,
    internal_relative = internal_relative,
    stringsAsFactors = FALSE
  )
  
  if(nrow(cluster_set) > 0) {
      cluster_set$La = mean(internal_relative)
      cluster_set$Lw = weighted.mean(internal_relative, w = row_flows)
      cluster_set$Lm = min(internal_relative)
      cluster_set$N = sum(row_flows)
  }

  if(save_k == TRUE) {
    return(list(unit_set = unit_set_order, cluster_set = cluster_set,
            F_matrix=F_matrix, F_matrix_history = F_matrix_history,
            C_matrix_history = C_matrix_history))
  } else {
    return(list(unit_set = unit_set_order, cluster_set = cluster_set,
            F_matrix=F_matrix))
  }
}