#' Calculate Modularity for each step of flowbca
#'
#' This function calculates the modularity for each step of the
#' Flow based-Clustering Algorithm (flowbca) clustering process.
#'
#' @param unit_set A data frame from the result of flowbca function.
#'        It must contain 'sourceunit', 'destinationunit', and 'round' columns.
#' @param F_matrix_history A list of flow matrices from the flowbca process.
#'        The first matrix in the list (the original OD matrix) is used for the calculation.
#' @param resolution A numeric value for the resolution parameter (gamma). Default is 1.0.
#' @return A data frame with two columns: 'round' and 'modularity',
#'         showing the modularity value for each clustering step.
#' @export
#' @examples
#' # Assuming 'bca_result' is an object from a flowbca function call
#' # modularity_result <- flowbca_modularity(bca_result$unit_set, bca_result$F_matrix_history)
#' # plot(rev(modularity_result$round$round), modularity_result$round$modularity)
flowbca_modularity <- function(unit_set, F_matrix_history, resolution = 1){
  # --- Input Validation ---
  if (!is.data.frame(unit_set) || !all(c("sourceunit", "destinationunit", "round") %in% names(unit_set))) {
    stop("`unit_set` must be a data frame with 'sourceunit', 'destinationunit', and 'round' columns.")
  }
  if (!is.list(F_matrix_history) || length(F_matrix_history) == 0) {
    stop("`F_matrix_history` must be a non-empty list of matrices.")
  }

  # --- Main Logic ---
  F_matrix_1 <- F_matrix_history[[1]]
  cluster_history <- cluster_set(unit_set[match(colnames(F_matrix_1),unit_set$sourceunit),])

  v_modularity <- unlist(lapply(cluster_history, function(x) {
    suppressMessages(calculate_modularity(input_matrix = F_matrix_1, communities = x$clusterid, resolution = resolution))
  }))

  modularity_df <- data.frame(
    round = names(v_modularity),
    modularity = v_modularity,
    row.names = NULL
  )
  return(modularity_df)
}

#' Create a set of cluster assignments for each merge event.
#' @noRd
cluster_set <- function(unit_set){
  unit_set <- as.data.frame(unit_set)
  round_set <- unit_set[!is.na(unit_set$round), ]
  round_set <- round_set[order(round_set$round, decreasing = TRUE), ]

  cluster_history <- list()
  cluster_history[[1]] <- data.frame(sourceunit = unit_set$sourceunit,
                                     clusterid = unit_set$sourceunit,
                                     stringsAsFactors = FALSE)
  rd <- round_set$round
  # Iterate through each merge event by its row index
  for(i in seq_len(nrow(round_set))) {
    merge_info <- round_set[i, ] # Ensures only one merge is processed
    previous_clusters <- cluster_history[[i]]
    new_clusters <- previous_clusters
    new_clusters$clusterid[new_clusters$clusterid == merge_info$sourceunit] <- merge_info$destinationunit
    cluster_history[[i+1]] <- new_clusters
  }

  max_rd <- if(length(rd) > 0) max(rd) else 0
  names(cluster_history) <- c(as.character(max_rd + 1), as.character(rd))
  return(cluster_history)
}


#' Calculate the Modularity Matrix
#'
#' This function calculates the modularity matrix from an input adjacency matrix.
#' It includes validation to ensure the input is a square matrix and that
#' row and column names match if they exist.
#'
#' @param input_matrix A numeric matrix representing the network's adjacency matrix.
#' @param is_directed A boolean indicating if the graph is directed (TRUE) or undirected (FALSE).
#' @param resolution A numeric value for the resolution parameter (gamma). Default is 1.0.
#' @return A modularity matrix with the same dimensions and names as the input matrix.
#' @export

CalcModMatrix <- function(input_matrix, is_directed = TRUE, resolution = 1.0) {

  # --- 1. Input Validation ---

  # Ensure the input is treated as a matrix
  mat <- as.matrix(input_matrix)

  # Check 1: The matrix must be a square matrix.
  if (nrow(mat) != ncol(mat)) {
    stop("Validation Error: The input_matrix is not a square matrix.")
  }

  # Check 2: If names exist, row and column names must be identical.
  # This check is performed only if the matrix has row names.
  if (!is.null(rownames(mat))) {
    if (!identical(rownames(mat), colnames(mat))) {
      stop("Validation Error: The row names and column names of the matrix must be identical.")
    }
  }

  # --- 2. Symmetrize Matrix if Undirected ---

  # If the graph is undirected, the matrix is symmetrized by adding its transpose.
  if (is_directed == FALSE) {
    mat <- mat + t(mat)
  }

  # --- 3. Calculate Modularity Matrix ---

  # Calculate the total weight (sum of all edge weights) in the matrix.
  total_weight <- sum(mat)

  # Edge case: If the matrix is empty (total_weight is 0),
  # return a zero matrix with the same dimensions and names to avoid division by zero.
  if (total_weight == 0) {
    return(matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)))
  }

  # Calculate the out-degree and in-degree for each node.
  out_degree <- rowSums(mat)
  in_degree <- colSums(mat)

  # Calculate the matrix of expected edge weights under the configuration model.
  # E_ij = (out_degree_i * in_degree_j) / total_weight
  E_degree <- outer(out_degree, in_degree) / total_weight

  # Calculate the modularity matrix (B). B_ij = A_ij - gamma * E_ij
  m_modularity_matrix <- mat - (resolution * E_degree)
  
  # Note: Matrix operations in R preserve dimension names, so there's no need
  # to manually re-assign rownames and colnames to m_modularity_matrix.
  # The result will automatically inherit the names from 'mat'.

  # --- 4. Return Result ---
  
  return(m_modularity_matrix)
}

#' Function to calculate modularity for a given community structure.
#'
#' @param input_matrix A numeric adjacency matrix representing the network.
#' @param communities A vector indicating the community membership of each node. The order or names must match the input_matrix.
#' @param is_directed A logical value indicating if the network is directed (TRUE) or should be treated as undirected (FALSE).
#' @param resolution A numeric value for the resolution parameter (gamma). Default is 1.0.
#' @return A single numeric value representing the modularity score (Q) for the given community partition.
#' @export
calculate_modularity <- function(input_matrix, communities, is_directed = TRUE, resolution = 1.0) {

  # --- Input validation logic start ---
  # 1. Check if the number of rows in input_matrix matches the length of the communities vector
  if (nrow(input_matrix) != length(communities)) {
    stop("Error: The number of rows in input_matrix must match the length of the communities vector.")
  }

  # 2. Check whether input_matrix and communities have names
  has_matrix_names <- !is.null(rownames(input_matrix))
  has_community_names <- !is.null(names(communities))

  # 3. Validate based on names or notify user
  if (has_matrix_names && has_community_names) {
    # If both objects have names, check if they match
    if (!all(rownames(input_matrix) == names(communities))) {
      stop("Error: Row names of input_matrix and names of communities do not match.")
    }
  } else {
    # If either of them does not have names, inform the user that order matters
    message("The order of the matrix must match the order of the objects in communities.")
  }
  # --- Input validation logic end ---

  # Original modularity calculation logic
  # Symmetrization happens here exactly once, so CalcModMatrix is called with
  # is_directed = TRUE to avoid symmetrizing the matrix a second time.
  if (is_directed == TRUE) {
    mat <- as.matrix(input_matrix)
  } else {
    mat <- as.matrix(input_matrix) + t(as.matrix(input_matrix))
  }

  M_matrix <- CalcModMatrix(input_matrix = mat, is_directed = TRUE, resolution = resolution)

  total_weight <- sum(mat)
  
  # Avoid division by zero if the graph has no edges
  if (total_weight == 0) {
    return(0)
  }
  
  same_community_matrix <- outer(communities, communities, "==")
  modularity <- sum(M_matrix[same_community_matrix]) / total_weight
  
  return(modularity)
}
