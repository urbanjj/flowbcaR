#' Find Hierarchy Path using Base R (Internal Helper)
#'
#' Traverses a parent-child map to find the full hierarchy path for a single node.
#' This is an internal helper function implemented using only base R.
#'
#' @param start_node The node to start the traversal from.
#' @param cluster_map A named vector where names are children and values are parents.
#' @return A string representing the full path.
#' @noRd
find_hierarchy <- function(start_node, cluster_map) {
  path <- c(start_node)
  current_node <- start_node
  visited <- c(start_node)

  while (current_node %in% names(cluster_map)) {
    parent_node <- cluster_map[[current_node]]

    # Stop if the parent is NA or empty
    if (is.na(parent_node)) {
      break
    }

    # Check for circular references
    if (parent_node %in% visited) {
      warning(paste("Circular reference detected:", paste(c(visited, parent_node), collapse=" -> ")))
      path <- c(path, "...") # Indicate a circular path
      break
    }

    path <- c(path, parent_node)
    visited <- c(visited, parent_node)
    current_node <- parent_node
  }

  return(paste(rev(path), collapse = "/"))
}

#' Build and Add Hierarchy Path Column using Base R
#'
#' This function constructs hierarchy paths from a data frame of parent-child
#' relationships using only base R functions. It is a replacement for the
#' original `find_hierarchy_db`.
#'
#' @param data A data frame, typically the `unit_set` from a `flowbca` result.
#' @param child_col A string, the name of the column with child units (defaults to "sourceunit").
#' @param parent_col A string, the name of the column with parent units (defaults to "destinationunit").
#' @return The input data frame with a new `hierarchy` column.
#' @importFrom stats setNames
#' @export
build_hierarchy <- function(data, child_col = "sourceunit", parent_col = "destinationunit") {
  # Ensure input is a data frame
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  # Create the parent-child map from the data frame
  valid_links <- data[!is.na(data[[parent_col]]), ]
  cluster_map <- setNames(valid_links[[parent_col]], valid_links[[child_col]])

  # Use sapply to apply the find_hierarchy function to each child unit
  data$hierarchy <- sapply(data[[child_col]], find_hierarchy, cluster_map = cluster_map)
  data$h_level <- sapply(data$hierarchy, \(x) nchar(x)-nchar(gsub('/','',x)))+1

  # Calculate the number of children for each node
  child_counts <- table(data[[parent_col]])
  data$child_count <- as.integer(child_counts[as.character(data[[child_col]])])
  data$child_count[is.na(data$child_count)] <- 0
  parentid <- data$sourceunit[data$child_count > 0]
  data$h_parent <- sapply(data$sourceunit, \(x) ifelse(x %in% parentid,1,0))

  return(data)
}

#' Build a Nested List Representing the Cluster Hierarchy
#'
#' This function converts the flat parent-child relationship in `unit_set` into a
#' nested list (tree) structure. This is useful for hierarchical analysis or
#' visualization with tree-based packages.
#'
#' @param unit_set A data frame, typically the `unit_set` from a `build_hierarchy` result.
#' @return A named, nested list representing the cluster hierarchy. Each name is a
#'   node, and its value is a list of its children. Leaf nodes are empty lists.
#' @export
#' @examples
#' \dontrun{
#'   # Assuming 'bca_result' is the output from flowbca()
#'   tree <- build_cluster_tree(bca_result$unit_set)
#'   # You can inspect the tree structure
#'   str(tree)
#'   # Access a specific sub-tree
#'   seoul_cluster <- tree$Seoul
#' }
build_cluster_tree <- function(unit_set) {

  # 1. Create a simple, reliable adjacency list (parent -> vector of children)
  valid_links <- unit_set[!is.na(unit_set$destinationunit), ]
  if (nrow(valid_links) == 0) {
    # Handle case with no merges by returning a list of all units as empty lists
    all_units <- unit_set$sourceunit[!is.na(unit_set$sourceunit)]
    return(setNames(lapply(all_units, function(x) list()), all_units))
  }
  adjacency_list <- split(as.character(valid_links$sourceunit), as.character(valid_links$destinationunit))

  # 2. Define the recursive function to build the nested list from the adjacency list
  build_recursive <- function(node_name) {
    children <- adjacency_list[[node_name]]
    if (is.null(children) || length(children) == 0) {
      return(list())  # leaf node: do not include itself as child
    }
    children <- sort(children)
    child_list <- setNames(lapply(children, build_recursive), children)
    # If the node has at least one child, add itself as a child
    child_list[[node_name]] <- list()
    # Sort by name (optional)
    child_list[sort(names(child_list))]
  }

  # 3. Identify the top-level nodes (roots of the trees)
  # These are nodes that are not children of any other node.
  all_units <- unique(c(unit_set$sourceunit, unit_set$destinationunit))
  all_units <- all_units[!is.na(all_units)]
  merged_units <- unit_set$sourceunit[!is.na(unit_set$destinationunit)]
  top_level_nodes <- sort(setdiff(all_units, merged_units))

  # 4. Build the final tree by calling the recursive function on each top-level node
  final_tree <- setNames(lapply(top_level_nodes, build_recursive), top_level_nodes)
  return(final_tree)
}

#' @title Convert a nested list to a dendrogram with labels
#' @description Converts each hierarchy (name) of a list into branch labels of a dendrogram.
#' @param unit_set A data frame, typically the `unit_set` from a `build_hierarchy` result.
#' @return An object of class 'dendrogram'
#' @export
build_dendrogram <- function(unit_set) {

  nested_list <- build_cluster_tree(unit_set)
  
  # Calculate the maximum depth of the list
  get_depth <- function(x) {
    if (!is.list(x) || length(x) == 0) {
      return(0L)
    } else {
      return(1L + max(vapply(x, get_depth, integer(1))))
    }
  }
  
  max_h <- get_depth(nested_list)
  
  # Core function to recursively build dendrogram nodes
  build_node <- function(l, h) {
    node_names <- names(l)
    if (is.null(node_names)) {
      stop("All elements of the list must have names.")
    }
    
    children <- lapply(node_names, function(name) {
      sub_list <- l[[name]]
      
      # Base Case: Handle leaf nodes
      if (length(sub_list) == 0) {
        leaf_node <- 1
        attr(leaf_node, "label") <- name
        attr(leaf_node, "members") <- 1
        attr(leaf_node, "height") <- 0.0
        attr(leaf_node, "leaf") <- TRUE
        class(leaf_node) <- "dendrogram"
        return(leaf_node)
        
      } else {
        # Recursive Step: Handle branch nodes
        subtree <- build_node(sub_list, h - 1)
        attr(subtree, "label") <- name
        attr(subtree, "height") <- h - 1
        return(subtree)
      }
    })
    
    node <- children
    attr(node, "members") <- sum(vapply(node, attr, "members", FUN.VALUE = numeric(1)))
    class(node) <- "dendrogram"
    
    return(node)
  }
  
  dendro <- build_node(nested_list, max_h)
  attr(dendro, "height") <- max_h
  
  return(dendro)
}

#' Analyze and Organize Clusters by Hierarchy
#'
#' @description
#' This function takes a cluster hierarchy (as a nested list) and analyzes its
#' structure at different levels. It converts the tree into a dendrogram, then
#' iteratively cuts it at different heights to identify clusters. For each level,
#' it produces a mapping of units to their respective clusters and identifies
#' "core" clusters (those with more than one member). It also provides information
#' about the parent-child relationships between clusters across hierarchical levels.
#'
#' @param unit_set A data frame, typically the `unit_set` from a `build_hierarchy` result.
#'
#' @return
#' A list containing two named elements:
#' \describe{
#'   \item{unit_set}{A list of data frames. Each data frame corresponds to a
#'     hierarchy level (from the top down) and contains `sourceunit` (the basic
#'     unit), `h_cl` (the cluster it belongs to at that level), and a `core`
#'     flag indicating if its cluster is a core cluster (has more than 1 member).}
#'   \item{cluster_info}{A list of data frames detailing the parent-child
#'     relationships for units that belong to core clusters. It includes the
#'     `sourceunit`, its cluster `h_cl`, and the parent cluster `upper_h_cl`
#'     from the level above.}
#' }
#'
#' @export
hierarchy_cluster <- function(unit_set){

  nested_list <- build_cluster_tree(unit_set)

  flowbca_dendrogram <- build_dendrogram(unit_set)

  h <- attributes(flowbca_dendrogram)$height
  h_nm <- paste0('hierarchy_',h:1)

  h_list <- list()
  for(i in 1:(h-1)) {
    temp <- cut(flowbca_dendrogram, h=i)
    lower_list <- lapply(temp$lower, labels)
    names(lower_list) <- lapply(temp$lower, function(x) attr(x, "label"))
    df <- data.frame(
          sourceunit = unlist(lower_list),
          h_cl = rep(names(lower_list), times = sapply(lower_list, length)),
          row.names = NULL,
          stringsAsFactors = FALSE)
    h_list[[i]] <- df
  }
  names(h_list) <- h_nm[-1]

  cluster_info <- list()
  core_info <- list()
  for(i in 1:(h-1)){
    temp_df <- h_list[[i]]
    
    # Clusters (cores) with more than one member are identified.
    h_cl_counts <- table(temp_df$h_cl)
    multi_member_clusters <- names(h_cl_counts[h_cl_counts > 1])
    df <- temp_df[temp_df$h_cl %in% multi_member_clusters, ]

    # Top-level clusters have no parent: get_parent_node() returns NULL there,
    # which must map to NA (unlist() would silently drop NULLs and misalign p1).
    p1 <- vapply(df$h_cl, function(x) {
      p <- get_parent_node(nested_list, x)
      if (is.null(p)) NA_character_ else p
    }, character(1), USE.NAMES = FALSE)
    cluster_info[[i]] <- cbind(df, upper_h_cl=NA)
    cluster_info[[i]]$upper_h_cl <- p1

    # A lookup table for core clusters is generated.
    core_counts_df <- as.data.frame(table(temp_df$h_cl), stringsAsFactors = FALSE)
    names(core_counts_df) <- c("sourceunit", "core_N")
    df_core <- core_counts_df[core_counts_df$core_N > 1, ]
    if(nrow(df_core) > 0) {
      df_core$core <- 1
    }
    core_info[[i]] <- df_core
  }
  names(cluster_info) <- h_nm[-1]
  names(core_info) <- h_nm[-1]
  
  # The data frames are merged based on the 'sourceunit' (unit ID) from h_list and the 'sourceunit' (cluster name) from core_info.
  unit_set <- Map(function(x,y) {
                                 rl <- merge(x, y, by = 'sourceunit', all.x = TRUE)
                                 rl[is.na(rl)] <- 0
                                 return(rl)
                                 }, h_list, core_info)
  
  hierarchy_cluster <- list(unit_set, cluster_info)
  names(hierarchy_cluster) <- c('unit_set','cluster_info')
  return(hierarchy_cluster)
}
