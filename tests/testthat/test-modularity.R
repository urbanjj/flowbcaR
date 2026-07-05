# Two triangles (1-2-3, 4-5-6) joined by a single bridge (3-4): a textbook
# two-community graph used as the reference network throughout.
two_triangles <- function() {
  A <- matrix(0, 6, 6)
  edges <- rbind(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6), c(4, 6), c(3, 4))
  for (k in seq_len(nrow(edges))) {
    A[edges[k, 1], edges[k, 2]] <- 1
    A[edges[k, 2], edges[k, 1]] <- 1
  }
  A
}

# Reference implementation of Q for a (possibly directed) weighted network:
# Q = (1/W) * sum_ij (W_ij - gamma * s_i^out * s_j^in / W) * delta(C_i, C_j)
reference_modularity <- function(A, comm, gamma = 1) {
  W <- sum(A)
  E <- outer(rowSums(A), colSums(A)) / W
  same <- outer(comm, comm, '==')
  sum((A - gamma * E)[same]) / W
}

test_that("CalcModMatrix validates its input", {
  expect_error(CalcModMatrix(matrix(1, 2, 3)), "square")
  m <- matrix(1, 2, 2, dimnames = list(c('a', 'b'), c('a', 'x')))
  expect_error(CalcModMatrix(m), "identical")
  z <- matrix(0, 3, 3)
  expect_equal(CalcModMatrix(z), z)
})

test_that("CalcModMatrix matches the configuration-model formula", {
  A <- two_triangles()
  B <- CalcModMatrix(A, is_directed = TRUE)
  expect_equal(B, A - outer(rowSums(A), colSums(A)) / sum(A))
  # resolution scales only the expected-weight term
  B2 <- CalcModMatrix(A, is_directed = TRUE, resolution = 2)
  expect_equal(B2, A - 2 * outer(rowSums(A), colSums(A)) / sum(A))
})

test_that("undirected modularity equals the Newman reference (no doubling)", {
  A <- two_triangles()
  comm <- c(1, 1, 1, 2, 2, 2)
  Q_ref <- reference_modularity(A, comm)
  expect_equal(suppressMessages(calculate_modularity(A, comm, is_directed = FALSE)),
               Q_ref)
  # for a symmetric matrix the directed formula must agree
  expect_equal(suppressMessages(calculate_modularity(A, comm, is_directed = TRUE)),
               Q_ref)
})

test_that("directed modularity matches the Arenas/Leicht-Newman formula", {
  A <- matrix(c(0, 3, 0, 1,
                0, 0, 2, 0,
                4, 0, 0, 0,
                0, 0, 1, 0), 4, 4, byrow = TRUE)
  comm <- c(1, 1, 1, 2)
  expect_equal(suppressMessages(calculate_modularity(A, comm, is_directed = TRUE)),
               reference_modularity(A, comm))
  # resolution parameter
  expect_equal(suppressMessages(calculate_modularity(A, comm, is_directed = TRUE,
                                                     resolution = 1.5)),
               reference_modularity(A, comm, gamma = 1.5))
})

test_that("calculate_modularity validates names and lengths", {
  A <- two_triangles()
  expect_error(suppressMessages(calculate_modularity(A, c(1, 2))), "length")
  rownames(A) <- colnames(A) <- letters[1:6]
  comm <- setNames(c(1, 1, 1, 2, 2, 2), letters[c(1:5, 26)])
  expect_error(calculate_modularity(A, comm), "do not match")
})

test_that("flowbca_modularity returns one value per clustering state", {
  fi <- make_flow_input(toy_flow_matrix())
  res <- suppressMessages(flowbca(fi, k = 2, save_k = TRUE))
  md <- flowbca_modularity(res$unit_set, res$F_matrix_history)
  n_merges <- sum(!is.na(res$unit_set$round))
  expect_equal(nrow(md), n_merges + 1)
  expect_true(all(md$modularity >= -1 & md$modularity <= 1))
  # the fully-singleton partition has (near-)zero modularity at most
  first <- md$modularity[1]
  expect_lt(abs(first), 1e-8 + 0.5)
})

test_that("cluster_set helper handles a unit_set with no merges", {
  us <- data.frame(sourceunit = c('a', 'b'),
                   destinationunit = NA_character_,
                   round = NA_integer_,
                   stringsAsFactors = FALSE)
  ch <- flowbcaR:::cluster_set(us)
  expect_length(ch, 1)
  expect_equal(ch[[1]]$clusterid, c('a', 'b'))
})
