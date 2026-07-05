test_that("input validation catches malformed inputs", {
  m <- toy_flow_matrix()

  dup <- make_flow_input(m)
  dup$id[2] <- 'a'
  expect_error(flowbca(dup), "unique")

  mismatch <- make_flow_input(m)
  mismatch$id <- c('a', 'b', 'c', 'x')
  expect_error(flowbca(mismatch), "identical")

  expect_error(flowbca(make_flow_input(m), k = 99), "k must be between")
})

test_that("only one stopping condition may be specified", {
  fi <- make_flow_input(toy_flow_matrix())
  expect_error(flowbca(fi, k = 2, q = 0.5), "Only one stopping condition")
  expect_error(flowbca(fi, la = 0.5, non_zero = TRUE), "Only one stopping condition")
})

test_that("smaller outside [0, 1] is an error", {
  fi <- make_flow_input(toy_flow_matrix())
  expect_error(flowbca(fi, smaller = 1.5), "between 0 and 1")
  expect_error(flowbca(fi, smaller = -0.1), "between 0 and 1")
})

test_that("default run merges down to a single cluster", {
  fi <- make_flow_input(toy_flow_matrix())
  res <- suppressMessages(flowbca(fi))
  expect_equal(nrow(res$cluster_set), 1)
  expect_equal(nrow(res$unit_set), 4)
  # exactly one core remains
  expect_equal(sum(res$unit_set$core), 1)
  # every unit resolves to the single final cluster
  expect_equal(unique(res$unit_set$clusterid), res$cluster_set$clusterid)
})

test_that("k stopping condition yields k clusters", {
  fi <- make_flow_input(toy_flow_matrix())
  res <- suppressMessages(flowbca(fi, k = 2))
  expect_equal(nrow(res$cluster_set), 2)
  expect_setequal(res$cluster_set$clusterid, c('b', 'd'))
  # merged units point at their cores
  us <- res$unit_set
  expect_equal(us$clusterid[us$sourceunit == 'a'], 'b')
  expect_equal(us$clusterid[us$sourceunit == 'c'], 'd')
})

test_that("q threshold stops before weak merges", {
  # relative flows: a->b is 10/11, c->d is 8/8 = 1, b->d weak
  fi <- make_flow_input(toy_flow_matrix())
  res <- suppressMessages(flowbca(fi, q = 0.95))
  # only merges with relative flow >= 0.95 happen (c -> d)
  expect_true(nrow(res$cluster_set) > 1)
  g <- res$unit_set$g
  expect_true(all(g[!is.na(g)] >= 0.95))
})

test_that("save_k controls whether histories are returned", {
  fi <- make_flow_input(toy_flow_matrix())
  res_no <- suppressMessages(flowbca(fi, k = 2))
  expect_null(res_no$F_matrix_history)
  expect_null(res_no$C_matrix_history)

  res_yes <- suppressMessages(flowbca(fi, k = 2, save_k = TRUE))
  # initial state + one state per merge
  n_merges <- sum(!is.na(res_yes$unit_set$round))
  expect_length(res_yes$F_matrix_history, n_merges + 1)
  expect_length(res_yes$C_matrix_history, n_merges)
  # history names run from (max round + 1) down to min round
  rounds <- res_yes$unit_set$round
  expect_equal(names(res_yes$F_matrix_history),
               as.character((max(rounds, na.rm = TRUE) + 1):min(rounds, na.rm = TRUE)))
})

test_that("results are deterministic", {
  fi <- make_flow_input(toy_flow_matrix())
  r1 <- suppressMessages(flowbca(fi, k = 2))
  r2 <- suppressMessages(flowbca(fi, k = 2))
  expect_identical(r1$unit_set, r2$unit_set)
  expect_identical(r1$cluster_set, r2$cluster_set)
})

test_that("Mata tie-breaking picks the candidate receiving flow from its rivals", {
  # a -> c and b -> c tie at 10 (absolute flows); a -> b = 5 means candidate b
  # receives inflow from the other candidate, so caveat 1 must select r = b.
  m <- matrix(0, 4, 4, dimnames = list(letters[1:4], letters[1:4]))
  m['a', 'c'] <- 10
  m['b', 'c'] <- 10
  m['a', 'b'] <- 5
  m['d', 'd'] <- 100
  res <- suppressMessages(flowbca(make_flow_input(m), opt_f = 3, k = 3))
  first_merge <- res$unit_set[which.max(res$unit_set$round), ]
  expect_equal(first_merge$sourceunit, 'b')
  expect_equal(first_merge$destinationunit, 'c')
})

test_that("non_zero stops once all diagonals are positive", {
  # after c merges into d, both remaining units have internal flow
  m <- matrix(0, 3, 3, dimnames = list(c('c', 'd', 'e'), c('c', 'd', 'e')))
  m['c', 'd'] <- 8
  m['d', 'd'] <- 15
  m['e', 'e'] <- 3
  m['e', 'd'] <- 1
  res <- suppressMessages(flowbca(make_flow_input(m), non_zero = TRUE))
  expect_true(all(diag(res$F_matrix) > 0))
  expect_equal(nrow(res$cluster_set), 2)
})

test_that("cluster_set statistics are internally consistent", {
  fi <- make_flow_input(toy_flow_matrix())
  res <- suppressMessages(flowbca(fi, k = 2))
  cs <- res$cluster_set
  expect_equal(cs$internal_relative,
               ifelse(cs$rowflows == 0, 0, cs$internal / cs$rowflows))
  expect_equal(unique(cs$N), sum(cs$rowflows))
  expect_equal(unique(cs$La), mean(cs$internal_relative))
})
