# Toy merge history: a and b merged into r1; r1 and r2 are final clusters.
toy_unit_set <- function() {
  data.frame(
    sourceunit = c('a', 'b', 'r1', 'r2'),
    destinationunit = c('r1', 'r1', NA, NA),
    clusterid = c('r1', 'r1', 'r1', 'r2'),
    round = c(4, 3, NA, NA),
    stringsAsFactors = FALSE
  )
}

test_that("build_hierarchy produces full paths and levels", {
  us <- data.frame(
    sourceunit = c('a', 'b', 'c'),
    destinationunit = c('b', 'c', NA),
    clusterid = 'c',
    stringsAsFactors = FALSE
  )
  h <- build_hierarchy(us)
  expect_equal(h$hierarchy[h$sourceunit == 'a'], 'c/b/a')
  expect_equal(h$h_level[h$sourceunit == 'a'], 3)
  expect_equal(h$hierarchy[h$sourceunit == 'c'], 'c')
  expect_equal(h$h_level[h$sourceunit == 'c'], 1)
  # b and c are parents, a is not
  expect_equal(h$h_parent[match(c('a', 'b', 'c'), h$sourceunit)], c(0, 1, 1))
})

test_that("build_hierarchy warns on circular references", {
  us <- data.frame(
    sourceunit = c('a', 'b'),
    destinationunit = c('b', 'a'),
    clusterid = NA_character_,
    stringsAsFactors = FALSE
  )
  # one warning is raised per unit on the cycle, so capture them all
  warns <- capture_warnings(build_hierarchy(us))
  expect_true(all(grepl("Circular reference", warns)))
  expect_length(warns, 2)
})

test_that("build_cluster_tree builds a nested list keyed by cluster roots", {
  tree <- build_cluster_tree(toy_unit_set())
  expect_setequal(names(tree), c('r1', 'r2'))
  # a parent with children lists itself as a leaf child
  expect_setequal(names(tree$r1), c('a', 'b', 'r1'))
  expect_length(tree$r2, 0)
})

test_that("build_cluster_tree handles a unit_set with no merges", {
  us <- data.frame(sourceunit = c('a', 'b'),
                   destinationunit = NA_character_,
                   stringsAsFactors = FALSE)
  tree <- build_cluster_tree(us)
  expect_setequal(names(tree), c('a', 'b'))
  expect_true(all(lengths(tree) == 0))
})

test_that("build_dendrogram returns a labelled dendrogram", {
  dend <- build_dendrogram(toy_unit_set())
  expect_s3_class(dend, 'dendrogram')
  expect_equal(attr(dend, 'members'), 4)  # a, b, r1 (as leaf), r2
  expect_equal(attr(dend, 'height'), 2)
  labs <- vapply(dend, attr, character(1), 'label')
  expect_setequal(labs, c('r1', 'r2'))
})

test_that("hierarchy_cluster maps top-level clusters to NA parents", {
  hc <- hierarchy_cluster(toy_unit_set())
  expect_named(hc, c('unit_set', 'cluster_info'))

  level1 <- hc$unit_set$hierarchy_1
  expect_setequal(level1$sourceunit, c('a', 'b', 'r1', 'r2'))
  expect_equal(unique(level1$h_cl[level1$sourceunit %in% c('a', 'b')]), 'r1')
  # r1 is a core cluster with 3 members; r2 is a singleton
  expect_equal(level1$core[level1$sourceunit == 'r1'], 1)
  expect_equal(level1$core_N[level1$sourceunit == 'r1'], 3)
  expect_equal(level1$core[level1$sourceunit == 'r2'], 0)

  # top-level clusters have no parent: upper_h_cl must exist and be NA
  info1 <- hc$cluster_info$hierarchy_1
  expect_true('upper_h_cl' %in% names(info1))
  expect_true(all(is.na(info1$upper_h_cl)))
})

test_that("hierarchy_cluster records parent clusters at deeper levels", {
  # three levels: q1, q2 -> p; p, x -> r
  us <- data.frame(
    sourceunit = c('q1', 'q2', 'p', 'x', 'r'),
    destinationunit = c('p', 'p', 'r', 'r', NA),
    clusterid = 'r',
    round = c(5, 4, 3, 2, NA),
    stringsAsFactors = FALSE
  )
  hc <- hierarchy_cluster(us)
  info2 <- hc$cluster_info$hierarchy_2
  expect_equal(unique(info2$upper_h_cl[info2$h_cl == 'p']), 'r')
})
