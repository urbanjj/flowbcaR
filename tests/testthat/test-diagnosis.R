test_that("flowbca_diagnosis returns merged statistics for both strategies", {
  fi <- make_flow_input(toy_flow_matrix())
  grDevices::png(tempfile(fileext = '.png'))
  on.exit(grDevices::dev.off(), add = TRUE)
  ds <- suppressMessages(flowbca_diagnosis(fi, is_directed = TRUE, data_name = 'toy'))

  expect_named(ds, c('relative', 'absolute'))
  needed <- c('round', 'mean', 'intra_flow_ratio', 'modularity', 'g', 'relative_g')
  expect_true(all(needed %in% names(ds$relative)))
  expect_true(all(needed %in% names(ds$absolute)))
  expect_true(all(ds$relative$relative_g <= 1))
})

test_that("flowbca_diagnosis symmetrizes the input when undirected", {
  # asymmetric input: a->b only; symmetrized it must also flow b->a
  m <- matrix(0, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
  m['a', 'b'] <- 10
  m['c', 'b'] <- 2
  fi <- make_flow_input(m)

  grDevices::png(tempfile(fileext = '.png'))
  on.exit(grDevices::dev.off(), add = TRUE)
  ds <- suppressMessages(flowbca_diagnosis(fi, is_directed = FALSE, data_name = 'toy'))
  # absolute g of the first merge reflects the symmetrized search
  # (undirected absolute: F + t(F), so a<->b carries 20)
  expect_equal(max(ds$absolute$g), 20)
})

test_that("OD_SiGun regression: non_zero run reproduces the documented clusters", {
  skip_on_cran()
  data(OD_SiGun, envir = environment())
  flow_input <- OD_SiGun[, -1]
  colnames(flow_input) <- c('SiGun_NM', flow_input[, 1])
  res <- suppressMessages(flowbca(flow_input, non_zero = TRUE, save_k = FALSE))
  expect_equal(nrow(res$unit_set), 159)
  expect_setequal(
    res$cluster_set$clusterid,
    c('Andong', 'Busan', 'Cheongsong', 'Daegu', 'Daejeon', 'Gwangju', 'Jeonju',
      'Jinju', 'Mokpo', 'Sangju', 'Seoul', 'Suncheon', 'Yeongju')
  )
  # first merge documented in the README
  first <- res$unit_set[which.max(res$unit_set$round), ]
  expect_equal(first$sourceunit, 'Gyeongsan')
  expect_equal(first$clusterid, 'Daegu')
})
