test_that("flowbca_stat computes internal-flow statistics", {
  m1 <- matrix(c(1, 1,
                 0, 2), 2, 2, byrow = TRUE)  # internal: 1/2 and 2/2
  m2 <- matrix(4, 1, 1)                      # single cluster: all internal
  stats <- flowbca_stat(list(`3` = m1, `2` = m2))

  expect_equal(stats$round, c(3, 2))
  expect_equal(stats$mean, c(mean(c(0.5, 1)), 1))
  expect_equal(stats$min, c(0.5, 1))
  expect_equal(stats$max, c(1, 1))
  # intra ratio = total internal / total flow
  expect_equal(stats$intra_flow_ratio, c(3 / 4, 1))
  expect_equal(stats$intra_flow_ratio + stats$inter_flow_ratio, c(1, 1))
})

test_that("flowbca_stat treats zero-outflow rows as zero internal flow", {
  m <- matrix(0, 2, 2)
  m[1, 1] <- 5
  stats <- flowbca_stat(list(`2` = m))
  expect_equal(stats$min, 0)
  expect_equal(stats$mean, 0.5)
})

test_that("find_intersection interpolates the crossing point", {
  df <- data.frame(
    round = c(10, 9),
    intra_flow_ratio = c(0.4, 0.6),
    inter_flow_ratio = c(0.6, 0.4)
  )
  hit <- flowbcaR:::find_intersection(df)
  expect_equal(hit$round, 9.5)
  expect_equal(hit$ratio, 0.5)
})

test_that("find_intersection handles degenerate inputs", {
  one_row <- data.frame(round = 5, intra_flow_ratio = 0.4, inter_flow_ratio = 0.6)
  expect_null(flowbcaR:::find_intersection(one_row))

  no_cross <- data.frame(
    round = c(10, 9),
    intra_flow_ratio = c(0.1, 0.2),
    inter_flow_ratio = c(0.9, 0.8)
  )
  expect_null(flowbcaR:::find_intersection(no_cross))
})

test_that("flowbca_plot runs without error on real statistics", {
  fi <- make_flow_input(toy_flow_matrix())
  res <- suppressMessages(flowbca(fi, k = 2, save_k = TRUE))
  stats <- flowbca_stat(res$F_matrix_history)

  grDevices::png(tempfile(fileext = '.png'))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(flowbca_plot(stats))
})
