# A 2x2 grid of unit squares with a 'geometry'-named sf column, deliberately
# different from the 'geom' name used by the bundled KR_SiGun data so that
# hardcoded geometry-column names would fail these tests.
square <- function(x0, y0) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 1, y0), c(x0 + 1, y0 + 1),
                            c(x0, y0 + 1), c(x0, y0))))
}

toy_gis <- function() {
  sf::st_sf(
    id = c('a', 'b', 'c', 'd'),
    geometry = sf::st_sfc(square(0, 0), square(1, 0), square(0, 1), square(1, 1))
  )
}

test_that("remove_holes_sf strips interior rings", {
  skip_if_not_installed("sf")
  outer <- rbind(c(0, 0), c(3, 0), c(3, 3), c(0, 3), c(0, 0))
  hole <- rbind(c(1, 1), c(2, 1), c(2, 2), c(1, 2), c(1, 1))
  with_hole <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(outer, hole))))
  expect_equal(as.numeric(sf::st_area(with_hole)), 8)
  filled <- remove_holes_sf(with_hole)
  expect_equal(as.numeric(sf::st_area(filled)), 9)
})

test_that("flowbca_gis_layer dissolves polygons round by round", {
  skip_if_not_installed("sf")
  unit_set <- data.frame(
    sourceunit = c('a', 'c', 'b', 'd'),
    destinationunit = c('b', 'd', NA, NA),
    round = c(3, 2, NA, NA),
    stringsAsFactors = FALSE
  )
  layers <- suppressMessages(
    flowbca_gis_layer(unit_set, toy_gis(), join_col = c('sourceunit' = 'id'))
  )
  expect_equal(names(layers), c('4', '3', '2'))
  expect_equal(vapply(layers, nrow, integer(1)), c(`4` = 4L, `3` = 3L, `2` = 2L))
  # total area is preserved by dissolving
  expect_equal(sum(as.numeric(sf::st_area(layers[['2']]))), 4)
  # merged unit ids are replaced by their destination
  expect_setequal(layers[['2']]$sourceunit, c('b', 'd'))
})

test_that("flowbca_map works with a 'geometry'-named sf column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("magick")
  unit_set <- data.frame(
    sourceunit = c('a', 'b', 'c', 'd'),
    clusterid = c('b', 'b', 'd', 'd'),
    core = c(0, 1, 0, 1),
    stringsAsFactors = FALSE
  )
  out_png <- tempfile(fileext = '.png')
  grDevices::png(tempfile(fileext = '.png'))  # sink for the preview plot
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- suppressWarnings(
    flowbca_map(unit_set, toy_gis(), join_col = c('sourceunit' = 'id'),
                file_nm = out_png, width = 200)
  )
  expect_true(file.exists(out_png))
  expect_named(res, c('img', 'cluster_gis', 'core_gis'))
  expect_equal(nrow(res$cluster_gis), 2)
  expect_setequal(res$core_gis$id, c('b', 'd'))
})

test_that("flowbca_map accepts hierarchy_cluster level data (h_cl)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("magick")
  level_df <- data.frame(
    sourceunit = c('a', 'b', 'c', 'd'),
    h_cl = c('b', 'b', 'c', 'd'),
    core_N = c(2, 2, 0, 0),
    core = c(0, 1, 0, 0),
    stringsAsFactors = FALSE
  )
  out_png <- tempfile(fileext = '.png')
  grDevices::png(tempfile(fileext = '.png'))
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- suppressWarnings(
    flowbca_map(level_df, toy_gis(), join_col = c('sourceunit' = 'id'),
                file_nm = out_png, width = 200)
  )
  # only the multi-member (core) cluster 'b' is mapped
  expect_equal(res$cluster_gis$clusterid, 'b')
})
