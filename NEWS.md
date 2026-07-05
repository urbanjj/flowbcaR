# flowbcaR 1.1.0

## Breaking changes

* flowbcaR now requires R >= 4.6, with version floors on Imports
  (magick >= 2.8.0, sf >= 1.0-0).
* `flowbca()`: the Mata tie-breaking rules (caveats 1-4 of the original
  Stata implementation) are now actually applied. Previously an internal
  `-Inf` diagonal made the tie-breaking code unreachable, so ties always
  fell back to the first candidate. Results may change for data with tied
  maximum flows (typically absolute-flow or symmetric inputs); they now
  match the original `flowbca.ado` behaviour.
* `flowbca()`: a `smaller` value outside `[0, 1]` is now an error instead
  of a repeated warning.
* `flowbca_diagnosis()`: when `is_directed = FALSE`, the input flow matrix
  is now symmetrized (`F + t(F)`) as documented, and the result is returned
  invisibly.
* `calculate_modularity(is_directed = FALSE)` returned exactly twice the
  correct value because the matrix was symmetrized twice. Undirected
  modularity now matches the standard Newman formulation.

## New features

* `build_cluster_tree()` is now exported: it converts the merge history
  into a named nested list (tree), the building block shared by
  `build_dendrogram()` and `hierarchy_cluster()`.
* `flowbca_map()` and `flowbca_gis_layer()` now work with any geometry
  column name (previously hardcoded to `geom`) and with custom join
  column names.
* `hierarchy_cluster()`: top-level clusters (which have no parent) are now
  reported with `upper_h_cl = NA`. Previously the column could be silently
  dropped or misaligned.
* `flowbca()` no longer accumulates the flow-matrix history when
  `save_k = FALSE`, reducing memory use on large inputs.

## Bug fixes

* Edge cases no longer error: `flowbca_stat()`/`flowbca_plot()` with a
  single round, `flowbca_gis_layer()` and the modularity helpers with a
  merge-free `unit_set`, and `flowbca_ani()` with a single layer (which now
  gives an informative error).
* `flowbca_map()` no longer fails when more clusters than named R colors
  exist, and its id-column handling follows `join_col` correctly.
* Parameter documentation for `flowbca()` now reflects the actual defaults
  (all stopping conditions default to `NULL`; exactly one may be specified;
  `k = 1` is used when none is given).

## Documentation and infrastructure

* README and vignette fully rewritten sections for the Hierarchy
  (`build_cluster_tree()`, `build_dendrogram()`, `hierarchy_cluster()`) and
  Visualization (`flowbca_map()`) tools, with a redrawn workflow diagram.
* New testthat suite (105 tests) covering the clustering algorithm,
  modularity reference values, hierarchy tools, spatial layers, and an
  `OD_SiGun` regression snapshot. `R CMD check`: Status OK.

# flowbcaR 1.0.0

* Initial release: R implementation of the flow-based clustering algorithm
  of Meekes and Hassink (2018) with diagnosis (`flowbca_stat()`,
  `flowbca_plot()`, `flowbca_modularity()`, `flowbca_diagnosis()`),
  hierarchy (`build_hierarchy()`, `build_dendrogram()`,
  `hierarchy_cluster()`), and visualization (`flowbca_gis_layer()`,
  `flowbca_map()`, `flowbca_ani()`) tools, plus the `OD_SiGun` and
  `KR_SiGun` example datasets.
