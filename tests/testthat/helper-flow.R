# Build a flowbca-ready data.frame from a named square flow matrix.
make_flow_input <- function(mat) {
  data.frame(id = rownames(mat), mat, check.names = FALSE,
             stringsAsFactors = FALSE)
}

# A small 4-unit flow matrix with an unambiguous merge order:
# strongest relative flow is a -> b, then c -> d.
toy_flow_matrix <- function() {
  m <- matrix(0, 4, 4, dimnames = list(letters[1:4], letters[1:4]))
  m['a', 'b'] <- 10
  m['b', 'a'] <- 4
  m['b', 'b'] <- 20   # internal flow so b stays the core
  m['c', 'd'] <- 8
  m['d', 'd'] <- 15
  m['b', 'd'] <- 1
  m
}
