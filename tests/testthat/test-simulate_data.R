test_that("generate_sigma returns correct dimensions", {
  mat <- generate_sigma(n_vars = 5, mode = "low")
  print(mat)
  expect_equal(dim(mat), c(5, 5))
  expect_true(is.matrix(mat))
})
