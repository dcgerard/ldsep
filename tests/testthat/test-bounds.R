test_that("find_bounds() works", {

  expect_equal(c(1, 2), c(find_bounds_cpp(1.5, 1:10)))
  expect_equal(c(2, 3), c(find_bounds_cpp(2.5, 1:10)))
  expect_equal(c(3, 4), c(find_bounds_cpp(3.5, 1:10)))
  expect_equal(c(4, 5), c(find_bounds_cpp(4.5, 1:10)))
  expect_equal(c(5, 6), c(find_bounds_cpp(5.5, 1:10)))
  expect_equal(c(6, 7), c(find_bounds_cpp(6.5, 1:10)))
  expect_equal(c(7, 8), c(find_bounds_cpp(7.5, 1:10)))
  expect_equal(c(8, 9), c(find_bounds_cpp(8.5, 1:10)))
  expect_equal(c(9, 10), c(find_bounds_cpp(9.5, 1:10)))

  expect_equal(c(1, 2), c(find_bounds_cpp(1.5, 1:11)))
  expect_equal(c(2, 3), c(find_bounds_cpp(2.5, 1:11)))
  expect_equal(c(3, 4), c(find_bounds_cpp(3.5, 1:11)))
  expect_equal(c(4, 5), c(find_bounds_cpp(4.5, 1:11)))
  expect_equal(c(5, 6), c(find_bounds_cpp(5.5, 1:11)))
  expect_equal(c(6, 7), c(find_bounds_cpp(6.5, 1:11)))
  expect_equal(c(7, 8), c(find_bounds_cpp(7.5, 1:11)))
  expect_equal(c(8, 9), c(find_bounds_cpp(8.5, 1:11)))
  expect_equal(c(9, 10), c(find_bounds_cpp(9.5, 1:11)))
  expect_equal(c(10, 11), c(find_bounds_cpp(10.5, 1:11)))

  ## find bounds not worth it for 200
  # bench::mark(
  #   find_bounds(33.4, 1:200),
  #   c(find_bounds_cpp(33.4, 1:200))
  # )

})
