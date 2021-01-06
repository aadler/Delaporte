nanWarn <- 'NaNs produced'
nlErr <- 'invalid arguments'
set.seed(4175L)
DP1 <- rdelap(5e5, alpha = 10, beta = 2, lambda = 10) 
DP2 <- rdelap(5e6, alpha = 2, beta = 14, lambda = 2, exact = FALSE)
DP3 <- rdelap(5e4, alpha = c(5, 5), beta = c(5, 5), lambda = c(5, 5))
DP4 <- rdelap(1e5, alpha = c(5, 5), beta = c(5, 5), lambda = c(5, 5),
              exact = FALSE)

test_that("Singleton exact function accuracy", {
  expect_equal(mean(DP1), 30, tolerance = 5e-4)
})
test_that("Singleton approximate function accuracy", {
  expect_equal(mean(DP2), 30, tolerance = 5e-4)
})
test_that("Singleton NaN", {
  expect_warning(rdelap(1, 0, 1, 2), nanWarn)
  expect_warning(rdelap(1, 1, 0, 2), nanWarn)
  expect_warning(rdelap(1, 1, 4, 0), nanWarn)
  expect_warning(rdelap(1, -3, 1, 2), nanWarn)
  expect_warning(rdelap(1, 1, -5e5, 2), nanWarn)
  expect_warning(rdelap(1, 1, 4, -1e-5), nanWarn)
})
test_that("Singleton size", {
  expect_length(rdelap(8, 4, 1, 2), 8)
  expect_length(rdelap(0, 4, 1, 2), 0)
  expect_error(rdelap(-4, 4, 1, 2), nlErr)
})  
test_that("Vector exact function accuracy", {
  expect_equal(mean(DP3), 30, tolerance = 5e-4)
})  
test_that("Vector approximate function accuracy", {
  expect_equal(mean(DP4), 30, tolerance = 5e-4)
})  
test_that("Vector NaN", {
  expect_warning(rdelap(2, 0, 1, 2), nanWarn)
  expect_warning(rdelap(2, -3, 1, 2), nanWarn)
  expect_warning(rdelap(2, 4, 0, 2), nanWarn)
  expect_warning(rdelap(2, 4, -5e7, 2), nanWarn)
  expect_warning(rdelap(2, 4, 2, 0), nanWarn)
  expect_warning(rdelap(2, 4, 2, -1e-6), nanWarn)
  expect_warning(rdelap(3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0)), nanWarn)
  expect_equal(sum(is.nan(suppressWarnings(rdelap(3, c(0, 1, 2), c(1, 0, 2),
                                                  c(1, 2, 1))))), 2)
})
test_that("Vector size", {
  expect_length(rdelap(8, c(4, 2), c(1, 2, 3, 4), 2), 8)
  expect_length(rdelap(0, c(4, 2), c(1, 2, 3, 4), 2), 0)
  expect_error(rdelap(-1, c(4, 2), c(1, 2, 3, 4), 2), nlErr)
})
test_that("Approximate throws error when nonpositive is passed", {
  expect_warning(rdelap(8, 0, 2, 3, exact = FALSE), nanWarn)
  expect_warning(rdelap(8, 1, -0.8, 3, exact = FALSE), nanWarn)
  expect_warning(rdelap(8, 1, 2, 0, exact = FALSE), nanWarn)
})
test_that("Non-double parameters converted", {
  set.seed(17L)
  INTG <- rdelap(3, 1L, 2L, 3L)
  set.seed(17L)
  DOUBL <- rdelap(3, 1, 2, 3)
  expect_equal(INTG, DOUBL)
})
