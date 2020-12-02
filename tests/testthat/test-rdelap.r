zeroErr <- 'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.'
negLenErr <- 'negative length vectors are not allowed'
set.seed(4175L)
DP1 <- rdelap(1e6, alpha = 10, beta = 2, lambda = 10) 
DP2 <- rdelap(3e6, alpha = 2, beta = 14, lambda = 2, exact = FALSE)
DP3 <- rdelap(1.75e5, alpha = c(5, 5), beta = c(5, 5), lambda = c(5, 5))
DP4 <- rdelap(1e6, alpha = c(5, 5), beta = c(5, 5), lambda = c(5, 5),
              exact = FALSE)

test_that("Singleton exact function accuracy", {
  expect_true(abs((mean(DP1) / 30 - 1)) < 5e-4)
})
test_that("Singleton approximate function accuracy", {
  expect_true(abs((mean(DP2) / 30 - 1)) < 5e-4)
})
test_that("Singleton NaN", {
  expect_true(is.nan(rdelap(1, 0, 1, 2)))
  expect_true(is.nan(rdelap(1, 1, 0, 2)))
  expect_true(is.nan(rdelap(1, 1, 4, 0)))
  expect_true(is.nan(rdelap(1, -3, 1, 2)))
  expect_true(is.nan(rdelap(1, 1, -5e5, 2)))
  expect_true(is.nan(rdelap(1, 1, 4, -1e-5)))
})
test_that("Singleton size", {
  expect_length(rdelap(8, 4, 1, 2), 8)
  expect_length(rdelap(0, 4, 1, 2), 0)
  expect_error(rdelap(-4, 4, 1, 2), negLenErr)
})  
test_that("Vector exact function accuracy", {
  expect_true(abs((mean(DP3) / 30 - 1)) < 5e-4)
})  
test_that("Vector approximate function accuracy", {
  expect_true(abs((mean(DP4) / 30 - 1)) < 5e-4)
})  
test_that("Vector NaN", {
  expect_identical(is.nan(rdelap(2, 0, 1, 2)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, -3, 1, 2)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, 4, 0, 2)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, 4, -5e7, 2)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, 4, 2, 0)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, 4, 2, -1e-6)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0))),
                   rep(TRUE, 3))
  expect_equal(sum(is.nan(rdelap(3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 1)))), 2)
})
test_that("Vector size", {
  expect_length(rdelap(8, c(4, 2), c(1, 2, 3, 4), 2), 8)
  expect_length(rdelap(0, c(4, 2), c(1, 2, 3, 4), 2), 0)
  expect_error(rdelap(-1, c(4, 2), c(1, 2, 3, 4), 2), negLenErr)
})
test_that("Approximate throws error when nonpositive is passed", {
  expect_error(rdelap(8, 0, 2, 3, exact = FALSE), zeroErr)
  expect_error(rdelap(8, 1, -0.8, 3, exact = FALSE), zeroErr)
  expect_error(rdelap(8, 1, 2, 0, exact = FALSE), zeroErr)
})

test_that("Non-double parameters converted", {
  set.seed(17L)
  INTG <- rdelap(3, 1L, 2L, 3L)
  set.seed(17L)
  DOUBL <- rdelap(3, 1, 2, 3)
  expect_equal(INTG, DOUBL)
})