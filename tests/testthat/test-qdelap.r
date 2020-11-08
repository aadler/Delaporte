test_that("Singleton exact function accuracy", {
  expect_equal(qdelap(.4, 1, 4, 2), 4)
})
test_that("Singleton exact lower.tail", {
  expect_equal(qdelap(0.49, 4, 6, 3, lower.tail = TRUE), 25)
})
test_that("Singleton exact log.p", {
  expect_equal(qdelap(-0.255, 20, 15, 50, log.p = TRUE), 400)
})
test_that("Singleton exact lower.tail & log.p", {
  expect_equal(qdelap(-0.7, 4, 6, 3, lower.tail = TRUE, log.p = TRUE), 25)
})
test_that("Singleton exact NaN", {
  expect_true(is.nan(qdelap(.05, 0, 1, 2)))
  expect_true(is.nan(qdelap(.05, -2, 1, 2)))
  expect_true(is.nan(qdelap(0.1, 1, 0, 2)))
  expect_true(is.nan(qdelap(0.1, 1, -4e5, 2)))
  expect_true(is.nan(qdelap(1, 1, 4, 0)))
  expect_true(is.nan(qdelap(1, 1, 4, -9e-4)))
  expect_true(is.nan(qdelap(-1, 2, 3, 4)))
})
test_that("Singleton exact Inf", {
  expect_true(is.infinite(qdelap(1, 3, 1, 2)))
  expect_true(is.infinite(qdelap(5, 1, 2, 3)))
})
test_that("Singleton approx function accuracy lower tail", {
  expect_equal(qdelap(.4, 1, 4, 2, exact = FALSE), 4)
})
test_that("Singleton approx function accuracy upper tail", {
  expect_equal(qdelap(.4, 1, 4, 2, exact = FALSE, lower.tail = FALSE), 6)
})
test_that("Singleton approx lower.tail & log.p", {
  expect_equal(qdelap(-0.7, 4, 6, 3, lower.tail = TRUE, log.p = TRUE,
                      exact = FALSE), 25)
})

test_that("Vector exact function accuracy", {
  expect_equal(qdelap(c(.4, .07), c(1, 2), c(4, 1), c(2, 5)), c(4, 3))
})
test_that("Vector exact lower.tail", {
  expect_equal(qdelap(c(0.49, 0.131), c(4, 1), c(6, 9), c(3, 12),
                      lower.tail = FALSE), c(25, 31))
})
test_that("Vector exact log.p", {
  expect_equal(qdelap(c(-0.9162907318741550, -2.6592600369327779), c(1, 2),
                      c(4, 1), c(2, 5), log.p = TRUE), c(4, 3))
})  
test_that("Vector exact lower.tail & log.p", {
  expect_equal(qdelap(c(-0.69895775020315487, -1.98413706125967337), c(4, 1),
                      c(6, 9), c(3, 12), lower.tail = FALSE, log.p = TRUE),
               c(25, 31))
})
test_that("Vector Nan", {
  expect_identical(is.nan(qdelap(seq_len(2)/10, 0, 1, 2)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, -1, 1, 2)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, 1, 0, 2)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, 1, -8, 2)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, 3, 1, 0)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, 3, 1, -4e-5)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(c(-1, -5), 3, 1, 6)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(3)/10, c(0, 1, 2), c(1, 0, 2),
                                 c(1, 2, 0))), rep(TRUE, 3))
})
test_that("Vector Inf", {
  expect_identical(is.infinite(qdelap(seq_len(2), 3, 1, 2)), rep(TRUE, 2))
  expect_identical(is.infinite(qdelap(seq_len(3), c(2, 1, 2),
                                      c(1, 6, 2), c(1, 2, 0.4))), rep(TRUE, 3))
})
test_that("Approximate throws error when 0 is passed", {
  expect_error(qdelap(0.1, 0, 2, 3, exact = FALSE),
               'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
  expect_error(qdelap(0.1, 1, 0, 3, exact = FALSE),
               'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
  expect_error(qdelap(0.1, 1, 2, 0, exact = FALSE),
               'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
})
test_that("Approximate throws error when parameter vectors are passed", {
  expect_error(qdelap(c(.4, .07), c(1, 2), c(4, 1), c(2, 5), exact = FALSE),
               'Quantile approximation relies on pooling and thus is not accurate when passing vector-valued parameters. Please use exact version.')
})

test_that("Non-double parameters converted", {
  expect_equal(qdelap(0.25, 1L, 2L, 3L), qdelap(0.25, 1, 2, 3))
})
