VAL <- data.frame(read.csv(file = "./RawTest.csv", header = TRUE))
SEQUENCE <- seq_len(37) - 1

test_that("Singleton function accuracy", {
  expect_equal(ddelap(SEQUENCE, 1, 4, 2), VAL$DDELAP_1)
})
test_that("Singleton log", {
  expect_equal(ddelap(SEQUENCE, 5, 3, 8, log = TRUE), log(ddelap(SEQUENCE, 5, 3, 8)))
})
test_that("Singleton NaN", {
  expect_true(is.nan(ddelap(1, 0, 1, 2)))
  expect_true(is.nan(ddelap(1, -2, 1, 2)))
  expect_true(is.nan(ddelap(0, 1, 0, 2)))
  expect_true(is.nan(ddelap(0, 1, -4, 2)))
  expect_true(is.nan(ddelap(0, 1, 4, 0)))
  expect_true(is.nan(ddelap(0, 1, 4, -3)))
})
test_that("Vector function accuracy", {
  expect_equal(ddelap(SEQUENCE, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7)), VAL$DDELAP_Triple)
})
test_that("Vector log", {
  expect_equal(ddelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9), log = TRUE), 
              log(ddelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9))))
})
test_that("Vector NaN", {
  expect_identical(is.nan(ddelap(seq_len(10), 0, 1, 2)), rep(TRUE, 10))
  expect_identical(is.nan(ddelap(seq_len(3), c(0, 1, 2), c(1, 0, 2), c(1, 2, 0))), rep(TRUE, 3))
  expect_identical(is.nan(ddelap(seq_len(3), c(-3, 1, 2), c(1, -5, 2), c(1, 2, -0.1))), rep(TRUE, 3))
})
test_that("Non-integer warning", {
  expect_warning(ddelap(1.1, 1, 2, 3), "Non-integers passed to ddelap. These will have 0 probability.")
  expect_warning(ddelap(c(1, 1.1, 1.2, 3), c(1, 1), 2, 3), "Non-integers passed to ddelap. These will have 0 probability.")
  expect_warning(ddelap(seq(2, 3, .1), c(1, 1), 2, 3), "Non-integers passed to ddelap. These will have 0 probability.")
})

test_that("Non-double parameters converted", {
  expect_equal(ddelap(2L, 1L, 2L, 3L), ddelap(2L, 1, 2, 3))
})