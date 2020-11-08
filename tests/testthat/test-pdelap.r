VAL <- data.frame(read.csv(file = "./RawTest.csv", header = TRUE))
SEQUENCE <- seq_len(37) - 1

test_that("Singleton function accuracy", {
  expect_equal(pdelap(SEQUENCE, 2, 1, 5), VAL$PDELAP_2)
})  
test_that("Singleton log.p", {
  expect_equal(pdelap(SEQUENCE, 4, 5, 1, log.p = TRUE), log(pdelap(SEQUENCE, 4, 5, 1)))
})
test_that("Singleton lower.tail", {
  expect_equal(pdelap(seq_len(101) - 1, 8, 10, 6, lower.tail = FALSE), 1 - pdelap(seq_len(101) - 1, 8, 10, 6))
  expect_equal(pdelap(6, 2.9647, 0.005/2.9647, 0.0057, lower.tail = FALSE), 0)
})
test_that("Singleton lower.tail & log.p", {
  expect_equal(pdelap(seq_len(101) - 1, 8, 10, 6, lower.tail = FALSE, log.p = TRUE),
               log(1 - pdelap(seq_len(101) - 1, 8, 10, 6)))
})
test_that("Singleton NaN", {
  expect_true(is.nan(pdelap(1, 0, 1, 2)))
  expect_true(is.nan(pdelap(1, -85, 1, 2)))
  expect_true(is.nan(pdelap(0, 1, 0, 2)))
  expect_true(is.nan(pdelap(0, 1, -1e4, 2)))
  expect_true(is.nan(pdelap(0, 1, 4, 0)))
  expect_true(is.nan(pdelap(0, 1, 4, -1e-3)))
})  
test_that("Vector function accuracy", {
  expect_equal(pdelap(SEQUENCE, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7)), VAL$PDELAP_Triple)
})
test_that("Vector log.p", {
  expect_equal(pdelap(SEQUENCE, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7), log.p = TRUE), log(VAL$PDELAP_Triple))
})
test_that("Vector lower.tail", {
  expect_equal(pdelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9), lower.tail = FALSE),
               1 - pdelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9)))
})
test_that("Vector lower.tail & log.p", {
  expect_equal(pdelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9), lower.tail = FALSE, log.p = TRUE),
               log(1 - pdelap(seq_len(101) - 1, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9))))
})
test_that("Vector NaN", {
  expect_identical(is.nan(pdelap(seq_len(10), 0, 1, 2)), rep(TRUE, 10))
  expect_identical(is.nan(pdelap(seq_len(3), c(0, 1, 2), c(1, 0, 2), c(1, 2, 0))), rep(TRUE, 3))
  expect_identical(is.nan(pdelap(seq_len(3), c(-5e5, 1, 2), c(1, -2, 2), c(1, 2, -8e-4))), rep(TRUE, 3))
})

test_that("Negative values due to floating point issues are 0", {
  if (R.Version()$arch == "x86_64") {
    expect_identical(pdelap(500, 13.08251, 0.02414521, 0.04421658, FALSE, FALSE), 0)
  } else {
    expect_equal(pdelap(500, 13.08251, 0.02414521, 0.04421658, FALSE, FALSE), 0)
  }
})

test_that("Non-double parameters converted", {
  expect_equal(pdelap(2L, 1L, 2L, 3L), pdelap(2L, 1, 2, 3))
})

test_that("Floating point issues do not lead to CDF > 1", {
  expect_true(pdelap(1000, 8, 15, 100) <= 1) # print(pdelap(1000, 8, 15, 100), digits = 17) used to be 1.0000000000001035
})
