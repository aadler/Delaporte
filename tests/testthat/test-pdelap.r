VAL <- data.frame(read.csv(file = "./RawTest.csv", header = TRUE))
nanWarn <- 'NaNs produced'

test_that("Singleton function accuracy", {
  expect_equal(pdelap(0:36, 2, 1, 5), VAL$PDELAP_2)
})  
test_that("Singleton log.p", {
  expect_equal(pdelap(0:36, 4, 5, 1, log.p = TRUE), log(pdelap(0:36, 4, 5, 1)))
})
test_that("Singleton lower.tail", {
  expect_equal(pdelap(0:100, 8, 10, 6, lower.tail = FALSE), 
               1 - pdelap(0:100, 8, 10, 6))
  expect_equal(pdelap(6, 2.9647, 0.005/2.9647, 0.0057, lower.tail = FALSE), 0)
})
test_that("Singleton lower.tail & log.p", {
  expect_equal(pdelap(0:100, 8, 10, 6, lower.tail = FALSE, log.p = TRUE), 
               log(1 - pdelap(0:100, 8, 10, 6)))
})
test_that("Singleton NaN", {
  expect_warning(pdelap(1, 0, 1, 2), nanWarn)
  expect_warning(pdelap(1, -85, 1, 2), nanWarn)
  expect_warning(pdelap(0, 1, 0, 2), nanWarn)
  expect_warning(pdelap(0, 1, -1e4, 2), nanWarn)
  expect_warning(pdelap(0, 1, 4, 0), nanWarn)
  expect_warning(pdelap(0, 1, 4, -1e-3), nanWarn)
  expect_warning(ddelap(NaN, 1, 4, 6), nanWarn)
  expect_warning(ddelap(NA, 1, 4, 12), nanWarn)
  tst <- suppressWarnings(pdelap(c(NA, 4, NaN), 0.5, 4, 0.2))
  expect_equal(tst[is.finite(tst)], pdelap(4, 0.5, 4, 0.2))
})  
test_that("Vector function accuracy", {
  expect_equal(pdelap(0:36, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7)), VAL$PDELAP_Triple)
})
test_that("Vector log.p", {
  expect_equal(pdelap(0:36, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7), log.p = TRUE),
               log(VAL$PDELAP_Triple))
})
test_that("Vector lower.tail", {
  expect_equal(pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9),
                      lower.tail = FALSE),
               1 - pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9)))
})
test_that("Vector lower.tail & log.p", {
  expect_equal(pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9),
                      lower.tail = FALSE, log.p = TRUE),
               log(1 - pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9))))
})
test_that("Vector NaN", {
  expect_warning(pdelap(1:10, 0, 1, 2), nanWarn)
  expect_warning(pdelap(1:3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0)), nanWarn)
  expect_warning(pdelap(1:3, c(-5e5, 1, 2), c(1, -2, 2), c(1, 2, -8e-4)), nanWarn)
  expect_warning(pdelap(c(NA, 2), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
  expect_warning(pdelap(c(3, NaN), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
  tst <- suppressWarnings(pdelap(c(NA, 0, NaN), c(1, 2, 3), c(4, 1, 2), 
                                 c(2, 5, 7)))
  expect_equal(tst[is.finite(tst)], pdelap(0, 2, 1, 5))
})
test_that("Negative values due to floating point issues are 0", {
  if (R.Version()$arch == "x86_64") {
    expect_identical(pdelap(500, 13.08251, 0.02414521, 0.04421658, FALSE,
                            FALSE), 0)
  } else {
    expect_equal(pdelap(500, 13.08251, 0.02414521, 0.04421658, FALSE, FALSE), 0)
  }
})
test_that("Non-double parameters converted", {
  expect_equal(pdelap(2L, 1L, 2L, 3L), pdelap(2L, 1, 2, 3))
})

test_that("Floating point issues do not lead to CDF > 1", {
  # print(pdelap(1000, 8, 15, 100), digits = 17) used to be 1.0000000000001035
  expect_true(pdelap(1000, 8, 15, 100) <= 1) 
})
test_that("Infinite values", {
  expect_identical(pdelap(Inf, 1L, 2L, 3L), 1)
  expect_identical(pdelap(c(Inf, Inf), c(1L, 2L), 2L, 3L), c(1, 1))
  expect_warning(pdelap(-Inf, 1L, 2L, 3L), nanWarn)
})
