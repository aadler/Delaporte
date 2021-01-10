nonIntErr <- 'Non-integers passed to ddelap. These will have 0 probability.'
nanWarn <- 'NaNs produced'
VAL <- data.frame(read.csv(file = "./RawTest.csv", header = TRUE))

test_that("Singleton function accuracy", {
  expect_equal(ddelap(0:36, 1, 4, 2), VAL$DDELAP_1)
})
test_that("alpha < 0.8", {
  expect_equal(ddelap(4, 0.5, 4, 0.2), 0.0547024400602606)
})
test_that("Singleton log", {
  expect_equal(ddelap(0:36, 5, 3, 8, log = TRUE), log(ddelap(0:36, 5, 3, 8)))
})
test_that("Singleton NA", {
  expect_warning(ddelap(1, NA, 2, 3), nanWarn)
  expect_identical(suppressWarnings(ddelap(1:3, 4, NA, 3)), rep(NaN, 3))
})
test_that("Singleton NaN", {
  expect_warning(ddelap(1, 0, 1, 2), nanWarn)
  expect_warning(ddelap(1:10, 0, 1, 2), nanWarn)
  expect_warning(ddelap(1, -2, 1, 2), nanWarn)
  expect_warning(ddelap(0, 1, 0, 2), nanWarn)
  expect_warning(ddelap(0, 1, -4, 2), nanWarn)
  expect_warning(ddelap(0, 1, 4, 0), nanWarn)
  expect_warning(ddelap(0, 1, 4, -3), nanWarn)
  expect_warning(ddelap(0, 1, 4, -3), nanWarn)
  expect_warning(ddelap(NaN, 1, 4, 6), nanWarn)
  expect_warning(ddelap(NA, 1, 4, 12), nanWarn)
  tst <- suppressWarnings(ddelap(c(NA, 4, NaN), 0.5, 4, 0.2))
  expect_equal(tst, c(NaN, ddelap(4, 0.5, 4, 0.2), NaN))
})
test_that("Vector function accuracy", {
  expect_equal(ddelap(0:36, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7)), VAL$DDELAP_Triple)
})
test_that("Vector log", {
  expect_equal(ddelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9), log = TRUE), 
              log(ddelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9))))
})
test_that("Vector NA", {
  expect_identical(suppressWarnings(ddelap(1:3, c(4, 1, 2), c(1, 5, 3), NA)),
                   rep(NaN, 3))
  tst <- suppressWarnings(ddelap(c(4, 4, 4), c(1, 0.5, NA), 4, c(NaN, 0.2, 4)))
  expect_equal(tst, c(NaN, ddelap(4, 0.5, 4, 0.2), NaN))
})
test_that("Vector NaN", {
  expect_warning(ddelap(1:3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0)), nanWarn)
  expect_warning(ddelap(1:3, c(-3, 1, 2), c(1, -5, 2), c(1, 2, -0.1)), nanWarn)
  expect_warning(ddelap(c(NA, 2), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
  expect_warning(ddelap(c(3, NaN), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
  tst <- suppressWarnings(ddelap(c(NA, 0, NaN), c(1, 2, 3), c(4, 1, 2),
                                 c(2, 5, 7)))
  expect_equal(tst, c(NaN, ddelap(0, 2, 1, 5), NaN))
})
test_that("Non-integer warning", {
  expect_warning(ddelap(1.1, 1, 2, 3), nonIntErr)
  expect_warning(ddelap(c(1, 1.1, 1.2, 3), c(1, 1), 2, 3), nonIntErr)
  expect_warning(ddelap(seq(2, 3, .1), c(1, 1), 2, 3), nonIntErr)
})
test_that("Non-double parameters converted", {
  expect_equal(ddelap(2L, 1L, 2L, 3L), ddelap(2L, 1, 2, 3))
})
test_that("Infinite values", {
  expect_identical(ddelap(Inf, 1L, 2L, 3L), 0)
  expect_identical(ddelap(c(Inf, Inf), c(1L, 2L), 2L, 3L), c(0, 0))
  expect_warning(ddelap(-Inf, 1L, 2L, 3L), nanWarn)
})
