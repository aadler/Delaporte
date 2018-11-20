VAL <- data.frame(read.csv(file = "./RawTest.csv", header = TRUE))
SEQUENCE <- seq_len(37) - 1

context("Testing ddelap")
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
  
context("Testing pdelap")
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
  expect_identical(pdelap(500, 13.08251, 0.02414521, 0.04421658, FALSE, FALSE), 0)
})

context("Testing qdelap")
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
test_that("Singleton approx function accuracy", {
  expect_equal(qdelap(.4, 1, 4, 2, exact = FALSE), 4)
})
test_that("Singleton approx lower.tail & log.p", {
  expect_equal(qdelap(-0.7, 4, 6, 3, lower.tail = TRUE, log.p = TRUE, exact = FALSE), 25)
})
test_that("Vector exact function accuracy", {
  expect_equal(qdelap(c(.4, .07), c(1, 2), c(4, 1), c(2, 5)), c(4, 3))
})
test_that("Vector exact lower.tail", {
  expect_equal(qdelap(c(0.49, 0.131), c(4, 1), c(6, 9), c(3, 12), lower.tail = FALSE), c(25, 31))
})
test_that("Vector exact log.p", {
  expect_equal(qdelap(c(-0.9162907318741550, -2.6592600369327779), c(1, 2), c(4, 1), c(2, 5), log.p = TRUE), c(4, 3))
})  
test_that("Vector exact lower.tail & log.p", {
  expect_equal(qdelap(c(-0.69895775020315487, -1.98413706125967337), c(4, 1), c(6, 9), c(3, 12),
                      lower.tail = FALSE, log.p = TRUE), c(25, 31))
})
test_that("Vector Nan", {
  expect_identical(is.nan(qdelap(seq_len(2)/10, 0, 1, 2)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, -1, 1, 2)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, 1, 0, 2)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, 1, -8, 2)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, 3, 1, 0)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(2)/10, 3, 1, -4e-5)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(c(-1, -5), 3, 1, 6)), rep(TRUE, 2))
  expect_identical(is.nan(qdelap(seq_len(3)/10, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0))), rep(TRUE, 3))
})
test_that("Vector Inf", {
  expect_identical(is.infinite(qdelap(seq_len(2), 3, 1, 2)), rep(TRUE, 2))
  expect_identical(is.infinite(qdelap(seq_len(3), c(2, 1, 2), c(1, 6, 2), c(1, 2, 0.4))), rep(TRUE, 3))
})
test_that("Approximate throws error when 0 is passed", {
  expect_error(qdelap(0.1, 0, 2, 3, exact = FALSE), 'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
  expect_error(qdelap(0.1, 1, 0, 3, exact = FALSE), 'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
  expect_error(qdelap(0.1, 1, 2, 0, exact = FALSE), 'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
})
test_that("Approximate throws error when parameter vectors are passed", {
  expect_error(qdelap(c(.4, .07), c(1, 2), c(4, 1), c(2, 5), exact = FALSE), 'Quantile approximation relies on pooling and thus is not accurate when passing vector-valued parameters. Please use exact version.')
})

context("Testing rdelap")
DP1 <- rdelap(1e6, alpha = 10, beta = 2, lambda = 10) 
DP2 <- rdelap(1e6, alpha = 2, beta = 14, lambda = 2, exact = FALSE)
DP3 <- rdelap(1e4, alpha = c(5, 5), beta = c(5, 5), lambda = c(5, 5))
DP4 <- rdelap(1e6, alpha = c(5, 5), beta = c(5, 5), lambda = c(5, 5), exact = FALSE)

test_that("Singleton exact function accuracy", {
  expect_that(abs((mean(DP1) / 30 - 1)) < 5e-3, is_true())
})
test_that("Singleton approximate function accuracy", {
  expect_that(abs((mean(DP2) / 30 - 1)) < 5e-3, is_true())
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
  expect_error(rdelap(-4, 4, 1, 2), "negative length vectors are not allowed")
})  
test_that("Vector exact function accuracy", {
  expect_true(abs((mean(DP3) / 30 - 1)) < 5e-2)
})  
test_that("Vector approximate function accuracy", {
  expect_true(abs((mean(DP4) / 30 - 1)) < 5e-3)
})  
test_that("Vector NaN", {
  expect_identical(is.nan(rdelap(2, 0, 1, 2)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, -3, 1, 2)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, 4, 0, 2)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, 4, -5e7, 2)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, 4, 2, 0)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(2, 4, 2, -1e-6)), rep(TRUE, 2))
  expect_identical(is.nan(rdelap(3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0))), rep(TRUE, 3))
  expect_equal(sum(is.nan(rdelap(3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 1)))), 2)
})
test_that("Vector size", {
  expect_length(rdelap(8, c(4, 2), c(1, 2, 3, 4), 2), 8)
  expect_length(rdelap(0, c(4, 2), c(1, 2, 3, 4), 2), 0)
  expect_error(rdelap(-1, c(4, 2), c(1, 2, 3, 4), 2), "negative length vectors are not allowed")
})
test_that("Approximate throws error when 0 is passed", {
  expect_error(rdelap(8, 0, 2, 3, exact = FALSE), 'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
  expect_error(rdelap(8, 1, 0, 3, exact = FALSE), 'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
  expect_error(rdelap(8, 1, 2, 0, exact = FALSE), 'Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results')
})

context("Testing MoMdelap")
TestData <- c(5, 7, 9, 9, 10, 11, 11, 13, 17, 24)
MTD <- mean(TestData)
VTD <- var(TestData)
VmMTD <- VTD - MTD
SK1 <- 1.1944183061921252 # Three values from skewness in e1071 package
SK2 <- 1.4164058724628172
SK3 <- 1.0198122281732285

test_that("Function accuracy - type 1 explicit", {
  MoM <- MoMdelap(TestData, type = 1L)
  P2 <- 0.5 * (SK1 * (VTD ^ 1.5) - MTD - 3 * VmMTD) / VmMTD
  P1 <- VmMTD / (P2 ^ 2)
  P3 <- MTD - P1 * P2
  expect_equal(MoM[[1]], P1)
  expect_equal(MoM[[2]], P2)
  expect_equal(MoM[[3]], P3)
})
test_that("Function accuracy - type 2 implicit", {
  MoM <- MoMdelap(TestData)
  P2 <- 0.5 * (SK2 * (VTD ^ 1.5) - MTD - 3 * VmMTD) / VmMTD
  P1 <- VmMTD / (P2 ^ 2)
  P3 <- MTD - P1 * P2
  expect_equal(MoM[[1]], P1)
  expect_equal(MoM[[2]], P2)
  expect_equal(MoM[[3]], P3)
})
test_that("Function accuracy - type 2 explicit", {
  MoM <- MoMdelap(TestData, type = 2)
  expect_equal(MoM[[1]], 0.88342721893491116)
  expect_equal(MoM[[2]], 4.51388888888888928)
  expect_equal(MoM[[3]], 7.61230769230769155)
})
test_that("Function accuracy - type 3 explicit", {
  MoM <- MoMdelap(TestData, type = 3)
  P2 <- 0.5 * (SK3 * (VTD ^ 1.5) - MTD - 3 * VmMTD) / VmMTD
  P1 <- VmMTD / (P2 ^ 2)
  P3 <- MTD - P1 * P2
  expect_equal(MoM[[1]], P1)
  expect_equal(MoM[[2]], P2)
  expect_equal(MoM[[3]], P3)
})
test_that("MoMdelap traps bad types", {
  expect_error(MoMdelap(TestData, type = 4),
               'Skew type must be one of 1, 2, or 3.')
})
test_that("MoMdelap traps bad parameters", {
  TestData <- c(3,  2, 12, 11,  1,  7,  1,  4,  0, 4)
  expect_error(MoMdelap(TestData), 'Method of moments not appropriate for this data; results include non-positive parameters.')
})