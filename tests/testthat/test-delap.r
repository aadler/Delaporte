context("Mathematical Functionality")

test_that("ddelap pdelap, and qdelap are working", {
  VAL <- data.frame(read.csv(file = "./RawTest.csv", header = TRUE, colClasses = 'numeric'))
  SEQUENCE <- seq_len(37) - 1
  expect_that(ddelap(SEQUENCE, 1, 4, 2), equals(VAL$DDELAP))
  expect_that(pdelap(SEQUENCE, 1, 4, 2), equals(VAL$PDELAP))
  expect_that(qdelap(VAL$PDELAP - .0001, 1, 4, 2), equals(SEQUENCE))
})

test_that("rdelap generates properly", {
  DP1 <- rdelap(1e7, alpha = 10, beta = 2, lambda = 10) 
  DP2 <- rdelap(1e7, alpha = 6, beta = 4, lambda = 6) 
  DP3 <- rdelap(1e7, alpha = 3, beta = 9, lambda = 3)
  expect_that(abs((mean(DP1) / 30 - 1)) < 1e-3, is_true())
  expect_that(abs((mean(DP1) / 30 - 1)) < 1e-3, is_true())
  expect_that(abs((mean(DP1) / 30 - 1)) < 1e-3, is_true())
})

test_that("MoMdelap works", {
  TestData <- c(5, 7, 9, 9, 10, 11, 11, 13, 17, 24)
  MoM <- MoMdelap(TestData)
  expect_that(MoM[[1]], equals(0.769626413810781))
  expect_that(MoM[[2]], equals(4.83611111111112))
  expect_that(MoM[[3]], equals(7.87800114876508))
})

context("Error Trapping")
test_that("MoMdelap traps bad parameters", {
  TestData <- c(3,  2, 12, 11,  1,  7,  1,  4,  0, 4)
  expect_that(MoMdelap(TestData), throws_error("Method of moments not appropriate for this data; results include negative parameters."))
})