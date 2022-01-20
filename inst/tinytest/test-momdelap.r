TestData <- c(5, 7, 9, 9, 10, 11, 11, 13, 17, 24)
MTD <- mean(TestData)
VTD <- var(TestData)
VmMTD <- VTD - MTD
SK1 <- 1.1944183061921252 # Three values from skewness in e1071 package
SK2 <- 1.4164058724628172
SK3 <- 1.0198122281732285

# Function accuracy - type 1 explicit
MoM <- MoMdelap(TestData, type = 1L)
P2 <- 0.5 * (SK1 * (VTD ^ 1.5) - MTD - 3 * VmMTD) / VmMTD
P1 <- VmMTD / (P2 ^ 2)
P3 <- MTD - P1 * P2
expect_equal(MoM[[1]], P1)
expect_equal(MoM[[2]], P2)
expect_equal(MoM[[3]], P3)

# Function accuracy - type 2 implicit
MoM <- MoMdelap(TestData)
P2 <- 0.5 * (SK2 * (VTD ^ 1.5) - MTD - 3 * VmMTD) / VmMTD
P1 <- VmMTD / (P2 ^ 2)
P3 <- MTD - P1 * P2
expect_equal(MoM[[1]], P1)
expect_equal(MoM[[2]], P2)
expect_equal(MoM[[3]], P3)

# Function accuracy - type 2 explicit
MoM <- MoMdelap(TestData, type = 2)
expect_equal(MoM[[1]], 0.88342721893491116)
expect_equal(MoM[[2]], 4.51388888888888928)
expect_equal(MoM[[3]], 7.61230769230769155)

# Function accuracy - type 3 explicit
MoM <- MoMdelap(TestData, type = 3)
P2 <- 0.5 * (SK3 * (VTD ^ 1.5) - MTD - 3 * VmMTD) / VmMTD
P1 <- VmMTD / (P2 ^ 2)
P3 <- MTD - P1 * P2
expect_equal(MoM[[1]], P1)
expect_equal(MoM[[2]], P2)
expect_equal(MoM[[3]], P3)

# MoMdelap traps bad types
expect_error(MoMdelap(TestData, type = 4),
             'Skew type must be one of 1, 2, or 3.')

# MoMdelap traps bad parameters
TestData <- c(3,  2, 12, 11,  1,  7,  1,  4,  0, 4)
expect_error(MoMdelap(TestData),
             'Method of moments not appropriate for this data; results include non-positive parameters.')

# Non-double vector converted
expect_equal(MoMdelap(c(30L, 32L, 39L, 50L), type = 2L),
             MoMdelap(c(30, 32, 39, 50), type = 2L))
