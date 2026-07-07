# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# For CRAN
oldThreads <- getDelapThreads()
setDelapThreads(2L)

tol <- 1e-12
nanWarn <- "NaNs produced"
inpWarn <- "Quantile approximation relies on pooling"

# Singleton exact function accuracy
expect_equal(qdelap(0.4, 1, 4, 2), 4, tolerance = tol)
testV <- c(3, 7, 23)
expect_equal(qdelap(pdelap(testV, 5L, 2L, 3L), 5L, 2L, 3L), testV,
             tolerance = tol)

# Singleton exact lower.tail
expect_equal(qdelap(0.49, 4, 6, 3, lower.tail = TRUE), 25, tolerance = tol)

# Singleton exact log.p
expect_equal(qdelap(-0.255, 20, 15, 50, log.p = TRUE), 400, tolerance = tol)

# Singleton exact lower.tail & log.p
expect_equal(qdelap(-0.7, 4, 6, 3, lower.tail = TRUE, log.p = TRUE), 25,
             tolerance = tol)

# Singleton exact bad parameters
expect_warning(qdelap(0.05, 0, 1, 2), nanWarn)
expect_warning(qdelap(0.05, -2, 1, 2), nanWarn)
expect_warning(qdelap(0.1, 1, 0, 2), nanWarn)
expect_warning(qdelap(0.1, 1, -4e5, 2), nanWarn)
expect_warning(qdelap(1, 1, 4, 0), nanWarn)
expect_warning(qdelap(1, 1, 4, -9e-4), nanWarn)

# Singleton exact bad inputs
expect_warning(qdelap(-1, 2, 3, 4), nanWarn)
expect_warning(qdelap(NaN, 2, 3, 4), nanWarn)
expect_warning(qdelap(c(0.3, NA), 2, 3, 4), nanWarn)

# Singleton approx function accuracy lower tail
expect_equal(qdelap(0.4, 1, 4, 2, exact = FALSE), 4, tolerance = tol)

# Singleton approx function accuracy upper tail
expect_equal(qdelap(0.4, 1, 4, 2, exact = FALSE, lower.tail = FALSE), 6,
             tolerance = tol)

# Singleton approx lower.tail & log.p
expect_equal(qdelap(-0.7, 4, 6, 3, lower.tail = TRUE, log.p = TRUE,
                    exact = FALSE), 25, tolerance = tol)

# Singleton approx bad parameters
expect_warning(qdelap(0.05, 0, 1, 2, exact = FALSE), nanWarn)
expect_warning(qdelap(0.05, -2, 1, 2, exact = FALSE), nanWarn)
expect_warning(qdelap(0.1, 1, 0, 2, exact = FALSE), nanWarn)
expect_warning(qdelap(0.1, 1, -4e5, 2, exact = FALSE), nanWarn)
expect_warning(qdelap(1, 1, 4, 0, exact = FALSE), nanWarn)
expect_warning(qdelap(1, 1, 4, -9e-4, exact = FALSE), nanWarn)

# Singleton approx bad inputs
expect_warning(qdelap(-1, 2, 3, 4, exact = FALSE), nanWarn)
expect_warning(qdelap(c(0.2, NaN), 2, 3, 4, exact = FALSE), inpWarn)
expect_warning(qdelap(c(0.3, NA), 2, 3, 4, exact = FALSE), inpWarn)

# Vector exact function accuracy
expect_equal(qdelap(c(0.4, 0.07), c(1, 2), c(4, 1), c(2, 5)), c(4, 3),
             tolerance = tol)

# Vector exact lower.tail
expect_equal(qdelap(c(0.49, 0.131), c(4, 1), c(6, 9), c(3, 12),
                    lower.tail = FALSE), c(25, 31), tolerance = tol)

# Vector exact log.p
expect_equal(qdelap(c(-0.9162907318741550, -2.6592600369327779), c(1, 2),
                    c(4, 1), c(2, 5), log.p = TRUE), c(4, 3),
             tolerance = tol)

# Vector exact lower.tail & log.p
expect_equal(qdelap(c(-0.69895775020315487, -1.98413706125967337), c(4, 1),
                    c(6, 9), c(3, 12), lower.tail = FALSE, log.p = TRUE),
             c(25, 31), tolerance = tol)

# Vector exact bad parameters
t2 <- 1:2 / 10
t3 <- 1:3 / 10
expect_warning(qdelap(t2, c(0, 1), 1, 2), nanWarn)
expect_warning(qdelap(t2, c(1, -1), 1, 2), nanWarn)
expect_warning(qdelap(t2, 1, c(2, 0), 2), nanWarn)
expect_warning(qdelap(t2, 1, c(-8, 3), 2), nanWarn)
expect_warning(qdelap(t2, 3, 1, c(2, 0)), nanWarn)
expect_warning(qdelap(t2, 3, 1, c(-4e-5, 12)), nanWarn)
expect_warning(qdelap(t3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0)), nanWarn)
expect_warning(qdelap(t3 / 10, c(6, 1, 2), c(1, 4, 2), c(1, 2, -1)), nanWarn)

# Vector exact bad inputs
expect_warning(qdelap(c(-1, 3), c(1, 3), 1, 6), nanWarn)
expect_warning(qdelap(c(NA, 4), c(1, 3), 1, 6), nanWarn)
expect_warning(qdelap(c(5, NaN), c(1, 3), 1, 6), nanWarn)

# Vector approx bad parameters
t2 <- 1:2 / 10
t3 <- 1:3 / 10
expect_warning(qdelap(t2, c(0, 1), 1, 2, exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(t2, c(0, 1), 1, 2, exact = FALSE)),
                 suppressWarnings(qdelap(t2, c(0, 1), 1, 2)))
expect_warning(qdelap(t2, c(1, -1), 1, 2, exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(t2, c(1, -1), 1, 2, exact = FALSE)),
                 suppressWarnings(qdelap(t2, c(1, -1), 1, 2)))
expect_warning(qdelap(t2, 1, c(2, 0), 2, exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(t2, 1, c(2, 0), 2, exact = FALSE)),
                 suppressWarnings(qdelap(t2, 1, c(2, 0), 2)))
expect_warning(qdelap(t2, 1, c(-8, 3), 2, exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(t2, 1, c(-8, 3), 2, exact = FALSE)),
                 suppressWarnings(qdelap(t2, 1, c(-8, 3), 2)))
expect_warning(qdelap(t2, 3, 1, c(2, 0), exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(t2, 3, 1, c(2, 0), exact = FALSE)),
                 suppressWarnings(qdelap(t2, 3, 1, c(2, 0))))
expect_warning(qdelap(t2, 3, 1, c(-4e-5, 12), exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(t2, 3, 1, c(-4e-5, 12),
                                         exact = FALSE)),
                 suppressWarnings(qdelap(t2, 3, 1, c(-4e-5, 12))))
expect_warning(qdelap(t3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0), exact = FALSE),
               inpWarn)
expect_identical(suppressWarnings(qdelap(t3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0),
                                         exact = FALSE)),
                 suppressWarnings(qdelap(t3, c(0, 1, 2), c(1, 0, 2),
                                         c(1, 2, 0))))
expect_warning(qdelap(t3 / 10, c(6, 1, 2), c(1, 4, 2), c(1, 2, -1),
                      exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(t3 / 10, c(6, 1, 2), c(1, 4, 2),
                                         c(1, 2, -1), exact = FALSE)),
                 suppressWarnings(qdelap(t3 / 10, c(6, 1, 2), c(1, 4, 2),
                                         c(1, 2, -1))))

# Vector exact bad inputs
expect_warning(qdelap(c(-1, 3), c(1, 3), 1, 6, exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(c(-1, 3), c(1, 3), 1, 6,
                                         exact = FALSE)),
                 suppressWarnings(qdelap(c(-1, 3), c(1, 3), 1, 6)))
expect_warning(qdelap(c(NA, 4), c(1, 3), 1, 6, exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(c(NA, 4), c(1, 3), 1, 6,
                                         exact = FALSE)),
                 suppressWarnings(qdelap(c(NA, 4), c(1, 3), 1, 6)))
expect_warning(qdelap(c(5, NaN), c(1, 3), 1, 6, exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(c(5, NaN), c(1, 3), 1, 6,
                                         exact = FALSE)),
                 suppressWarnings(qdelap(c(5, NaN), c(1, 3), 1, 6)))

# Singleton Inf
expect_true(is.infinite(qdelap(1, 3, 1, 2)))
expect_true(is.infinite(qdelap(1, 3, 1, 2, exact = FALSE)))
expect_true(is.infinite(qdelap(5, 1, 2, 3)))
expect_true(is.infinite(qdelap(5, 1, 2, 3, exact = FALSE)))
expect_identical(is.infinite(qdelap(c(1, 3), 3, 1, 2)), rep(TRUE, 2))
expect_identical(is.infinite(qdelap(c(1, 3), 3, 1, 2, exact = FALSE)),
                 rep(TRUE, 2))
# Vector Inf
expect_identical(is.infinite(qdelap(1:2, 3, c(1, 1), 2)), rep(TRUE, 2))
expect_identical(is.infinite(qdelap(1:3, c(2, 1, 2), c(1, 6, 2), c(1, 2, 0.4))),
                 rep(TRUE, 3))
expect_warning(qdelap(1:2, 3, c(1, 1), 2, exact = FALSE), inpWarn)
expect_identical(suppressWarnings(qdelap(1:2, 3, c(1, 1), 2, exact = FALSE)),
                 qdelap(1:2, 3, c(1, 1), 2))
expect_warning(qdelap(1:3, c(2, 1, 2), c(1, 6, 2), c(1, 2, 0.4), exact = FALSE),
               inpWarn)
expect_identical(suppressWarnings(qdelap(1:3, c(2, 1, 2), c(1, 6, 2),
                                         c(1, 2, 0.4), exact = FALSE)),
                 qdelap(1:3, c(2, 1, 2), c(1, 6, 2), c(1, 2, 0.4)))

# Approximate throws error when nonpositive is passed
expect_warning(qdelap(0.1, 0, 2, 3, exact = FALSE), nanWarn)
expect_warning(qdelap(0.1, 1, 0, 3, exact = FALSE), nanWarn)
expect_warning(qdelap(0.1, 1, 2, -3, exact = FALSE), nanWarn)

# Approximate throws warning when parameter vectors are passed and is equal to
# exact.
expect_warning(qdelap(c(0.4, 0.07), c(1, 2), c(4, 1), c(2, 5), exact = FALSE),
               "Quantile approximation relies on pooling")
expect_identical(suppressWarnings(qdelap(c(0.4, 0.07), c(1, 2), c(4, 1),
                                         c(2, 5), exact = FALSE)),
                 qdelap(c(0.4, 0.07), c(1, 2), c(4, 1), c(2, 5)))

# Non-double parameters converted
expect_equal(qdelap(0.25, 1L, 2L, 3L), qdelap(0.25, 1, 2, 3), tolerance = tol)

# Zero-length inputs return numeric(0) in BOTH exact and approximate branches.
expect_identical(qdelap(0.5, numeric(0), 2, 3), numeric(0))
expect_identical(qdelap(numeric(0), 1, 2, 3), numeric(0))
expect_identical(qdelap(0.5, numeric(0), 2, 3, exact = FALSE), numeric(0))
expect_identical(qdelap(numeric(0), 1, 2, 3, exact = FALSE), numeric(0))

# Approximate branch preserves input order for mixed edge/valid probabilities
# (Issue #3: previously bucketed as c(qNeg, q0, qValid, qInf), scrambling the
# result whenever categories were interleaved).
set.seed(9174L)
mixedApprox <- qdelap(c(0.9, 0, 0.1), 4, 6, 10, exact = FALSE)
mixedExact <- qdelap(c(0.9, 0, 0.1), 4, 6, 10, exact = TRUE)
expect_identical(mixedApprox[2L], 0)
expect_true(mixedApprox[1L] > mixedApprox[3L])
expect_equal(mixedApprox, mixedExact, tolerance = 0.02)

# Edge values land in the right positions even with p < 0 and p >= 1 present.
set.seed(9174L)
edgeApprox <- suppressWarnings(
  qdelap(c(1, 0.5, -0.2, 0), 4, 6, 10, exact = FALSE)
)
expect_identical(edgeApprox[1L], Inf)
expect_identical(edgeApprox[3L], NaN)
expect_identical(edgeApprox[4L], 0)
expect_true(is.finite(edgeApprox[2L]))
expect_true(edgeApprox[2L] > 0)
expect_warning(qdelap(c(1, 0.5, -0.2, 0), 4, 6, 10, exact = FALSE), nanWarn)

# All-edge-case input skips simulation but still returns correct positions.
expect_identical(suppressWarnings(qdelap(c(0, 1, -1), 4, 6, 10, exact = FALSE)),
                 c(0, Inf, NaN))

# qdelap builds its CDF lookup with the same recurrence and accumulation
# order as pdelap, so the round trip must hold exactly, including in the
# underflow-scaled regime and for deep percentiles.
# (Round trips are only possible while the CDF is strictly below 1 in double
# precision and still strictly increasing; past saturation every k maps to
# the same stored CDF, in any implementation.)
testV <- c(0, 1, 13, 14, 25, 60, 90)
expect_equal(qdelap(pdelap(testV, 5, 2, 3), 5, 2, 3), testV, tolerance = tol)
testV <- c(700, 800, 900)
expect_equal(qdelap(pdelap(testV, 5, 3, 800), 5, 3, 800), testV,
             tolerance = tol)

# A quantile at survival 1e-12 with K ~ 11700 support points is
# ill-conditioned for any linear-space CDF accumulation (total rounding
# ~ K * eps ~ 2.6e-12 exceeds the target), so cross-path equality cannot be
# asserted there. Assert accuracy instead, using the directly-summed upper
# tail as the oracle: the returned quantile must land where the true
# survival crosses 1e-12.
qDeep <- qdelap(1 - 1e-12, 50, 100, 100)
expect_true(pdelap(qDeep, 50, 100, 100, lower.tail = FALSE) <= 2e-12)
expect_true(pdelap(qDeep - 100, 50, 100, 100, lower.tail = FALSE) > 1e-12)

# The geometrically-doubled lookup table must return the identical quantiles
# as the elemental search path (forced by a vector-valued parameter).
set.seed(495L)
ptst <- c(runif(200), 1e-14, 1 - 1e-13)
expect_identical(qdelap(ptst, 4, 6, 10), qdelap(ptst, c(4, 4), 6, 10))

# NaN singleton parameters previously slipped past the <= 0 screens and hung
# the build loop while eating memory; they must now return NaN immediately.
expect_warning(vapply(list(c(NaN, 2, 3), c(2, NaN, 3), c(2, 3, NaN)),
                      function(g) qdelap(0.5, g[1], g[2], g[3]), double(1)),
               nanWarn)
expect_identical(suppressWarnings(qdelap(c(0.1, 0.5), NaN, 2, 3)),
                 rep(NaN, 2))

# Parameters beyond TBLMAXCOEF route to the legacy incremental build and must
# agree with the elemental search path.
expect_equal(qdelap(c(0.05, 0.5, 0.95), 1e-28, 1e31, 2),
             qdelap(c(0.05, 0.5, 0.95), c(1e-28, 1e-28), 1e31, 2),
             tolerance = tol)

# Specialty test to bring coverage to 100%
# expect_identical(qdelap(0.1, 1e10, 1e10, 1e10), Inf)
# qdelap(0.01, 1e-20, 1e6, 1e24) HANGS

# Restore original thread count
setDelapThreads(oldThreads)
