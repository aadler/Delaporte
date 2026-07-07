# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# For CRAN
oldThreads <- getDelapThreads()
setDelapThreads(2L)

tol <- sqrt(.Machine$double.eps)
nonIntErr <- "Non-integers passed to ddelap. These will have 0 probability."
nanWarn <- "NaNs produced"
VAL <- data.frame(read.csv(file = file.path(".", "RawTest.csv"), header = TRUE))

# Singleton function accuracy
expect_equal(ddelap(0:36, 1, 4, 2), VAL$DDELAP_1, tolerance = tol)

# alpha < 0.8
expect_equal(ddelap(4L, 0.5, 4, 0.2), 0.0547024400602606, tolerance = tol)

# Singleton log
expect_equal(ddelap(0:36, 5, 3, 8, log = TRUE), log(ddelap(0:36, 5, 3, 8)),
             tolerance = tol)

# Singleton NA
expect_warning(ddelap(1L, NA, 2, 3), nanWarn)
expect_identical(suppressWarnings(ddelap(1:3, 4, NA, 3)), rep(NaN, 3))

# Singleton NaN
expect_warning(ddelap(1L, 0, 1, 2), nanWarn)
expect_warning(ddelap(1:10, 0, 1, 2), nanWarn)
expect_warning(ddelap(1L, -2, 1, 2), nanWarn)
expect_warning(ddelap(0L, 1, 0, 2), nanWarn)
expect_warning(ddelap(0L, 1, -4, 2), nanWarn)
expect_warning(ddelap(0L, 1, 4, 0), nanWarn)
expect_warning(ddelap(0L, 1, 4, -3), nanWarn)
expect_warning(ddelap(0L, 1, 4, -3), nanWarn)
expect_warning(ddelap(NaN, 1, 4, 6), nanWarn)
expect_warning(ddelap(NA, 1, 4, 12), nanWarn)
tst <- suppressWarnings(ddelap(c(NA, 4, NaN), 0.5, 4, 0.2))
expect_equal(tst, c(NaN, ddelap(4L, 0.5, 4, 0.2), NaN), tolerance = tol)

# Vector function accuracy
expect_equal(ddelap(0:36, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7)),
             VAL$DDELAP_Triple, tolerance = tol)
# Vector log
expect_equal(ddelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9), log = TRUE),
             log(ddelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9))),
             tolerance = tol)

# Vector NA
expect_identical(suppressWarnings(ddelap(1:3, c(4, 1, 2), c(1, 5, 3), NA)),
                 rep(NaN, 3))
tst <- suppressWarnings(ddelap(c(4, 4, 4), c(1, 0.5, NA), 4, c(NaN, 0.2, 4)))
expect_equal(tst, c(NaN, ddelap(4, 0.5, 4, 0.2), NaN), tolerance = tol)

# Vector NaN
expect_warning(ddelap(1:3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0)), nanWarn)
expect_warning(ddelap(1:3, c(-3, 1, 2), c(1, -5, 2), c(1, 2, -0.1)), nanWarn)
expect_warning(ddelap(c(NA, 2), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
expect_warning(ddelap(c(3, NaN), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
tst <- suppressWarnings(ddelap(c(NA, 0, NaN), c(1, 2, 3), c(4, 1, 2),
                               c(2, 5, 7)))
expect_equal(tst, c(NaN, ddelap(0, 2, 1, 5), NaN), tolerance = tol)

# Non-integer warning
expect_warning(ddelap(1.1, 1, 2, 3), nonIntErr)
expect_warning(ddelap(c(1, 1.1, 1.2, 3), c(1, 1), 2, 3), nonIntErr)
expect_warning(ddelap(seq(2, 3, 0.1), c(1, 1), 2, 3), nonIntErr)

# Non-double parameters converted
expect_equal(ddelap(2L, 1L, 2L, 3L), ddelap(2L, 1, 2, 3), tolerance = tol)

# Infinite values
expect_identical(ddelap(Inf, 1L, 2L, 3L), 0)
expect_identical(ddelap(c(Inf, Inf), c(1L, 2L), 2L, 3L), c(0, 0))
expect_warning(ddelap(-Inf, 1L, 2L, 3L), nanWarn)

# Test log1p using Taylor branch; only used by beta parameter.
expect_equal(ddelap(1, 1, 1e-10, 2), 0.270670566459692, tolerance = tol)

# Zero-length parameters must return numeric(0), not crash R (Issue: SIGFPE
# via integer division by zero in imk when a parameter has length 0).
expect_identical(ddelap(0:3, numeric(0), 1, 2), numeric(0))
expect_identical(ddelap(0:3, 1, numeric(0), 2), numeric(0))
expect_identical(ddelap(0:3, 1, 2, numeric(0)), numeric(0))
expect_identical(ddelap(numeric(0), 1, 2, 3), numeric(0))
expect_identical(ddelap(integer(0), 1, 2, 3, log = TRUE), numeric(0))

# log = TRUE is computed in log space via log-sum-exp, so deep-tail
# log-probabilities remain finite even where the linear-space PMF underflows
# below the smallest representable double (~1e-308). The oracle assembles the
# log-PMF from the definition of the Delaporte as a negative binomial
# convolved with a Poisson, keeping every term in log space.
lddelapOracle <- function(k, a, b, l) {
  j <- seq.int(0L, k)
  lt <- dnbinom(j, size = a, prob = 1 / (1 + b), log = TRUE) +
    dpois(k - j, l, log = TRUE)
  m <- max(lt)
  m + log(sum(exp(lt - m)))
}

# Singleton deep tail: linear space underflows here, log space must not
expect_equal(ddelap(2000, 1, 1, 1, log = TRUE), lddelapOracle(2000, 1, 1, 1),
             tolerance = tol)
expect_equal(ddelap(500, 0.5, 4, 0.2, log = TRUE),
             lddelapOracle(500, 0.5, 4, 0.2), tolerance = tol)

# Vector-parameter deep tail
expect_equal(ddelap(c(2000, 1500), c(1, 2), c(1, 3), c(1, 2), log = TRUE),
             c(lddelapOracle(2000, 1, 1, 1), lddelapOracle(1500, 2, 3, 2)),
             tolerance = tol)

# Closed form at x = 0: log P(0) = -lambda - alpha * log1p(beta)
expect_equal(ddelap(0, 4, 5, 6, log = TRUE), -6 - 4 * log1p(5),
             tolerance = tol)

# Structure zero cases keep their log-space image of -Inf
expect_warning(ddelap(1.5, 1, 1, 1, log = TRUE), nonIntErr)
expect_identical(suppressWarnings(ddelap(1.5, 1, 1, 1, log = TRUE)), -Inf)
expect_identical(ddelap(Inf, 1, 2, 3, log = TRUE), -Inf)

# Test NaN Return
expect_identical(suppressWarnings(ddelap(3, -1, 2, 3, log = TRUE)), NaN)
expect_warning(ddelap(3, -1, 2, 3, log = TRUE), nanWarn)

# Scalar-parameter fast path (O(K) recurrence table + lookup) must agree with
# the per-element summation path (forced by a vector-valued parameter) in
# both linear and log space, including non-integer observations and the
# deep-tail log values whose linear-space masses underflow to 0.
xtst <- c(0:400, 1000, 2000, 3.5, 7.2)
expect_equal(suppressWarnings(ddelap(xtst, 4, 6, 10)),
             suppressWarnings(ddelap(xtst, c(4, 4), 6, 10)), tolerance = 1e-12)
lFast <- suppressWarnings(ddelap(xtst, 4, 6, 10, log = TRUE))
lElem <- suppressWarnings(ddelap(xtst, c(4, 4), 6, 10, log = TRUE))
expect_equal(lFast[is.finite(lElem)], lElem[is.finite(lElem)],
             tolerance = 1e-10)
expect_identical(is.finite(lFast), is.finite(lElem))

# Log fast path must preserve the log-space guarantee: finite log-PMF where
# the linear PMF underflows (regression for the downward-rescaling branch of
# ddelap_table; without it these return -Inf).
expect_equal(ddelap(2000, 1, 1, 1, log = TRUE),
             ddelap(2000, c(1, 1), 1, 1, log = TRUE), tolerance = 1e-10)
expect_true(is.finite(ddelap(5000, 1, 1, 1, log = TRUE)))

# Underflow-scaled regime through the fast path: normalization and one point
# against the defining NB (*) Poisson convolution.
expect_equal(sum(ddelap(0:4000, 5, 3, 800)), 1, tolerance = 1e-9)
expect_equal(ddelap(850, 5, 3, 800),
             sum(dnbinom(0:850, size = 5, prob = 0.25) *
                   dpois(850 - (0:850), 800)), tolerance = 1e-9)

# Parameters beyond TBLMAXCOEF route around the table to the per-element
# summation and must agree with the vector-parameter path.
expect_equal(suppressWarnings(ddelap(0:5, 1e-28, 1e31, 2)),
             suppressWarnings(ddelap(0:5, c(1e-28, 1e-28), 1e31, 2)),
             tolerance = tol)

# Degenerate parameters whose PMF can fall by more than the table's rescue
# band in one recurrence step (here the guaranteed step ratio is
# lambda / (k + 1) ~ 1e-54) must route the log-space request to the elemental
# path. The log-PMF is finite - alpha ~ 0 collapses the NB factor to a point
# mass at zero, leaving essentially a Poisson(lambda) - not -Inf.
expect_equal(ddelap(1e4, 1e-50, 1e-100, 1e-50, log = TRUE),
             ddelap(1e4, c(1e-50, 1e-50), 1e-100, 1e-50, log = TRUE),
             tolerance = tol)
expect_equal(ddelap(1e4, 1e-50, 1e-100, 1e-50, log = TRUE),
             dpois(1e4, 1e-50, log = TRUE), tolerance = 1e-9)

# Parameters beyond TBLMAXCOEF route the log-space request around the table
# as well and must agree with the per-element path.
expect_equal(suppressWarnings(ddelap(0:5, 1e-28, 1e31, 2, log = TRUE)),
             suppressWarnings(ddelap(0:5, c(1e-28, 1e-28), 1e31, 2,
                                     log = TRUE)), tolerance = tol)

# Invalid parameters through the per-element log path return NaN.
expect_true(all(is.nan(suppressWarnings(
  ddelap(0:2, c(-1, -1), 2, 3, log = TRUE)))))

# Specialty tests to bring coverage to 100%
# Trigger lpmf = ieee_value(x, ieee_quiet_nan) in ddelap_f_s_log
expect_warning(ddelap(2e30, 1, 1, NaN, log = TRUE), nanWarn)

# Non-finite parameters are invalid: any infinite parameter implies an
# infinite-mean distribution with no mass at any finite point, so the result
# is NaN with a warning - previously Inf-driven NaNs were laundered into
# exactly 0 and 1 by the floating-point clamps. Note that 1e3000 overflows
# double precision and IS Inf.
expect_warning(ddelap(1, Inf, Inf, Inf, log = TRUE), nanWarn)
expect_true(is.nan(suppressWarnings(ddelap(1, Inf, Inf, Inf, log = TRUE))))
expect_true(is.nan(suppressWarnings(ddelap(1, Inf, Inf, Inf))))
expect_true(is.nan(suppressWarnings(ddelap(1, 4, Inf, 0.1))))
expect_true(is.nan(suppressWarnings(ddelap(1, c(Inf, Inf), 2, 3))))

# x at or beyond MAXD (huge 64-bit integer), including x = +Inf, carries no
# reportable mass by the package's summation cap: 0 in linear space, -Inf in
# log space, computed without converting x to an integer kind it overflows.
expect_identical(ddelap(Inf, 1, 2, 3), 0)
expect_identical(ddelap(2e30, 1, 1, 1, log = TRUE), -Inf)

# Restore original thread count
setDelapThreads(oldThreads)
