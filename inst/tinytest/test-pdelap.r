# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# For CRAN
oldThreads <- getDelapThreads()
setDelapThreads(2L)

tol <- sqrt(.Machine$double.eps)
VAL <- data.frame(read.csv(file = file.path(".", "RawTest.csv"), header = TRUE))
nanWarn <- "NaNs produced"

# Singleton function accuracy
expect_equal(pdelap(0:36, 2, 1, 5), VAL$PDELAP_2, tolerance = tol)

# Singleton log.p
expect_equal(pdelap(0:36, 4, 5, 1, log.p = TRUE), log(pdelap(0:36, 4, 5, 1)),
             tolerance = tol)

# Singleton lower.tail
expect_equal(pdelap(0:100, 8, 10, 6, lower.tail = FALSE),
             1 - pdelap(0:100, 8, 10, 6), tolerance = tol)
expect_equal(pdelap(6, 2.9647, 0.005 / 2.9647, 0.0057, lower.tail = FALSE), 0,
             tolerance = tol)

# Singleton lower.tail & log.p
expect_equal(pdelap(0:100, 8, 10, 6, lower.tail = FALSE, log.p = TRUE),
             log(1 - pdelap(0:100, 8, 10, 6)), tolerance = tol)

# Singleton NaN
expect_warning(pdelap(1, 0, 1, 2), nanWarn)
expect_warning(pdelap(1, -85, 1, 2), nanWarn)
expect_warning(pdelap(0, 1, 0, 2), nanWarn)
expect_warning(pdelap(0, 1, -1e4, 2), nanWarn)
expect_warning(pdelap(0, 1, 4, 0), nanWarn)
expect_warning(pdelap(0, 1, 4, -1e-3), nanWarn)
expect_warning(ddelap(NaN, 1, 4, 6), nanWarn)
expect_warning(ddelap(NA, 1, 4, 12), nanWarn)
tst <- suppressWarnings(pdelap(c(NA, 4, NaN), 0.5, 4, 0.2))
expect_equal(tst, c(NaN, pdelap(4, 0.5, 4, 0.2), NaN), tolerance = tol)
expect_identical(suppressWarnings(pdelap(c(NA, NaN), 0.5, 4, 0.2)), rep(NaN, 2))

# Vector function accuracy
expect_equal(pdelap(0:36, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7)),
             VAL$PDELAP_Triple, tolerance = tol)

# Vector log.p
expect_equal(pdelap(0:36, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7), log.p = TRUE),
             log(VAL$PDELAP_Triple), tolerance = tol)

# Vector lower.tail
expect_equal(pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9),
                    lower.tail = FALSE),
             1 - pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9)),
             tolerance = tol)

# Vector lower.tail & log.p
expect_equal(pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9),
                    lower.tail = FALSE, log.p = TRUE),
             log(1 - pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9))),
             tolerance = tol)

# Vector NaN
expect_warning(pdelap(1:10, 0, 1, 2), nanWarn)
expect_warning(pdelap(1:3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0)), nanWarn)
expect_warning(pdelap(1:3, c(-5e5, 1, 2), c(1, -2, 2), c(1, 2, -8e-4)), nanWarn)
expect_warning(pdelap(c(NA, 2), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
expect_warning(pdelap(c(3, NaN), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
tst <- suppressWarnings(pdelap(c(NA, 0, NaN), c(1, 2, 3), c(4, 1, 2),
                               c(2, 5, 7)))
expect_equal(tst, c(NaN, pdelap(0, 2, 1, 5), NaN), tolerance = tol)
tst <- suppressWarnings(pdelap(c(0, 0, 0), c(NA, 2, 3), c(4, 1, 2),
                               c(2, 5, NaN)))
expect_equal(tst, c(NaN, pdelap(0, 2, 1, 5), NaN), tolerance = tol)

# Negative values due to floating point issues are 0
expect_equal(pdelap(500, 13.08251, 0.02414521, 0.04421658, FALSE, FALSE), 0,
             tolerance = tol)

# Non-double parameters converted
expect_equal(pdelap(2L, 1L, 2L, 3L), pdelap(2L, 1, 2, 3), tolerance = tol)

# Floating point issues do not lead to CDF > 1
# print(pdelap(1000, 8, 15, 100), digits = 17) used to be 1.0000000000001035
expect_true(pdelap(1000, 8, 15, 100) <= 1)

# Positive infinite arguments
expect_identical(pdelap(Inf, 1L, 2L, 3L), 1)
expect_identical(pdelap(c(Inf, Inf), c(1L, 2L), 2L, 3L), c(1, 1))

# Negative and -Inf arguments follow base R: below the support the CDF is 0
# (-Inf on the log scale) and the survival is 1 (0 on the log scale), all
# returned silently. Previously NaN with a warning.
expect_identical(pdelap(-Inf, 1L, 2L, 3L), 0)
expect_identical(pdelap(-1, 1, 2, 3), 0)
expect_identical(pdelap(-Inf, 1L, 2L, 3L, log.p = TRUE), -Inf)
expect_identical(pdelap(-1, 1, 2, 3, log.p = TRUE), -Inf)
expect_identical(pdelap(-1, 1, 2, 3, lower.tail = FALSE), 1)
expect_identical(pdelap(-Inf, 1, 2, 3, lower.tail = FALSE), 1)
expect_identical(pdelap(-1, 1, 2, 3, lower.tail = FALSE, log.p = TRUE), 0)
expect_identical(pdelap(-Inf, 1, 2, 3, lower.tail = FALSE, log.p = TRUE), 0)

# Mixed vector keeps input order and stays silent. A negative entry forces the
# per-element path, whose values differ from the singleton fast-path table by at
# most an ulp, so compare with tolerance. Below-support entries are exactly 0.
expect_equal(pdelap(c(-1, 5, -Inf, 10), 1, 2, 3),
             c(0, pdelap(5, 1, 2, 3), 0, pdelap(10, 1, 2, 3)), tolerance = tol)
expect_identical(pdelap(c(-1, 5, -Inf, 10), 1, 2, 3)[c(1, 3)], c(0, 0))
expect_silent(pdelap(c(-1, -Inf, 0, 5), 1, 2, 3))

# Zero-length inputs return numeric(0) (SIGFPE regression guard).
expect_identical(pdelap(0:3, numeric(0), 1, 2), numeric(0))
expect_identical(pdelap(numeric(0), 1, 2, 3), numeric(0))
expect_identical(pdelap(0:3, 1, 2, numeric(0), lower.tail = FALSE), numeric(0))

# Deep upper tail is computed by direct summation, not 1 - CDF, so survival
# probabilities below machine epsilon retain full relative accuracy. The oracle
# uses the definition of the Delaporte as the convolution of a negative binomial
# with a Poisson: P(X > k) = sum_j dnbinom(j) * P(Pois > k - j), where ppois
# computes its own upper tail accurately.
sdelapOracle <- function(k, a, b, l, J = 5000L) {
  j <- seq.int(0L, J)
  sum(dnbinom(j, size = a, prob = 1 / (1 + b)) *
        ppois(k - j, l, lower.tail = FALSE))
}

# Singleton (fast path) accuracy across shallow, deep, and very deep tails.
# Relative tolerance is meaningful here since targets are nonzero.
qDeep <- c(0, 5, 30, 60, 120)
expect_equal(pdelap(qDeep, 1, 1, 1, lower.tail = FALSE),
             vapply(qDeep, sdelapOracle, double(1), a = 1, b = 1, l = 1),
             tolerance = tol)

# Vector-parameter (slow path) accuracy in the deep tail
expect_equal(pdelap(c(60, 80), c(1, 2), c(1, 3), c(1, 2), lower.tail = FALSE),
             c(sdelapOracle(60, 1, 1, 1), sdelapOracle(80, 2, 3, 2)),
             tolerance = tol)

# log.p on the deep upper tail no longer collapses to log of rounding noise
expect_equal(pdelap(60, 1, 1, 1, lower.tail = FALSE, log.p = TRUE),
             log(sdelapOracle(60, 1, 1, 1)), tolerance = tol)

# log.p in the astronomically deep tail, where the linear-space value
# underflows to exactly 0 (survival ~ exp(-1386), far below the smallest
# representable double): log() of that would be -Inf, but the true log-CDF/
# log-survival is finite. Oracle built independently via log-sum-exp over
# already-validated ddelap(..., log = TRUE) points, not via pdelap itself.
logSumExpOracle <- function(lx) {
  m <- max(lx)
  m + log(sum(exp(lx - m)))
}
deepUpperOracle <- logSumExpOracle(ddelap(2001:20000, 1, 1, 1, log = TRUE))
expect_equal(pdelap(2000, 1, 1, 1, lower.tail = FALSE, log.p = TRUE),
             deepUpperOracle, tolerance = tol)
expect_false(is.infinite(pdelap(0, 1, 1, 1000, log.p = TRUE)))
expect_equal(pdelap(0, 1, 1, 1000, log.p = TRUE), log(0.5) - 1000,
             tolerance = tol)

# Same two deep-tail cases via the vector-recycling fallback path (na > 1
# forces this branch instead of the singleton table build)
expect_equal(pdelap(2000, c(1, 1), 1, 1, lower.tail = FALSE, log.p = TRUE),
             deepUpperOracle, tolerance = tol)
expect_equal(pdelap(0, c(1, 1), 1, 1000, log.p = TRUE), log(0.5) - 1000,
             tolerance = tol)

# Moderate tail (survival comfortably above sqrt(EPS), where the complement
# 1 - CDF is used instead of direct summation after the TAILSWITCH threshold
# change): linear and log.p agree with each other and with the oracle.
qMod <- c(50, 55, 60, 65, 70)
modOracle <- vapply(qMod, sdelapOracle, double(1), a = 4, b = 6, l = 10)
expect_equal(pdelap(qMod, 4, 6, 10, lower.tail = FALSE), modOracle,
             tolerance = tol)
expect_equal(pdelap(qMod, 4, 6, 10, lower.tail = FALSE, log.p = TRUE),
             log(modOracle), tolerance = tol)

# Performance guard: sdelap_f_s's remainder-bound floor forces roughly
# 36 * beta iterations regardless of q once invoked, an O(beta^2) cost per
# call; pdelap(694, 1, 1000, 1, lower.tail = FALSE) took ~27.5s under the old
# CDF > 0.5 invocation threshold (confirmed by execution) because its CDF is
# only ~0.5, nowhere near the tiny-survival regime the direct sum exists for.
# Generous 5s bound avoids flakiness on slow runners while still catching a
# regression back to the old threshold (which would take 20-30s).
perfElapsed <- system.time(
  perfVal <- pdelap(694, 1, 1000, 1, lower.tail = FALSE)
)[["elapsed"]]
expect_true(perfElapsed < 5)
expect_equal(perfVal, sdelapOracle(694, 1, 1000, 1, J = 20000L),
             tolerance = tol)
perfElapsedLog <- system.time(
  perfValLog <- pdelap(694, 1, 1000, 1, lower.tail = FALSE, log.p = TRUE)
)[["elapsed"]]
expect_true(perfElapsedLog < 5)
expect_equal(perfValLog, log(sdelapOracle(694, 1, 1000, 1, J = 20000L)),
             tolerance = tol)

# Survival function is nonincreasing and nonnegative over a long range
survivalCheck <- pdelap(0:400, 1, 1, 1, lower.tail = FALSE)
expect_true(all(diff(survivalCheck) <= 0))
expect_true(all(survivalCheck >= 0))

# Upper-tail edge cases: infinite q survives the direct-summation branch
expect_identical(pdelap(Inf, 1, 1, 1, lower.tail = FALSE), 0)
expect_identical(suppressWarnings(pdelap(NaN, 1, 1, 1, lower.tail = FALSE)),
                 NaN)

# Fast-path upper tail where even the largest q is below the median
# (CDF(max(q)) <= 0.5, here about 0.33): the survival anchor at floor(max(q))
# is computed as the complement 1 - CDF, which is safe from cancellation
# because the result is at least one half, and the backward accumulation
# builds the remaining survival values on top of that anchor.
qLow <- 0:5
expect_equal(pdelap(qLow, 2, 1, 5, lower.tail = FALSE),
             vapply(qLow, sdelapOracle, double(1), a = 2, b = 1, l = 5),
             tolerance = tol)
expect_equal(pdelap(qLow, 2, 1, 5, lower.tail = FALSE),
             1 - pdelap(qLow, 2, 1, 5), tolerance = tol)

# O(K) recurrence table build (ddelap_table): agreement with an independent
# base R oracle built from the defining NB (*) Poisson convolution, across
# ordinary, skewed, and underflow-scaled (lambda + alpha * log1p(beta) > 745)
# parameter regimes. The last two regimes exercise the log-space scaling
# branch, whose bookkeeping error grows like lambda * .Machine$double.eps,
# hence the looser (but still tight) tolerance.
pOracle <- function(q, a, b, l) {
  vapply(q, function(k) {
    i <- 0:k
    sum(dnbinom(i, size = a, prob = 1 / (1 + b)) * ppois(k - i, l))
  }, double(1))
}
qtst <- c(0:20, 50, 100, 400)
for (prm in list(c(4, 6, 10), c(0.001, 1000, 5), c(50, 0.02, 3),
                 c(0.5, 0.5, 0.5))) {
  expect_equal(pdelap(qtst, prm[1], prm[2], prm[3]),
               pOracle(qtst, prm[1], prm[2], prm[3]), tolerance = 1e-14)
}
for (prm in list(c(5, 3, 800), c(2, 7, 1500))) {
  ora <- pOracle(qtst, prm[1], prm[2], prm[3])
  got <- pdelap(qtst, prm[1], prm[2], prm[3])
  keep <- ora > 1e-290                # below this both should agree on ~0
  expect_equal(got[keep], ora[keep], tolerance = 1e-9)
  expect_true(all(got[!keep] < 1e-280))
}

# The recurrence fast path and the per-element summation path (forced by a
# vector-valued parameter) must agree to rounding.
qtst <- 0:400
expect_equal(pdelap(qtst, 4, 6, 10), pdelap(qtst, c(4, 4), 6, 10),
             tolerance = tol)
expect_equal(pdelap(qtst, 4, 6, 10, lower.tail = FALSE),
             pdelap(qtst, c(4, 4), 6, 10, lower.tail = FALSE), tolerance = tol)

# Parameters beyond TBLMAXCOEF route around the recurrence to the legacy
# summation build and must agree with the per-element path.
expect_equal(pdelap(0:6, 1e-28, 1e31, 2), pdelap(0:6, c(1e-28, 1e-28), 1e31, 2),
             tolerance = tol)

# Values far past the old 2^15 gate now compute on the fast path; spot-check
# monotonicity, the CDF limit, and a tail point against the direct survival
# summation.
bigP <- pdelap(c(1e5, 2e5, 5e5), 4, 6, 10)
expect_true(all(diff(bigP) >= 0))
expect_true(bigP[3] == 1)
expect_equal(pdelap(500, 4, 6, 10, lower.tail = FALSE),
             pdelap(500, c(4, 4), 6, 10, lower.tail = FALSE), tolerance = tol)

expect_warning(pdelap(1, Inf, 2, 3), nanWarn)
expect_true(is.nan(suppressWarnings(pdelap(1, Inf, 2, 3))))
expect_true(is.nan(suppressWarnings(pdelap(1, 2, Inf, 3, lower.tail = FALSE))))
expect_true(is.nan(suppressWarnings(pdelap(1, 2, 3, c(Inf, Inf)))))

# Restore original thread count
setDelapThreads(oldThreads)
