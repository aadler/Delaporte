# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# For CRAN
oldThreads <- getDelapThreads()
setDelapThreads(2L)

# Test the work done to ensure the qdelap(pedelap(x)) == x as much as possible.

tol <- sqrt(.Machine$double.eps)

# ---- F1: native-space searches -------------------------------------------

# (F,F) round trip on the shared table (previously Inf / order-corrupted)
sfd <- pdelap(0:5, 1, 0.1, 0.05, lower.tail = FALSE)
expect_equal(qdelap(sfd, 1, 0.1, 0.05, lower.tail = FALSE), as.numeric(0:5))

# (F,F) deep tail resolves finitely (previously Inf); 17 verified by
# pdelap bracketing: sf(17) <= 1e-18 < sf(16)
expect_identical(qdelap(1e-18, 1, 0.1, 0.05, lower.tail = FALSE), 17)

# (F,T) deep log-survival target stays in log space (previously exp()
# collapsed it); 417 verified by bracketing
expect_identical(qdelap(-1000, 1, 0.1, 0.05, lower.tail = FALSE,
                        log.p = TRUE), 417)

# (T,T) deep lower-tail log target
q <- qdelap(-900, 4, 5, 5, log.p = TRUE)
expect_true(pdelap(q, 4, 5, 5, log.p = TRUE) >= -900 &&
              pdelap(q - 1, 4, 5, 5, log.p = TRUE) < -900)

# (T,T) round trip including the resolvable region
x <- 0:150
lp <- pdelap(x, 2, 0.2, 80, log.p = TRUE)
expect_equal(qdelap(lp, 2, 0.2, 80, log.p = TRUE)[lp < log(1 - tol)],
             as.numeric(x)[lp < log(1 - tol)])

# ---- Near-1 CDF-side routing (mirror of pdelap's TAILSWITCH rule) --------

# Survival targets within sqrt(EPS) of 1 resolve through the forward CDF
# instead of returning 0/garbage; verify sane monotone finite answers
sf <- pdelap(0:60, 20, 40, 0.1, lower.tail = FALSE)
qq <- qdelap(sf, 20, 40, 0.1, lower.tail = FALSE)
expect_true(all(is.finite(qq)) && all(diff(qq) >= 0))
lsf <- pdelap(0:60, 20, 40, 0.1, lower.tail = FALSE, log.p = TRUE)
qq2 <- qdelap(lsf, 20, 40, 0.1, lower.tail = FALSE, log.p = TRUE)
expect_true(all(is.finite(qq2)) && all(diff(qq2) >= 0))

# ---- Boundary matrix (qpois convention, all four modes) ------------------

expect_identical(qdelap(0, 2, 1, 1, lower.tail = FALSE), Inf)
expect_identical(qdelap(1, 2, 1, 1, lower.tail = FALSE), 0)
expect_warning(expect_identical(qdelap(-0.1, 2, 1, 1, lower.tail = FALSE),
                                NaN))
expect_warning(expect_identical(qdelap(1.1, 2, 1, 1, lower.tail = FALSE),
                                NaN))
expect_identical(qdelap(0, 2, 1, 1, log.p = TRUE), Inf)
expect_identical(qdelap(-Inf, 2, 1, 1, log.p = TRUE), 0)
expect_warning(expect_identical(qdelap(0.5, 2, 1, 1, log.p = TRUE), NaN))
expect_identical(qdelap(0, 2, 1, 1, lower.tail = FALSE, log.p = TRUE), 0)
expect_identical(qdelap(-Inf, 2, 1, 1, lower.tail = FALSE, log.p = TRUE),
                 Inf)
expect_warning(expect_identical(
  qdelap(0.5, 2, 1, 1, lower.tail = FALSE, log.p = TRUE), NaN))
expect_identical(qdelap(numeric(0), 2, 1, 1, lower.tail = FALSE),
                 numeric(0))

# ---- Elemental (vector-parameter) paths ----------------------------------

av <- c(1, 2); bv <- c(0.5, 3); lv <- c(1, 4); xv <- c(3, 7)
expect_equal(qdelap(pdelap(xv, av, bv, lv, lower.tail = FALSE),
                    av, bv, lv, lower.tail = FALSE), as.numeric(xv))
expect_equal(qdelap(pdelap(xv, av, bv, lv, log.p = TRUE),
                    av, bv, lv, log.p = TRUE), as.numeric(xv))
expect_equal(qdelap(pdelap(xv, av, bv, lv, lower.tail = FALSE,
                           log.p = TRUE),
                    av, bv, lv, lower.tail = FALSE, log.p = TRUE),
             as.numeric(xv))

# ---- Structural branches -------------------------------------------------

# Heavy tail: S at the moment-based table start exceeds sqrt(EPS), so the
# top seed takes the 1 - CDF flip instead of the direct tail sum
expect_true(is.finite(qdelap(0.3, 0.1, 100, 0.1, lower.tail = FALSE)))
expect_true(is.finite(qdelap(-1.2, 0.1, 100, 0.1, lower.tail = FALSE,
                             log.p = TRUE)))
expect_true(is.finite(qdelap(-1.2, 0.1, 100, 0.1, log.p = TRUE)))

# TBLMINRATIO degenerate route inside the lg == 1 table build
expect_true(is.finite(qdelap(-5, 1, 1e-37, 1e-37, log.p = TRUE)))
expect_true(is.finite(qdelap(-5, 1, 1e-37, 1e-37, lower.tail = FALSE,
                             log.p = TRUE)))

# TBLMAXCOEF corner without the extreme-parameter cost: huge beta with a
# tiny alpha keeps the mean small, so the elemental fallback returns
# instantly; verified identical to the legacy (T,F) path on these params
expect_identical(qdelap(0.5, 1e-30, 1e31, 0.1, lower.tail = FALSE), 0)
expect_identical(qdelap(c(0.9, 0.1), 1e-30, 1e31, 0.1,
                        lower.tail = FALSE), c(0, 0))

# ---- F4: log cumsum caps -------------------------------------------------

# No positive log-probabilities anywhere, either tail, any table length
for (K in c(15, 60, 250)) {
  expect_true(all(pdelap(0:K, 2, 0.2, 80, log.p = TRUE) <= 0))
  expect_true(all(pdelap(0:K, 2, 0.2, 80, log.p = TRUE,
                         lower.tail = FALSE) <= 0))
}

# ---- (T,F) regression sentinel -------------------------------------------

# The default path must remain byte-identical; spot-check against values
# computed from the CDF definition (full-grid identity vs the pre-fix build
# was verified during the review)
p <- pdelap(0:20, 3, 2, 4)
expect_identical(qdelap(p, 3, 2, 4), as.numeric(0:20))

# ---- F5: saturated-band log values from the survival side ----------------

# (T,T) band values are accurate on the log scale (previously accumulated
# noise ~K*EPS: a true ~-4e-24 surfaced as ~-5e-10) and round-trip exactly
x <- 180:230
lp <- pdelap(x, 2, 3, 0.1, log.p = TRUE)
expect_true(all(lp < 0))
expect_equal(qdelap(lp, 2, 3, 0.1, log.p = TRUE), as.numeric(x))

# (F,T) mirror: deep-left log-survival keeps sub-EPS resolution
# (previously the linear flip rounded it to exactly 0)
lsf <- pdelap(20:60, 20, 40, 0.1, lower.tail = FALSE, log.p = TRUE)
expect_true(all(lsf < 0))
expect_equal(qdelap(lsf, 20, 40, 0.1, lower.tail = FALSE, log.p = TRUE),
             as.numeric(20:60))

# All three routes agree in the band to relative ulps (bitwise identity
# between the table and elemental routes has never been a package
# guarantee -- summation orders differ): table build, TBLMINRATIO
# elemental fill, and the recycling elemental path
expect_equal(pdelap(200, 2, 3, 0.1, log.p = TRUE),
             pdelap(200, c(2, 2), 3, 0.1, log.p = TRUE)[1],
             tolerance = 1e-12)
expect_equal(pdelap(5, 1, 1e-37, 1e-37, log.p = TRUE),
             pdelap(5, c(1, 1), 1e-37, 1e-37, log.p = TRUE)[1],
             tolerance = 1e-12)
expect_equal(pdelap(40, 20, 40, 0.1, lower.tail = FALSE, log.p = TRUE),
             pdelap(40, c(20, 20), 40, 0.1, lower.tail = FALSE,
                    log.p = TRUE)[1], tolerance = 1e-12)

# ---- F3: extreme-parameter refusal ----------------------------------------

# Boundaries and resolvable targets still answer through the capped legacy
# route; the tiny-alpha corner (instant, mean ~ 0) is covered above under
# "Structural branches". The refusal itself costs O(cap^2) ~ tens of seconds
# by design, so it runs only locally, not on CRAN.
if (at_home()) {
  expect_warning(
    QF3 <- qdelap(c(0, 0.5, 1), 1, 1e31, 1),
    pattern = "quantile too large to compute exactly"
  )
  expect_identical(QF3, c(0, NaN, Inf))
}

# Restore original thread count
setDelapThreads(oldThreads)
