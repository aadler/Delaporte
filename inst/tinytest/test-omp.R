# Copyright (c) 2023, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

expect_silent(setDelapThreads(1L))
expect_identical(getDelapThreads(), 1L)
expect_error(setDelapThreads(0L), "Number of threads must be > 0.")
expect_error(setDelapThreads(-2L), "Number of threads must be > 0.")
expect_error(setDelapThreads(NA_integer_), "Number of threads must be > 0.")
expect_silent(setDelapThreads(1.9))
expect_identical(getDelapThreads(), 1L)

# State lives in the package namespace (per WRE r88406), never in the
# user's workspace and never in the global options registry.
expect_false(exists("DelaporteEnv", envir = globalenv(), inherits = FALSE))
expect_true(is.environment(Delaporte:::DelaporteEnv))
expect_false(any(startsWith(names(options()), "Delap")))

ncpus <- get("DLPCPU", envir = Delaporte:::DelaporteEnv)
if (!is.na(ncpus)) {
  capMsg <- paste0("Capping at system maximum of ", ncpus, ".")
  expect_message(setDelapThreads(1024 ^ 2), capMsg)
  expect_identical(getDelapThreads(), ncpus)
}

# Thread-count independence (added: this invariant was previously untested).
# ddelap_f/pdelap_f/qdelap_f all carry an OMP parallel-do region, but only in
# their *vector-parameter* branch (na > 1 .or. nb > 1 .or. nl > 1); the
# singleton fast-table path is inherently sequential and never touches OMP.
# rdelap_f has no OMP region of its own but delegates to qdelap_f, so it
# inherits the same exposure through vector parameters. Every other test file
# forces setDelapThreads(2L) for its whole run, which exercises the parallel
# branch, but none of them ever compared output at different thread counts
# for equality -- so a race or non-deterministic reduction order in the OMP
# region could have silently passed the whole suite.
oldThreads <- getDelapThreads()
xInd <- 0:200
aInd <- c(1.5, 2.75, 0.3)
bInd <- c(0.4, 1.1, 3.2)
lInd <- c(2, 0.05, 10)
pInd <- seq(0.01, 0.99, by = 0.01)

setDelapThreads(1L)
d1 <- ddelap(xInd, aInd, bInd, lInd)
p1 <- pdelap(xInd, aInd, bInd, lInd)
q1 <- qdelap(pInd, aInd, bInd, lInd)
set.seed(42)
r1 <- rdelap(50, aInd, bInd, lInd)

# 2 threads is enough to exercise the OMP-parallel branch vs. the threads = 1L
# serial baseline above, and stays within CRAN Repository Policy's cap of 2
# simultaneous threads/cores for check-farm jobs -- no need to request more
# and no need to special-case CRAN.
suppressMessages(setDelapThreads(2L)) # capping (if any) already tested above
d2 <- ddelap(xInd, aInd, bInd, lInd)
p2 <- pdelap(xInd, aInd, bInd, lInd)
q2 <- qdelap(pInd, aInd, bInd, lInd)
set.seed(42)
r2 <- rdelap(50, aInd, bInd, lInd)

expect_identical(d1, d2)
expect_identical(p1, p2)
expect_identical(q1, q2)
expect_identical(r1, r2)

setDelapThreads(oldThreads)
