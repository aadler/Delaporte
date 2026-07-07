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
