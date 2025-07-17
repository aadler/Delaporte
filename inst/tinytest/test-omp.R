# Copyright (c) 2023, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

expect_silent(setDelapThreads(1L))
expect_identical(getDelapThreads(), 1L)
expect_error(setDelapThreads(0L), "Number of threads must be > 0.")
expect_error(setDelapThreads(-2L), "Number of threads must be > 0.")
expect_silent(setDelapThreads(1.9))
expect_identical(getDelapThreads(), 1L)
capMsg <- paste0("Capping at system maximum of ",
                 get("DLPCPU", envir = DelaporteEnv), ".")
expect_message(setDelapThreads(1024^2), capMsg)
