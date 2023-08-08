# Copyright (c) 2023, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

expect_silent(setDelapThreads(4L))
expect_identical(getDelapThreads(), 4L)
expect_error(setDelapThreads(0L), "Number of threads must be >= 0.")
expect_error(setDelapThreads(-2L), "Number of threads must be >= 0.")
expect_silent(setDelapThreads(3.9))
expect_identical(getDelapThreads(), 3L)
capMsg <- paste0("Capping at system maximum of ", getOption("DLPCPU"), ".")
expect_message(setDelapThreads(1024^2), capMsg)
expect_identical(getDelapThreads(), getOption("DLPCPU"))
