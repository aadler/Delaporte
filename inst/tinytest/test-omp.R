# Copyright (c) 2023, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

expect_silent(setThreads(4L))
expect_identical(getThreads(), 4L)
expect_error(setThreads(0L), "Number of threads must be >= 0.")
expect_error(setThreads(-2L), "Number of threads must be >= 0.")
expect_silent(setThreads(3.9))
expect_identical(getThreads(), 3L)
capMsg <- paste0("Capping at system maximum of ", getOption("DLPCPU"), ".")
expect_message(setThreads(1024^2), capMsg)
expect_identical(getThreads(), getOption("DLPCPU"))
