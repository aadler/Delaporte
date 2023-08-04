# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

pV <- packageVersion("Delaporte")

# Test CITATION has most recent package version
expect_true(any(grepl(pV, toBibtex(citation("Delaporte")), fixed = TRUE)))

# Test NEWS has most recent package version
expect_true(any(grepl(pV, news(package = "Delaporte"), fixed = TRUE)))

# Test that NEWS has an entry with DESCRIPTION's Date
expect_true(any(grepl(packageDate("Delaporte"), news(package = "Delaporte"),
                      fixed = TRUE)))
