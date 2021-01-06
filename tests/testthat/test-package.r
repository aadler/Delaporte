test_that("Citation", {
  expect_match(toBibtex(citation('Delaporte')),
               as.character(packageVersion('Delaporte')),
               fixed = TRUE, all = FALSE)
})

skip("Only run manually; different file structure.")
test_that("README", {
  expect_match(scan("../../README.md", what = 'character', skip = 23L,
                    nlines = 1L, quiet = TRUE),
               as.character(packageVersion('Delaporte')),
               fixed = TRUE, all = FALSE)
  expect_match(scan("../../README.md", what = 'character', skip = 33L,
                    nlines = 1L, quiet = TRUE),
               as.character(packageVersion('Delaporte')),
               fixed = TRUE, all = FALSE)
})
