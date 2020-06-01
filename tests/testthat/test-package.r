context("Package Maintenance")
test_that("Version", {
  expect_match(toBibtex(citation('Delaporte')),
               as.character(packageVersion('Delaporte')),
               fixed = TRUE, all = FALSE)
  expect_match(scan("../../README.md", what = 'character'),
               as.character(packageVersion('Delaporte')),
               fixed = TRUE, all = FALSE)
})