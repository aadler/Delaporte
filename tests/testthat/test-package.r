test_that("Version", {
  expect_match(toBibtex(citation('Delaporte')),
               as.character(packageVersion('Delaporte')),
               fixed = TRUE, all = FALSE)
})