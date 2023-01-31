# Test Citation
expect_true(any(grepl(packageVersion("Delaporte"),
                      toBibtex(citation("Delaporte")), fixed = TRUE)))
