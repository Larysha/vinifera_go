# Run with:  Rscript tests/testthat.R   (from the project root)
# or:        Rscript -e 'testthat::test_dir("tests/testthat")'
library(testthat)
testthat::test_dir(file.path("tests", "testthat"))
