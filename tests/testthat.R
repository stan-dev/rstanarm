library(testthat)
library(rstanarm)
Sys.unsetenv("R_TESTS")
test_check("rstanarm")
