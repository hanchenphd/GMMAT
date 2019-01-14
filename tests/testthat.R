library(testthat)
library(GMMAT)
Sys.setenv(MKL_NUM_THREADS = 1)

test_check("GMMAT")
