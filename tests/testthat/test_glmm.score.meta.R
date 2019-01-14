context("single-variant score test meta-analysis")

test_that("single-variant meta-analysis", {
	infile1 <- system.file("extdata", "meta1.txt", package = "GMMAT")
	infile2 <- system.file("extdata", "meta2.txt", package = "GMMAT")
	infile3 <- system.file("extdata", "meta3.txt", package = "GMMAT")
	outfile <- tempfile()
	glmm.score.meta(files = c(infile1, infile2, infile3), outfile = outfile, SNP = rep("SNP", 3), A1 = rep("A1", 3), A2 = rep("A2", 3))
	out <- read.table(outfile, header = TRUE, as.is = TRUE)
	expect_equal(signif(out$PVAL, digits = 5), signif(c(3.076661e-01, 4.661635e-01, 2.996084e-01, 2.716879e-01, 4.661635e-01, 9.775175e-01, 1, 3.076661e-01, 3.076661e-01, 5.789150e-06), digits = 5))
	unlink(outfile)
})
