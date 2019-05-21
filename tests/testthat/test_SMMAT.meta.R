context("variant set mixed model association test (SMMAT)")

test_that("cross-sectional id le 400 binomial", {
	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
	pheno <- rbind(example$pheno, example$pheno[1:100, ])
	pheno$id <- 1:500
	pheno$disease[sample(1:500,20)] <- NA
	pheno$age[sample(1:500,20)] <- NA
	pheno$sex[sample(1:500,20)] <- NA
	pheno <- pheno[sample(1:500,450), ]
	pheno <- pheno[pheno$id <= 400, ]
	kins <- example$GRM
	obj1 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	outfile1 <- tempfile()
	out1 <- SMMAT(obj1, gdsfile, group.file, meta.file.prefix = outfile1, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out1.meta <- SMMAT.meta(outfile1, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out1[, -(1:8)]), signif(out1.meta[, -(1:2)]))
	unlink(c(paste0(outfile1, ".score.*"), paste0(outfile1, ".var.*")))

	skip_on_cran()

	outfile1.tmp <- tempfile()
	out1.tmp <- SMMAT(obj1, gdsfile, group.file, meta.file.prefix = outfile1.tmp, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"), ncores = 2)
	n.files <- 2
	if(Sys.info()["sysname"] == "Windows") n.files <- 1
	out1.meta.tmp <- SMMAT.meta(outfile1.tmp, n.files = n.files, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out1.tmp[, -(1:8)]), signif(out1.meta.tmp[, -(1:2)]))
	expect_equal(out1, out1.tmp)
	expect_equal(out1.meta, out1.meta.tmp)
	unlink(c(paste0(outfile1.tmp, ".score.*"), paste0(outfile1.tmp, ".var.*")))

	obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	outfile2 <- tempfile()
	out2 <- SMMAT(obj2, gdsfile, group.file, meta.file.prefix = outfile2, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out2.meta <- SMMAT.meta(outfile2, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out2[, -(1:8)]), signif(out2.meta[, -(1:2)]))
	unlink(c(paste0(outfile2, ".score.*"), paste0(outfile2, ".var.*")))
})

test_that("cross-sectional id gt 400 binomial", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
	pheno <- rbind(example$pheno, example$pheno[1:100, ])
	pheno$id <- 1:500
	pheno$disease[sample(1:500,20)] <- NA
	pheno$age[sample(1:500,20)] <- NA
	pheno$sex[sample(1:500,20)] <- NA
	pheno <- pheno[sample(1:500,450), ]
	kins <- diag(500)
	kins[1:400, 1:400] <- example$GRM
	rownames(kins) <- colnames(kins) <- 1:500
	obj1 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	outfile1 <- tempfile()
	out1 <- SMMAT(obj1, gdsfile, group.file, meta.file.prefix = outfile1, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out1.meta <- SMMAT.meta(outfile1, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out1[, -(1:8)], digits = 5), signif(out1.meta[, -(1:2)], digits = 5))
	obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	outfile2 <- tempfile()
	out2 <- SMMAT(obj2, gdsfile, group.file, meta.file.prefix = outfile2, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out2.meta <- SMMAT.meta(outfile2, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out2[, -(1:8)], digits = 5), signif(out2.meta[, -(1:2)], digits = 5))
	unlink(c(paste0(outfile1, ".score.*"), paste0(outfile1, ".var.*"), paste0(outfile2, ".score.*"), paste0(outfile2, ".var.*")))
})

test_that("cross-sectional id le 400 gaussian", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
	pheno <- rbind(example$pheno, example$pheno[1:100, ])
	pheno$id <- 1:500
	pheno$disease[sample(1:500,20)] <- NA
	pheno$age[sample(1:500,20)] <- NA
	pheno$sex[sample(1:500,20)] <- NA
	pheno <- pheno[sample(1:500,450), ]
	pheno <- pheno[pheno$id <= 400, ]
	kins <- example$GRM
	obj1 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	outfile1 <- tempfile()
	out1 <- SMMAT(obj1, gdsfile, group.file, meta.file.prefix = outfile1, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out1.meta <- SMMAT.meta(outfile1, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out1[, -(1:8)]), signif(out1.meta[, -(1:2)]))
	obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	outfile2 <- tempfile()
	out2 <- SMMAT(obj2, gdsfile, group.file, meta.file.prefix = outfile2, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out2.meta <- SMMAT.meta(outfile2, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out2[, -(1:8)], digits = 5), signif(out2.meta[, -(1:2)], digits = 5))
	unlink(c(paste0(outfile1, ".score.*"), paste0(outfile1, ".var.*"), paste0(outfile2, ".score.*"), paste0(outfile2, ".var.*")))
})

test_that("cross-sectional id gt 400 gaussian", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
	pheno <- rbind(example$pheno, example$pheno[1:100, ])
	pheno$id <- 1:500
	pheno$disease[sample(1:500,20)] <- NA
	pheno$age[sample(1:500,20)] <- NA
	pheno$sex[sample(1:500,20)] <- NA
	pheno <- pheno[sample(1:500,450), ]
	kins <- diag(500)
	kins[1:400, 1:400] <- example$GRM
	rownames(kins) <- colnames(kins) <- 1:500
	obj1 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	outfile1 <- tempfile()
	out1 <- SMMAT(obj1, gdsfile, group.file, meta.file.prefix = outfile1, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out1.meta <- SMMAT.meta(outfile1, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out1[, -(1:8)]), signif(out1.meta[, -(1:2)]))
	obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	outfile2 <- tempfile()
	out2 <- SMMAT(obj2, gdsfile, group.file, meta.file.prefix = outfile2, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out2.meta <- SMMAT.meta(outfile2, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out2[, -(1:8)]), signif(out2.meta[, -(1:2)]))
	unlink(c(paste0(outfile1, ".score.*"), paste0(outfile1, ".var.*"), paste0(outfile2, ".score.*"), paste0(outfile2, ".var.*")))
})

test_that("longitudinal repeated measures gaussian", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
        pheno <- example$pheno2
        kins <- example$GRM
        obj1 <- glmmkin(y.repeated ~ sex, data = pheno, kins = kins, id = "id",random.slope = NULL, family = gaussian(link = "identity"))
	outfile1 <- tempfile()
	out1 <- SMMAT(obj1, gdsfile, group.file, meta.file.prefix = outfile1, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out1.meta <- SMMAT.meta(outfile1, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out1[, -(1:8)], digits = 5), signif(out1.meta[, -(1:2)], digits = 5))
        obj2 <- glmmkin(y.repeated ~ sex, data = pheno, kins = NULL, id = "id",random.slope = NULL, family = gaussian(link = "identity"))
	outfile2 <- tempfile()
	out2 <- SMMAT(obj2, gdsfile, group.file, meta.file.prefix = outfile2, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out2.meta <- SMMAT.meta(outfile2, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out2[, -(1:8)]), signif(out2.meta[, -(1:2)]))
	unlink(c(paste0(outfile1, ".score.*"), paste0(outfile1, ".var.*"), paste0(outfile2, ".score.*"), paste0(outfile2, ".var.*")))
})

test_that("longitudinal random time trend gaussian", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
	suppressWarnings(RNGversion("3.5.0"))
	set.seed(123)
        pheno <- example$pheno2
        kins <- example$GRM
        obj1 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = kins, id = "id",random.slope = "time", family = gaussian(link = "identity"))
	outfile1 <- tempfile()
	out1 <- SMMAT(obj1, gdsfile, group.file, meta.file.prefix = outfile1, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out1.meta <- SMMAT.meta(outfile1, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out1[, -(1:8)]), signif(out1.meta[, -(1:2)]))
        obj2 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = NULL, id = "id",random.slope = "time", family = gaussian(link = "identity"))
	outfile2 <- tempfile()
	out2 <- SMMAT(obj2, gdsfile, group.file, meta.file.prefix = outfile2, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies", tests = c("O", "E"))
	out2.meta <- SMMAT.meta(outfile2, group.file = group.file, tests = c("O", "E"))
	expect_equal(signif(out2[, -(1:8)], digits = 5), signif(out2.meta[, -(1:2)], digits = 5))
	unlink(c(paste0(outfile1, ".score.*"), paste0(outfile1, ".var.*"), paste0(outfile2, ".score.*"), paste0(outfile2, ".var.*")))
})

