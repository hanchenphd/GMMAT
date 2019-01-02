context("variant set mixed model association test (SMMAT)")

test_that("cross-sectional id le 400 binomial", {
	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
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
	out1 <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out1$B.pval)), signif(c(0.1540811, 0.9500619)))
	expect_equal(signif(range(out1$E.pval)), signif(c(0.01675254, 0.76516832)))
	obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	out2 <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out2$B.pval)), signif(c(0.1401707, 0.9730224)))
	expect_equal(signif(range(out2$E.pval)), signif(c(0.01782418, 0.73027324)))

	skip_on_cran()

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
	obj1 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
	obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
	obj1 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
	obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)
})

test_that("cross-sectional id gt 400 binomial", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
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
	out1 <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out1$B.pval)), signif(c(0.7906281, 0.9818689)))
	expect_equal(signif(range(out1$E.pval), digits = 5), signif(c(0.0585866, 0.9423505), digits = 5))
	obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	out2 <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out2$B.pval)), signif(c(0.7086258, 0.9753463)))
	expect_equal(signif(range(out2$E.pval), digits = 5), signif(c(0.05865188, 0.96100550), digits = 5))

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
	obj1 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
	obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
	obj1 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
	obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)
})

test_that("cross-sectional id le 400 gaussian", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
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
	out1 <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out1$B.pval)), signif(c(0.1515303, 0.9184369)))
	expect_equal(signif(range(out1$E.pval), digits = 5), signif(c(0.1199585, 0.9946534), digits = 5))
	obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	out2 <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out2$B.pval)), signif(c(0.06160726, 0.95356346)))
	expect_equal(signif(range(out2$E.pval)), signif(c(0.04304936, 0.70569494)))

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
	obj1 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
	obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
	obj1 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
	obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)
})

test_that("cross-sectional id gt 400 gaussian", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
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
	out1 <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out1$B.pval)), signif(c(0.1882496, 0.6073079)))
	expect_equal(signif(range(out1$E.pval)), signif(c(0.08166496, 0.87474253)))
	obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	out2 <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out2$B.pval)), signif(c(0.02197351, 0.04946657)))
	expect_equal(signif(range(out2$E.pval)), signif(c(0.004167252, 0.151422542)))

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
	obj1 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
	obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
	obj1 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
	obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"), method = "REML", method.optim = "AI")
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)
})

test_that("longitudinal repeated measures gaussian", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
	set.seed(123)
        pheno <- example$pheno2
        kins <- example$GRM
        obj1 <- glmmkin(y.repeated ~ sex, data = pheno, kins = kins, id = "id",random.slope = NULL, family = gaussian(link = "identity"))
	out1 <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out1$B.pval)), signif(c(0.06520763, 0.76178619)))
	expect_equal(signif(range(out1$E.pval)), signif(c(0.08942398, 0.90207219)))
        obj2 <- glmmkin(y.repeated ~ sex, data = pheno, kins = NULL, id = "id",random.slope = NULL, family = gaussian(link = "identity"))
	out2 <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out2$B.pval)), signif(c(0.02074813, 0.89820104)))
	expect_equal(signif(range(out2$E.pval)), signif(c(0.03735127, 0.97354118)))

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
        obj1 <- glmmkin(y.repeated ~ sex, data = pheno, kins = kins, id = "id",random.slope = NULL, family = gaussian(link = "identity"))
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
        obj2 <- glmmkin(y.repeated ~ sex, data = pheno, kins = NULL, id = "id",random.slope = NULL, family = gaussian(link = "identity"))
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
        obj1 <- glmmkin(y.repeated ~ sex, data = pheno, kins = kins, id = "id",random.slope = NULL, family = gaussian(link = "identity"))
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
        obj2 <- glmmkin(y.repeated ~ sex, data = pheno, kins = NULL, id = "id",random.slope = NULL, family = gaussian(link = "identity"))
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)
})

test_that("longitudinal random time trend gaussian", {
	skip_on_cran()

	gdsfile <- system.file("extdata", "geno.gds", package = "GMMAT")
	group.file <- system.file("extdata", "SetID.withweights.txt", package = "GMMAT")
	data(example)
	set.seed(123)
        pheno <- example$pheno2
        kins <- example$GRM
        obj1 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = kins, id = "id",random.slope = "time", family = gaussian(link = "identity"))
	out1 <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out1$B.pval)), signif(c(0.1333868, 0.7679753)))
	expect_equal(signif(range(out1$E.pval)), signif(c(0.2490228, 0.8947057)))
        obj2 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = NULL, id = "id",random.slope = "time", family = gaussian(link = "identity"))
	out2 <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(signif(range(out2$B.pval)), signif(c(0.08870332, 0.94362906)))
	expect_equal(signif(range(out2$E.pval)), signif(c(0.1845990, 0.9835685)))

	idx <- sample(nrow(pheno))
	pheno <- pheno[idx, ]
        obj1 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = kins, id = "id",random.slope = "time", family = gaussian(link = "identity"))
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
        obj2 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = NULL, id = "id",random.slope = "time", family = gaussian(link = "identity"))
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)

	idx <- sample(nrow(kins))
	kins <- kins[idx, idx]
        obj1 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = kins, id = "id",random.slope = "time", family = gaussian(link = "identity"))
	tmpout <- SMMAT(obj1, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out1, tmpout)
        obj2 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = NULL, id = "id",random.slope = "time", family = gaussian(link = "identity"))
	tmpout <- SMMAT(obj2, gdsfile, group.file, MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
	expect_equal(out2, tmpout)
})

