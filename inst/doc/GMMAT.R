### R code from vignette source 'GMMAT.Rnw'

###################################################
### code chunk number 1: installation (eval = FALSE)
###################################################
## ## try http:// if https:// URLs are not supported
## ## remove "doMC" below if you are running Windows
## install.packages(c("devtools", "RcppArmadillo", "CompQuadForm", "doMC", 
##         "foreach", "Matrix", "BiocManager", "testthat"), 
## 	repos = "https://cran.r-project.org/")
## BiocManager::install(c("SeqArray", "SeqVarTools"))
## devtools::install_github("hanchenphd/GMMAT")


###################################################
### code chunk number 2: pheno
###################################################
pheno.file <- system.file("extdata", "pheno.txt", package = "GMMAT")
pheno <- read.table(pheno.file, header = TRUE)


###################################################
### code chunk number 3: pheno2 (eval = FALSE)
###################################################
## pheno <- read.table(pheno.file, header = TRUE, na.strings = ".")


###################################################
### code chunk number 4: GRM
###################################################
GRM.file <- system.file("extdata", "GRM.txt.bz2", package = "GMMAT")
GRM <- as.matrix(read.table(GRM.file, check.names = FALSE))


###################################################
### code chunk number 5: Mats (eval = FALSE)
###################################################
## Mats <- list(Mat1, Mat2, Mat3)


###################################################
### code chunk number 6: convert2GDS (eval = FALSE)
###################################################
## SeqArray::seqVCF2GDS("VCF_file_name", "GDS_file_name")
## SeqArray::seqBED2GDS("BED_file_name", "FAM_file_name", "BIM_file_name", 
##         "GDS_file_name")


###################################################
### code chunk number 7: loading
###################################################
library(GMMAT)


###################################################
### code chunk number 8: help (eval = FALSE)
###################################################
## ?glmmkin


###################################################
### code chunk number 9: GMMAT0
###################################################
model0 <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, 
        id = "id", family = binomial(link = "logit"))
model0$theta
model0$coefficients
model0$cov


###################################################
### code chunk number 10: GMMAT01 (eval = FALSE)
###################################################
## model0 <- glmmkin(fixed = disease ~ age + sex, data = pheno, kins = GRM, 
##         id = "id", family = binomial(link = "logit"))


###################################################
### code chunk number 11: GMMAT1
###################################################
model1 <- glmmkin(fixed = trait ~ age + sex, data = pheno, kins = GRM, 
        id = "id", family = gaussian(link = "identity"))


###################################################
### code chunk number 12: GMMAT2
###################################################
model2 <- glmmkin(fixed = trait ~ age + sex, data = pheno, kins = GRM, 
        id = "id", groups = "disease", 
        family = gaussian(link = "identity"))
model2$theta


###################################################
### code chunk number 13: GMMAT3
###################################################
M10 <- matrix(0, 400, 400)
for(i in 1:40) M10[(i-1)*10+(1:10), (i-1)*10+(1:10)] <- 1
rownames(M10) <- colnames(M10) <- 1:400
Mats <- list(GRM, M10)
model3 <- glmmkin(fixed = disease ~ age + sex, data = pheno, id = "id",
        kins = Mats, family = binomial(link = "logit"))
model3$theta


###################################################
### code chunk number 14: GMMAT4
###################################################
pheno2.file <- system.file("extdata", "pheno2.txt", package = "GMMAT")
pheno2 <- read.table(pheno2.file, header = TRUE)
model4 <- glmmkin(y.repeated ~ sex, data = pheno2, kins = GRM, id = "id",
        family = gaussian(link = "identity"))
model4$theta


###################################################
### code chunk number 15: GMMAT5
###################################################
model5 <- glmmkin(y.trend ~ sex + time, data = pheno2, kins = GRM, id = "id",
        random.slope = "time", family = gaussian(link = "identity"))
model5$theta


###################################################
### code chunk number 16: GMMATscoretxt
###################################################
geno.file <- system.file("extdata", "geno.txt", package = "GMMAT")
glmm.score(model0, infile = geno.file, outfile =
        "glmm.score.text.testoutfile.txt", infile.nrow.skip = 5, 
        infile.ncol.skip = 3, infile.ncol.print = 1:3,
        infile.header.print = c("SNP", "Allele1", "Allele2"))


###################################################
### code chunk number 17: GMMATscorebed
###################################################
geno.file <- strsplit(system.file("extdata", "geno.bed", 
        package = "GMMAT"), ".bed", fixed = TRUE)[[1]]
glmm.score(model0, infile = geno.file, outfile = 
        "glmm.score.bed.testoutfile.txt")


###################################################
### code chunk number 18: GMMATscoregds
###################################################
geno.file <- system.file("extdata", "geno.gds", package = "GMMAT")
glmm.score(model0, infile = geno.file, outfile = 
        "glmm.score.gds.testoutfile.txt")


###################################################
### code chunk number 19: GMMATscorebgen
###################################################
geno.file <- system.file("extdata", "geno.bgen", package = "GMMAT")
samplefile <- system.file("extdata", "geno.sample", package = "GMMAT")
glmm.score(model0, infile = geno.file, BGEN.samplefile = samplefile,
	outfile = "glmm.score.bgen.testoutfile.txt")


###################################################
### code chunk number 20: GMMATwaldtxt
###################################################
geno.file <- system.file("extdata", "geno.txt", package = "GMMAT")
snps <- c("SNP10", "SNP25", "SNP1", "SNP0")
glmm.wald(fixed = disease ~ age + sex, data = pheno, kins = GRM, id = "id",
        family = binomial(link = "logit"), infile = geno.file, snps = snps, 
	infile.nrow.skip = 5, infile.ncol.skip = 3, infile.ncol.print = 1:3, 
	infile.header.print = c("SNP", "Allele1", "Allele2"))


###################################################
### code chunk number 21: GMMATwaldbed
###################################################
geno.file <- strsplit(system.file("extdata", "geno.bed", 
        package = "GMMAT"), ".bed", fixed = TRUE)[[1]]
glmm.wald(fixed = disease ~ age + sex, data = pheno, kins = GRM, id = "id",
        family = binomial(link = "logit"), infile = geno.file, snps = snps)


###################################################
### code chunk number 22: GMMATwaldgds
###################################################
geno.file <- system.file("extdata", "geno.gds", package = "GMMAT")
glmm.wald(fixed = disease ~ age + sex, data = pheno, kins = GRM, id = "id",
        family = binomial(link = "logit"), infile = geno.file, snps = snps)


###################################################
### code chunk number 23: GMMATscoremeta
###################################################
meta1.file <- system.file("extdata", "meta1.txt", package = "GMMAT")
meta2.file <- system.file("extdata", "meta2.txt", package = "GMMAT")
meta3.file <- system.file("extdata", "meta3.txt", package = "GMMAT")
glmm.score.meta(files = c(meta1.file, meta2.file, meta3.file),
        outfile = "glmm.score.meta.testoutfile.txt",
        SNP = rep("SNP", 3), A1 = rep("A1", 3), A2 = rep("A2", 3))


###################################################
### code chunk number 24: SMMAT1
###################################################
group.file <- system.file("extdata", "SetID.withweights.txt", 
        package = "GMMAT")
geno.file <- system.file("extdata", "geno.gds", package = "GMMAT")
SMMAT(model0, group.file = group.file, geno.file = geno.file, 
        MAF.range = c(1e-7, 0.5), miss.cutoff = 1, method = "davies", 
        tests = c("O", "E"))


###################################################
### code chunk number 25: SMMAT2
###################################################
SMMAT(model0, group.file = group.file, geno.file = geno.file, 
        MAF.range = c(1e-7, 0.5), miss.cutoff = 1, method = "davies", 
        tests = "B", meta.file.prefix = "SMMAT.meta")


###################################################
### code chunk number 26: SMMATmeta
###################################################
SMMAT.meta(meta.files.prefix = "SMMAT.meta", n.files = 1,
        group.file = group.file, MAF.range = c(1e-7, 0.5), 
        miss.cutoff = 1, method = "davies", tests = "S")


###################################################
### code chunk number 27: MKL (eval = FALSE)
###################################################
## Sys.setenv(MKL_NUM_THREADS = 1)


###################################################
### code chunk number 28: select (eval = FALSE)
###################################################
## select <- match(geno_ID, pheno_ID[obj$id_include])
## select[is.na(select)] <- 0


###################################################
### code chunk number 29: select2 (eval = FALSE)
###################################################
## select <- match(geno_ID, unique(obj$id_include))
## select[is.na(select)] <- 0


###################################################
### code chunk number 30: select3 (eval = FALSE)
###################################################
## select <- match(geno_ID, unique(data[, id]))
## select[is.na(select)] <- 0


