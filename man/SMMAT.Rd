\name{SMMAT}
\alias{SMMAT}
\alias{SMMAT.prep}
\alias{SMMAT.lowmem}
\title{
Variant Set Mixed Model Association Tests (SMMAT)
}
\description{
Variant Set Mixed Model Association Tests (SMMAT-B, SMMAT-S, SMMAT-O and SMMAT-E) for multiple user-defined test units and a null generalized linear mixed model.
\code{SMMAT.prep} and \code{SMMAT.lowmem} are the two-step low-memory version of \code{SMMAT}. \code{SMMAT.lowmem} takes the returned R object from \code{SMMAT.prep} and uses less memory (if the returned R object from \code{SMMAT.prep} is saved to an R data file, the R session is terminated, and this R object is loaded into a new R session for running \code{SMMAT.lowmem}), especially when \code{group.file} contains only a subset of variants from \code{geno.file}.
}
\usage{
SMMAT(null.obj, geno.file, group.file, group.file.sep = "\t", 
	meta.file.prefix = NULL, MAF.range = c(1e-7, 0.5), 
	MAF.weights.beta = c(1, 25), miss.cutoff = 1, 
	missing.method = "impute2mean", method = "davies", 
	tests = "E", rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 
	0.5^2, 0.5, 1), use.minor.allele = FALSE, 
	auto.flip = FALSE, Garbage.Collection = FALSE,
	is.dosage = FALSE, ncores = 1, verbose = FALSE)
SMMAT.prep(null.obj, geno.file, group.file, group.file.sep = "\t", 
	auto.flip = FALSE)
SMMAT.lowmem(SMMAT.prep.obj, meta.file.prefix = NULL,
        MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(1, 25),
	miss.cutoff = 1, missing.method = "impute2mean",
	method = "davies", tests = "E", rho = c(0, 0.1^2,
	0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),
	use.minor.allele = FALSE, Garbage.Collection = FALSE,
	is.dosage = FALSE, ncores = 1, verbose = FALSE)
}
\arguments{
  \item{null.obj}{
a class glmmkin or a class glmmkin.multi object, returned by fitting the null GLMM using \code{glmmkin}.
}
  \item{geno.file}{
a .gds file for the full genotypes. The \code{sample.id} in \code{geno.file} should overlap \code{id_include} in \code{null.obj}. It is recommended that \code{sample.id} in \code{geno.file} include the full samples (at least all samples as specified in \code{id_include} of \code{null.obj}). It is not necessary for the user to take a subset of \code{geno.file} before running the analysis.
}
  \item{group.file}{
a plain text file with 6 columns defining the test units. There should be no headers in the file, and the columns are group name, chromosome, position, reference allele, alternative allele and weight, respectively.
}
  \item{group.file.sep}{
the delimiter in group.file (default = "\\t").
}
  \item{meta.file.prefix}{
prefix of intermediate files (.score.* and .var.*) required in a meta-analysis. If NULL, such intermediate files are not generated (default = NULL).
}
  \item{MAF.range}{
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
}
  \item{MAF.weights.beta}{
a numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies. This internal minor allele frequency weight is multiplied by the external weight given by the group.file. To turn off internal minor allele frequency weight and only use the external weight given by the group.file, use c(1, 1) to assign flat weights (default = c(1, 25)).
}
  \item{miss.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
}
  \item{missing.method}{
method of handling missing genotypes. Either "impute2mean" or "impute2zero" (default = "impute2mean").
}
  \item{method}{
a method to compute p-values for SKAT-type test statistics (default = "davies"). "davies" represents an exact method that computes a p-value by inverting the characteristic function of the mixture chisq distribution, with an accuracy of 1e-6. When "davies" p-value is less than 1e-5, it defaults to method "kuonen". "kuonen" represents a saddlepoint approximation method that computes the tail probabilities of the mixture chisq distribution. When "kuonen" fails to compute a p-value, it defaults to method "liu". "liu" is a moment-matching approximation method for the mixture chisq distribution.
}
  \item{tests}{
a character vector indicating which SMMAT tests should be performed ("B" for the burden test, "S" for SKAT, "O" for SKAT-O and "E" for the efficient hybrid test of the burden test and SKAT). The burden test and SKAT are automatically included when performing "O", and the burden test is automatically included when performing "E" (default = "E").
}
  \item{rho}{
a numeric vector defining the search grid used in SMMAT-O for SKAT-O (see the SKAT-O paper for details). Not used for SMMAT-B for the burden test, SMMAT-S for SKAT or SMMAT-E for the efficient hybrid test of the burden test and SKAT (default = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1)).
}
  \item{use.minor.allele}{
a logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SMMAT-S results, but SMMAT-B (as well as SMMAT-O and SMMAT-E) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction.
}
  \item{auto.flip}{
a logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
}
  \item{Garbage.Collection}{
a logical switch for whether to enable garbage collection in each test (default = FALSE). Pay for memory efficiency with slower computation speed.
}
  \item{is.dosage}{
a logical switch for whether imputed dosage should be used from \code{geno.file} (default = FALSE).
}
  \item{ncores}{
a positive integer indicating the number of cores to be used in parallel computing (default = 1).
}
  \item{verbose}{
a logical switch for whether a progress bar should be shown (default = FALSE).
}
  \item{SMMAT.prep.obj}{
a class SMMAT.prep object, returned by \code{SMMAT.prep}.
}
}
\value{
\code{SMMAT} and \code{SMMAT.lowmem} return a data frame with the following components:
  \item{group}{
name of the test unit group.
}
  \item{n.variants}{
number of variants in the test unit group that pass the missing rate and allele frequency filters.
}
  \item{miss.min}{
minimum missing rate for variants in the test unit group.
}
  \item{miss.mean}{
mean missing rate for variants in the test unit group.
}
  \item{miss.max}{
maximum missing rate for variants in the test unit group.
}
  \item{freq.min}{
minimum coding allele frequency for variants in the test unit group.
}
  \item{freq.mean}{
mean coding allele frequency for variants in the test unit group.
}
  \item{freq.max}{
maximum coding allele frequency for variants in the test unit group.
}
  \item{B.score}{
burden test score statistic.
}
  \item{B.var}{
variance of burden test score statistic.
}
  \item{B.pval}{
burden test p-value.
}
  \item{S.pval}{
SKAT p-value.
}
  \item{O.pval}{
SKAT-O p-value.
}
  \item{O.minp}{
minimum p-value in the SKAT-O search grid.
}
  \item{O.minp.rho}{
rho value at the minimum p-value in the SKAT-O search grid.
}
  \item{E.pval}{
SMMAT efficient hybrid test of the burden test and SKAT p-value.
}
\code{SMMAT.prep} return a list with the following components:
  \item{null.obj}{
a class glmmkin or a class glmmkin.multi object from the null model, after pre-processing.
}
  \item{geno.file}{
the name of the .gds file for the full genotypes.
}
  \item{group.file}{
the name of the plain text file with 6 columns defining the test units.
}
  \item{group.file.sep}{
the delimiter in group.file.
}
  \item{auto.flip}{
a logical indicator showing whether automatic allele flipping is enabled in pre-processing if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found.
}
  \item{residuals}{
residuals from the null model, after pre-processing.
}
  \item{sample.id}{
sample.id from geno.file, after pre-processing.
}
  \item{group.info}{
group.info read from group.file, after pre-processing.
}
  \item{groups}{
unique groups in group.info, after pre-processing.
}
  \item{group.idx.start}{
a vector of the start variant index for each group, after pre-processing.
}
  \item{group.idx.end}{
a vector of the end variant index for each group, after pre-processing.
}
}
\references{
Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., Lin, X. (2011) Rare-variant association testing for sequencing data with the sequence kernel association test. The American Journal of Human Genetics 89, 82-93.

Lee, S., Wu, M.C., Lin, X. (2012) Optimal tests for rare variant effects in sequencing association studies. Biostatistics 13, 762-775.

Sun, J., Zheng, Y., Hsu, L. (2013) A unified mixed-effects model for rare-variant association in sequencing studies. Genetic Epidemiology 37, 334-344.

Chen, H., Huffman, J.E., Brody, J.A., Wang, C., Lee, S., Li, Z., Gogarten, S.M., Sofer, T., Bielak, L.F., Bis, J.C., et al. (2019) Efficient variant set mixed model association tests for continuous and binary traits in large-scale whole-genome sequencing studies. The American Journal of Human Genetics 104, 260-274.
}
\author{
Han Chen
}
\seealso{
\code{\link{glmmkin}}, \code{\link{SMMAT.meta}}
}
\examples{
\donttest{
data(example)
attach(example)
model0 <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id",
       family = binomial(link = "logit"))
geno.file <- system.file("extdata", "geno.gds", package = "GMMAT")
group.file <- system.file("extdata", "SetID.withweights.txt",
	package = "GMMAT")
out <- SMMAT(model0, geno.file, group.file, MAF.range = c(0, 0.5), 
       		miss.cutoff = 1, method = "davies")
print(out)
}
\dontrun{
obj <- SMMAT.prep(model0, geno.file, group.file)
save(obj, file = "SMMAT.prep.tmp.Rdata")
# quit R session
# open a new R session
obj <- get(load("SMMAT.prep.tmp.Rdata"))
out <- SMMAT.lowmem(obj, MAF.range = c(0, 0.5), miss.cutoff = 1,
        method = "davies")
print(out)
unlink("SMMAT.prep.tmp.Rdata")
}
}
\keyword{ variant set-based test }
\keyword{ generalized linear mixed model }
