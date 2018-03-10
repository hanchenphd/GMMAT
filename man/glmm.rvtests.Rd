\name{glmm.rvtests}
\alias{glmm.rvtests}
\title{
Set-based rare variant tests for generalized linear mixed models
}
\description{
Set-based rare variant tests (burden test, SKAT, SKAT-O and SMMAT) for multiple user-defined test units and a null generalized linear mixed model.
}
\usage{
glmm.rvtests(null.obj, null.obj.id = NULL, geno.file, group.file,
		group.file.sep = "\t", meta.file.prefix = NULL,
		MAF.range = c(1e-7, 0.5), miss.cutoff = 1, 
		missing.method = "impute2mean", method = "davies", 
		tests = "SMMAT", rho = c(0, 0.1^2, 0.2^2, 0.3^2, 
		0.4^2, 0.5^2, 0.5, 1), use.minor.allele = FALSE, 
		auto.flip = FALSE, Garbage.Collection = FALSE, 
		ncores = 1)
}
\arguments{
  \item{null.obj}{
a class glmmkin object, returned by fitting the null GLMM using \code{glmmkin( )}.
}
  \item{null.obj.id}{
a vector with length equal to the sample size in null.obj indicating the sample ids in null.obj (default = NULL, then using null.obj$id_include). The ids should be a subset of sample ids in geno.file.
}
  \item{geno.file}{
a .gds file for the full genotypes. The associated sample ids should match the id format in null.obj.id, but sorting is not necessary and it can include more samples than in null.obj.id.
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
a character vector indicating which tests should be performed ("Burden", "SKAT", "SKAT-O" and "SMMAT"). Burden test is automatically included when performing "SMMAT" (default = "SMMAT").
}
  \item{rho}{
a numeric vector defining the search grid used in SKAT-O (see the SKAT-O paper for details). Not used for Burden, SKAT or SMMAT (default = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1)).
}
  \item{use.minor.allele}{
a logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SKAT results, but burden test (as well as hybrid tests such as SKAT-O and SMMAT) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction.
}
  \item{auto.flip}{
a logical switch for whether to enable automatic allele flipping if a variant with alleles ref/alt is not found at a position, but a variant at the same position with alleles alt/ref is found (default = FALSE). Use with caution for whole genome sequence data, as both ref/alt and alt/ref variants at the same position are not uncommon, and they are likely two different variants, rather than allele flipping.
}
  \item{Garbage.Collection}{
a logical switch for whether to enable garbage collection in each test (default = FALSE). Pay for memory efficiency with slower computation speed.
}
  \item{ncores}{
a positive integer indicating the number of cores to be used in parallel
computing (default = 1).
}
}
\value{
a data frame with the following components:
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
  \item{Burden.score}{
burden test score statistic.
}
  \item{Burden.var}{
variance of burden test score statistic.
}
  \item{Burden.pval}{
burden test p-value.
}
  \item{SKAT.pval}{
SKAT p-value.
}
  \item{SKATO.pval}{
SKAT-O p-value.
}
  \item{SKATO.minp}{
minimum p-value in the SKAT-O search grid.
}
  \item{SKATO.minp.rho}{
rho value at the minimum p-value in the SKAT-O search grid.
}
  \item{SMMAT.pval}{
SMMAT p-value.
}
}
\references{
Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., Lin, X. (2011) Rare-variant association testing for sequencing data with the sequence kernel association test. The American Journal of Human Genetics 89, 82-93.

Lee, S., Wu, M.C., Lin, X. (2012) Optimal tests for rare variant effects in sequencing association studies. Biostatistics 13, 762-775.

Sun, J., Zheng, Y., Hsu, L. (2013) A unified mixed-effects model for rare-variant association in sequencing studies. Genetic Epidemiology 37, 334-344.
}
\author{
Han Chen
}
\seealso{
\code{\link{glmmkin}}, \code{\link{glmm.rvtests.meta}}
}
\examples{
data(example)
attach(example)
model0 <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM,
       family = binomial(link = "logit"))
geno.file <- system.file("extdata", "geno.gds", package = "GMMAT")
group.file <- system.file("extdata", "SetID.withweights.txt",
	package = "GMMAT")
out <- glmm.rvtests(model0, null.obj.id = NULL, geno.file, group.file,
       		MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
print(out)
}
\keyword{ variant set-based test }
\keyword{ generalized linear mixed model }
