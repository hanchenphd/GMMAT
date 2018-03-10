\name{glmm.rvtests.meta}
\alias{glmm.rvtests.meta}
\title{
Set-based rare variant meta-analysis for generalized linear mixed models
}
\description{
Set-based rare variant tests (burden test, SKAT, SKAT-O and SMMAT) in the meta-analysis.
}
\usage{
glmm.rvtests.meta(meta.files.prefix, n.files = rep(1,
		length(meta.files.prefix)), cohort.group.idx = NULL, 
		group.file, group.file.sep = "\t", MAF.range = c(1e-7, 0.5), 
		miss.cutoff = 1, method = "davies", tests = "SMMAT", 
		rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1), 
		use.minor.allele = FALSE)
}
\arguments{
  \item{meta.files.prefix}{
a character vector for prefix of intermediate files (.score.* and .var.*) required in a meta-analysis. Each element represents the prefix of .score.* and .var.* from one cohort. The length of vector should be equal to the number of cohorts.
}
  \item{n.files}{
an integer vector with the same length as meta.files.prefix, indicating how many sets of intermediate files (.score.* and .var.*) are expected from each cohort, usually as the result of multi-threading in creating the intermediate files (default = rep(1, length(meta.files.prefix))).
}
  \item{cohort.group.idx}{
a vector with the same length as meta.files.prefix, indicating which cohorts should be grouped together in the meta-analysis assuming homogeneous genetic effects. For example, c("a","b","a","a","b") means cohorts 1, 3, 4 are assumed to have homogeneous genetic effects, and cohorts 2, 5 are in another group with homogeneous genetic effects (but possibly heterogeneous with group "a"). If NULL, all cohorts are in the same group (default = NULL).
}

  \item{group.file}{
a plain text file with 6 columns defining the test units. There should be no headers in the file, and the columns are group name, chromosome, position, reference allele, alternative allele and weight, respectively.
}
  \item{group.file.sep}{
the delimiter in group.file (default = "\\t").
}
  \item{MAF.range}{
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)). Filter applied to the combined samples.}
  \item{miss.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants). Filter applied to the combined samples.
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
a logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SKAT results, but burden test (as well as hybrid tests such as SKAT-O and SMMAT) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction. Use with caution, as major/minor alleles may flip in different cohorts. In that case, minor allele will be determined based on the allele frequency in the combined samples.
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
\code{\link{glmmkin}}, \code{\link{glmm.rvtests}}
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
       		meta.file.prefix = "glmm.rvtests", MAF.range = c(0, 0.5), 
		miss.cutoff = 1, method = "davies")
print(out)
out1 <- glmm.rvtests.meta("glmm.rvtests", group.file = group.file)
print(out1)
}
\keyword{ variant set-based test }
\keyword{ generalized linear mixed model }
