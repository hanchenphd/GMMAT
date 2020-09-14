\name{SMMAT.meta}
\alias{SMMAT.meta}
\title{
Meta-analysis for variant Set Mixed Model Association Tests (SMMAT)
}
\description{
Variant Set Mixed Model Association Tests (SMMAT-B, SMMAT-S, SMMAT-O and SMMAT-E) in the meta-analysis.
}
\usage{
SMMAT.meta(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), 
	cohort.group.idx = NULL, group.file, group.file.sep = "\t", 
	MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(1, 25), 
	miss.cutoff = 1, method = "davies", tests = "E", rho = c(0, 0.1^2, 
	0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1), use.minor.allele = FALSE,
	verbose = FALSE)
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
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)). Filter applied to the combined samples.
}
  \item{MAF.weights.beta}{
a numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies. This internal minor allele frequency weight is multiplied by the external weight given by the group.file. To turn off internal minor allele frequency weight and only use the external weight given by the group.file, use c(1, 1) to assign flat weights (default = c(1, 25)). Applied to the combined samples.
}
  \item{miss.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants). Filter applied to the combined samples.
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
a logical switch for whether to use the minor allele (instead of the alt allele) as the coding allele (default = FALSE). It does not change SMMAT-S results, but SMMAT-B (as well as SMMAT-O and SMMAT-E) will be affected. Along with the MAF filter, this option is useful for combining rare mutations, assuming rare allele effects are in the same direction. Use with caution, as major/minor alleles may flip in different cohorts. In that case, minor allele will be determined based on the allele frequency in the combined samples.
}
  \item{verbose}{
a logical switch for whether a progress bar should be shown (default = FALSE).
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
\code{\link{glmmkin}}, \code{\link{SMMAT}}
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
out <- SMMAT(model0, geno.file, group.file, meta.file.prefix = "SMMAT.meta", 
       		MAF.range = c(0, 0.5), miss.cutoff = 1, method = "davies")
print(out)
out1 <- SMMAT.meta("SMMAT.meta", group.file = group.file)
print(out1)
unlink(paste0("SMMAT.meta.", c("score", "var"), ".1"))
}
}
\keyword{ variant set-based test }
\keyword{ generalized linear mixed model }
