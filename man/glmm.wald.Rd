\name{glmm.wald}
\alias{glmm.wald}
\title{
Performing GLMM based Wald tests
}
\description{
Fit a GLMM under the alternative hypothesis to perform Wald tests for association with genotypes in a plink .bed file (binary genotypes), a GDS file .gds, or a plain text file (or compressed .gz or .bz2 file).
}
\usage{
glmm.wald(fixed, data = parent.frame(), kins = NULL, id, random.slope = NULL, 
	groups = NULL, family = binomial(link = "logit"), infile, snps, 
	method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, 
	taumin = 1e-5, taumax = 1e5, tauregion = 10, center = T, 
	select = NULL, missing.method = "impute2mean", infile.nrow = NULL, 
	infile.nrow.skip = 0, infile.sep = "\t", infile.na = "NA", 
	snp.col = 1, infile.ncol.skip = 1, infile.ncol.print = 1, 
	infile.header.print = "SNP", is.dosage = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{fixed}{
an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the fixed effects model to be fitted.
}
  \item{data}{
a data frame or list (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model.
}
  \item{kins}{
a known positive semi-definite relationship matrix (e.g. kinship matrix in genetic association studies) or a list of known positive semi-definite relationship matrices. The rownames and colnames of these matrices must at least include all samples as specified in the \code{id} column of the data frame \code{data}. If not provided, glmmkin will switch to the generalized linear model with no random effects (default = NULL).
}
  \item{id}{
a column in the data frame \code{data}, indicating the id of samples. When there are duplicates in id, the data is assumed to be longitudinal with repeated measures.
}
  \item{random.slope}{
an optional column indicating the random slope for time effect used in a mixed effects model for longitudinal data. It must be included in the names of \code{data}. There must be duplicates in \code{id} and \code{method.optim} must be "AI" (default = NULL).
}
\item{groups}{
an optional categorical variable indicating the groups used in a heteroscedastic linear mixed model (allowing residual variances in different groups to be different). This variable must be included in the names of \code{data}, and \code{family} must be "gaussian" and \code{method.optim} must be "AI" (default = NULL).
}
\item{family}{
a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
}
  \item{infile}{
the input file name. Note that for plink binary genotype files only the prefix without .bed, .bim or .fam should be used. Only SNP major mode recognized in the binary file. Alternatively, it can be the full name of a GDS file (including the suffix .gds) or a plain text file with some delimiters (comma, space, tab or something else), with one row for each SNP and one column for each individual. In that case, SNPs should be coded as numeric values (0/1/2 or dosages allowed, A/C/G/T coding is not recognized). There can be additional rows and columns to skip at the beginning. The order of individuals can be different from \code{obj} in the null GLMM (see the argument \code{select}). Some compressed files (.gz and .bz2) also allowed.
}
  \item{snps}{
a vector of SNP names to be tested.
}
  \item{method}{
method of fitting the generalized linear mixed model. Either "REML" or "ML" (default = "REML").
}
  \item{method.optim}{
optimization method of fitting the generalized linear mixed model. Either "AI", "Brent" or "Nelder-Mead" (default = "AI").
}
  \item{maxiter}{
a positive integer specifying the maximum number of iterations when fitting the generalized linear mixed model (default = 500).
}
  \item{tol}{
a positive number specifying tolerance, the difference threshold for parameter estimates below which iterations should be stopped. Also the threshold for determining monomorphism. If a SNP has value range less than the tolerance, it will be considered monomorphic and its association test p-value will be NA (default = 1e-5).
}
  \item{taumin}{
the lower bound of search space for the variance component parameter \eqn{\tau} (default = 1e-5), used when \code{method.optim = "Brent"}. See \code{\link{glmmkin}}.
}
  \item{taumax}{
the upper bound of search space for the variance component parameter \eqn{\tau} (default = 1e5), used when \code{method.optim = "Brent"}. See \code{\link{glmmkin}}.
}
  \item{tauregion}{
the number of search intervals for the REML or ML estimate of the variance component parameter \eqn{\tau} (default = 10), used when \code{method.optim = "Brent"}. See \code{\link{glmmkin}}.
}
  \item{center}{
a logical switch for centering genotypes before tests. If TRUE, genotypes will be centered to have mean 0 before tests, otherwise raw values will be directly used in tests (default = TRUE).
}
  \item{select}{
an optional vector indicating the order of individuals in \code{infile}. If supplied, the length must match the number of individuals in \code{infile} (default = NULL). Individuals to be excluded should be coded 0. For example, \code{select = c(2, 3, 1, 0)} means the 1st individual in \code{infile} corresponds to the 2nd individual in \code{data}, the 2nd individual in \code{infile} corresponds to the 3rd individual in \code{data}, the 3rd individual in \code{infile} corresponds to the 1st individual in \code{data}, the 4th individual in \code{infile} is not included in \code{data}. If there are any duplicated \code{id} in \code{data} (longitudinal data analysis), indices in \code{select} should match the order of individuals with unique \code{id} in \code{data}. For plink binary genotype files and GDS files, this argument is not required and the sample ID's are automatically matched.
}
  \item{missing.method}{
method of handling missing genotypes. Either "impute2mean" or "omit" (default = "impute2mean").
}
  \item{infile.nrow}{
number of rows to read in \code{infile}, including number of rows to skip at the beginning. If NULL, the program will determine how many rows there are in \code{infile} automatically and read all rows (default = NULL). Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{infile.nrow.skip}{
number of rows to skip at the beginning of \code{infile}. Must be nonnegative integers. Useful when header or comment lines are present (default = 0). Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{infile.sep}{
delimiter in \code{infile} (default = "\\t"). Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{infile.na}{
symbol in \code{infile} to denote missing genotypes (default = "NA"). Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{snp.col}{
a positive integer specifying which column in \code{infile} is SNP names. Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{infile.ncol.skip}{
number of columns to skip before genotype data in \code{infile}. These columns can be SNP name, alleles and/or quality measures and should be placed at the beginning in each line. After skipping these columns, the program will read in genotype data and perform Wald tests. Must be positive integers. It is recommended that SNP name should be included as the first column in \code{infile} and genotype data should start from the second column or later (default = 1). Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{infile.ncol.print}{
a vector indicating which column(s) in \code{infile} should be shown in the results. These columns can be SNP name, alleles and/or quality measures placed at the beginning in each line. Must be positive integers, no greater than \code{infile.ncol.skip} and sorted numerically in ascending order. By default, it is assumed that the first column is SNP name and genotype data start from the second column, and SNP name should be carried over to the results (default = 1). Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{infile.header.print}{
a character vector indicating column name(s) of column(s) selected to print by \code{infile.ncol.print} (default = "SNP"). Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{is.dosage}{
a logical switch for whether imputed dosage should be used from a GDS \code{infile} (default = FALSE).
}
  \item{verbose}{
a logical switch for printing a progress bar and detailed information (parameter estimates in each iteration) for testing and debugging purpose (default = FALSE).
}
  \item{\dots}{
additional arguments that could be passed to \code{\link{glm}}.
}
}
\value{
if \code{infile} is a plain text file, a data frame containing variables included in \code{infile.header.print} and the following:
  \item{N}{number of individuals with non-missing genotypes for each SNP.}
  \item{AF}{effect allele frequency for each SNP.}
  \item{BETA}{effect size estimate for each SNP from the GLMM under the alternative hypothesis.}
  \item{SE}{standard error of the effect size estimate for each SNP.}
  \item{PVAL}{Wald test p-value for each SNP.}
  \item{converged}{a logical indicator for convergence for each SNP.}
if \code{infile} is the prefix of plink binary files (.bed, .bim and .fam), a data frame containing the following:
  \item{CHR}{Chromosome, copied from .bim file.}
  \item{SNP}{SNP name, as supplied in \code{snps}.}
  \item{cM}{genetic location in centi Morgans, copied from .bim file.}
  \item{POS}{physical position in base pairs, copied from .bim file.}
  \item{A1}{allele 1, copied from .bim file.}
  \item{A2}{allele 2, copied from .bim file.}
  \item{N}{number of individuals with non-missing genotypes for each SNP.}
  \item{AF}{effect allele frequency for each SNP.}
  \item{BETA}{effect size estimate for each SNP from the GLMM under the alternative hypothesis.}
  \item{SE}{standard error of the effect size estimate for each SNP.}
  \item{PVAL}{Wald test p-value for each SNP.}
  \item{converged}{a logical indicator for convergence for each SNP.}
if \code{infile} is a GDS file (.gds), a data frame containing the following:
  \item{SNP}{SNP name, as supplied in \code{snps}.}
  \item{CHR}{Chromosome, copied from .gds file.}
  \item{POS}{physical position in base pairs, copied from .gds file.}
  \item{REF}{reference allele, copied from .gds file.}
  \item{ALT}{alternate allele, copied from .gds file.}
  \item{N}{number of individuals with non-missing genotypes for each SNP.}
  \item{AF}{ALT allele frequency for each SNP.}
  \item{BETA}{effect size estimate for each SNP from the GLMM under the alternative hypothesis.}
  \item{SE}{standard error of the effect size estimate for each SNP.}
  \item{PVAL}{Wald test p-value for each SNP.}
  \item{converged}{a logical indicator for convergence for each SNP.}
}
\references{
Brent, R.P. (1973) "Chapter 4: An Algorithm with Guaranteed Convergence for Finding a Zero of a Function", Algorithms for Minimization without Derivatives, Englewood Cliffs, NJ: Prentice-Hall, ISBN 0-13-022335-2.

Breslow, N.E. and Clayton, D.G. (1993) Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association 88, 9-25.

Chen, H., Wang, C., Conomos, M.P., Stilp, A.M., Li, Z., Sofer, T., Szpiro, A.A., Chen, W., Brehm, J.M., Celedón, J.C., Redline, S., Papanicolaou, G.J., Thornton, T.A., Laurie, C.C., Rice, K. and Lin, X. (2016) Control forpopulation structure and relatedness for binary traits in genetic association studies via logistic mixed models. The American Journal of Human Genetics 98, 653-666.

Gilmour, A.R., Thompson, R. and Cullis, B.R. (1995) Average Information REML: An Efficient Algorithm for Variance Parameter Estimation in Linear Mixed Models. Biometrics 51, 1440-1450.

Nelder, J.A. and Mead, R. (1965) A simplex algorithm for function minimization. Computer Journal 7, 308-313.

Yang, J., Lee, S.H., Goddard, M.E. and Visscher, P.M. (2011) GCTA: A Tool for Genome-wide Complex Trait Analysis. The American Journal of Human Genetics 88, 76-82.

Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis for association studies. Nature Genetics 44, 821-824.
}
\author{
Han Chen, Matthew P. Conomos
}
\seealso{
\code{\link{glmmkin}}, \code{\link{glmm.score}}
}
\examples{
\donttest{
data(example)
attach(example)
snps <- c("SNP10", "SNP25", "SNP1", "SNP0")
plinkfiles <- strsplit(system.file("extdata", "geno.bed", package = "GMMAT"), 
       ".bed", fixed = TRUE)[[1]]
glmm.wald(disease ~ age + sex, data = pheno, kins = GRM, id = "id",
	family = binomial(link = "logit"), infile = plinkfiles, snps = snps) 
infile <- system.file("extdata", "geno.gds", package = "GMMAT")
glmm.wald(disease ~ age + sex, data = pheno, kins = GRM, id = "id",
	family = binomial(link = "logit"), infile = infile, snps = snps) 
infile <- system.file("extdata", "geno.txt", package = "GMMAT")
glmm.wald(disease ~ age + sex, data = pheno, kins = GRM, id = "id",
	family = binomial(link = "logit"), infile = infile, snps = snps, 
	infile.nrow.skip = 5, infile.ncol.skip = 3, infile.ncol.print = 1:3, 
	infile.header.print = c("SNP", "Allele1", "Allele2"))
}
}
\keyword{ generalized linear mixed model }
\keyword{ Wald test }
