\name{glmm.score}
\alias{glmm.score}
\title{
Performing GLMM based score tests
}
\description{
Use a glmmkin class object from the null GLMM to perform score tests for association with genotypes in a plink .bed file (binary genotypes), a GDS file .gds, or a plain text file (or compressed .gz or .bz2 file).
}
\usage{
glmm.score(obj, infile, outfile, BGEN.samplefile = NULL, center = T, select = NULL, 
	MAF.range = c(1e-7, 0.5), miss.cutoff = 1, 
	missing.method = "impute2mean", nperbatch = 100, tol = 1e-5,
	infile.nrow = NULL, infile.nrow.skip = 0, infile.sep = "\t",
	infile.na = "NA", infile.ncol.skip = 1, infile.ncol.print = 1,
	infile.header.print = "SNP", is.dosage = FALSE, ncores = 1, verbose = FALSE)
}
\arguments{
  \item{obj}{
a class glmmkin or class glmmkin.multi object, returned by fitting the null GLMM using \code{glmmkin}.
}
  \item{infile}{
the input file name. Note that for plink binary genotype files only the prefix without .bed, .bim or .fam should be used. Only SNP major mode recognized in the binary file. Alternatively, it can be the full name of a BGEN file (including the suffix .bgen), a GDS file (including the suffix .gds) or a plain text file with some delimiters (comma, space, tab or something else), with one row for each SNP and one column for each individual. In that case, SNPs should be coded as numeric values (0/1/2 or dosages allowed, A/C/G/T coding is not recognized). There can be additional rows and columns to skip at the beginning. The order of individuals can be different from \code{obj} in the null GLMM (see the argument \code{select}). Some compressed files (.gz and .bz2) also allowed.
}
  \item{outfile}{
the output file name.
}
  \item{BGEN.samplefile}{
path to the BGEN sample file. Required when the BGEN file does not contain sample identifiers or the \code{select} parameter is NULL (default = NULL).
}
  \item{center}{
a logical switch for centering genotypes before tests. If TRUE, genotypes will be centered to have mean 0 before tests, otherwise raw values will be directly used in tests (default = TRUE).
}
  \item{select}{
an optional vector indicating the order of individuals in \code{infile}. If supplied, the length must match the number of individuals in \code{infile} (default = NULL). Individuals to be excluded should be coded 0. For example, \code{select = c(2, 3, 1, 0)} means the 1st individual in \code{infile} corresponds to the 2nd individual in \code{obj}, the 2nd individual in \code{infile} corresponds to the 3rd individual in \code{obj}, the 3rd individual in \code{infile} corresponds to the 1st individual in \code{obj}, the 4th individual in \code{infile} is not included in \code{obj}. If there are any duplicated \code{id_include} in \code{obj} (longitudinal data analysis), indices in \code{select} should match the order of individuals with unique \code{id_include} in \code{obj}. For plink binary genotype files and GDS files, this argument is not required and the sample ID's are automatically matched.
}
  \item{MAF.range}{
a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(1e-7, 0.5)).
}
  \item{miss.cutoff}{
the maximum missing rate allowed for a variant to be included (default = 1, including all variants).
}
  \item{missing.method}{
method of handling missing genotypes. Either "impute2mean" or "omit" (default = "impute2mean").
}
  \item{nperbatch}{
an integer for how many SNPs should be tested in a batch (default = 100). The computational time can increase dramatically if this value is either small or large. The optimal value for best performance depends on the user's system.
}
  \item{tol}{
the threshold for determining monomorphism. If a SNP has value range less than the tolerance, it will be considered monomorphic and its association test p-value will be NA (default = 1e-5). Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
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
  \item{infile.ncol.skip}{
number of columns to skip before genotype data in \code{infile}. These columns can be SNP name, alleles and/or quality measures and should be placed at the beginning in each line. After skipping these columns, the program will read in genotype data and perform score tests. Must be nonnegative integers. It is recommended that SNP name should be included as the first column in \code{infile} and genotype data should start from the second column or later (default = 1). Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{infile.ncol.print}{
a vector indicating which column(s) in \code{infile} should be printed to the output directly. These columns can be SNP name, alleles and/or quality measures placed at the beginning in each line. Must be nonnegative integers, no greater than \code{infile.ncol.skip} and sorted numerically in ascending order. By default, it is assumed that the first column is SNP name and genotype data start from the second column, and SNP name should be carried over to the output (default = 1). Should be set to NULL if \code{infile.ncol.skip} is 0. Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{infile.header.print}{
a character vector indicating column name(s) of column(s) selected to print by \code{infile.ncol.print} (default = "SNP"). Should be set to NULL if \code{infile.ncol.skip} is 0. Only used when \code{infile} is a plain text file (or compressed .gz or .bz2 file).
}
  \item{is.dosage}{
a logical switch for whether imputed dosage should be used from a GDS \code{infile} (default = FALSE).
}
  \item{ncores}{
a positive integer indicating the number of cores to be used in parallel computing (default = 1).
}
  \item{verbose}{
a logical switch for whether a progress bar should be shown for a GDS \code{infile} (default = FALSE).
}
}
\value{
NULL if \code{infile} is a BGEN file (.bgen) or a GDS file (.gds), otherwise computational time in seconds, excluding I/O time.
}
\references{
Chen, H., Wang, C., Conomos, M.P., Stilp, A.M., Li, Z., Sofer, T., Szpiro, A.A., Chen, W., Brehm, J.M., Celedón, J.C., Redline, S., Papanicolaou, G.J., Thornton, T.A., Laurie, C.C., Rice, K. and Lin, X. (2016) Control forpopulation structure and relatedness for binary traits in genetic association studies via logistic mixed models. The American Journal of Human Genetics 98, 653-666.
}
\author{
Han Chen, Duy T. Pham
}
\seealso{
\code{\link{glmmkin}}, \code{\link{glmm.wald}}
}
\examples{
\donttest{
data(example)
attach(example)
model0 <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id",
       family = binomial(link = "logit"))
plinkfiles <- strsplit(system.file("extdata", "geno.bed", package = "GMMAT"), 
       ".bed", fixed = TRUE)[[1]]
outfile <- "glmm.score.bed.testoutfile.txt"
glmm.score(model0, infile = plinkfiles, outfile = outfile)
infile <- system.file("extdata", "geno.gds", package = "GMMAT")
outfile <- "glmm.score.gds.testoutfile.txt"
glmm.score(model0, infile = infile, outfile = outfile)
infile <- system.file("extdata", "geno.txt", package = "GMMAT")
outfile <- "glmm.score.text.testoutfile.txt"
glmm.score(model0, infile = infile, outfile = outfile, infile.nrow.skip = 5, 
	infile.ncol.skip = 3, infile.ncol.print = 1:3, 
	infile.header.print = c("SNP", "Allele1", "Allele2"))
infile <- system.file("extdata", "geno.bgen", package = "GMMAT")
samplefile <- system.file("extdata", "geno.sample", package = "GMMAT")
outfile <- "glmm.score.bgen.testoutfile.txt"
glmm.score(model0, infile = infile, BGEN.samplefile = samplefile, outfile = outfile)
unlink(paste0("glmm.score.", c("bed", "gds", "text","bgen"), ".testoutfile.txt"))
}
}
\keyword{ generalized linear mixed model }
\keyword{ score test }
