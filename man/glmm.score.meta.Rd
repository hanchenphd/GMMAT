\name{glmm.score.meta}
\alias{glmm.score.meta}
\title{
Performing meta-analysis for GLMM based score test results
}
\description{
Use output files from GLMM based score tests to perform meta-analysis.
}
\usage{
glmm.score.meta(files, outfile, SNP = rep("SNP", length(files)), 
	A1 = rep("A1", length(files)), A2 = rep("A2", length(files)))
}
\arguments{
  \item{files}{
a vector of input file names. The input files should be the output files of \code{glmm.score( )}, or customized tab or space delimited files that include at least 8 columns: SNP, effect allele, noneffect allele, N, AF, SCORE, VAR and PVAL. The column names of SNP, effect allele and noneffect allele can be customized and provided in \code{SNP}, \code{A1} and \code{A2}.
}
  \item{outfile}{
the output file name.
}
  \item{SNP}{
a character vector of SNP column names in each input file. The length and order must match the length and order of \code{files} (default = rep("SNP", length(files))).
}
  \item{A1}{
a character vector of allele 1 column names in each input file. The length and order must match the length and order of \code{files} (default = rep("A1", length(files))). Note that \code{glmm.score.meta( )} does not define A1 as the effect allele or noneffect allele: it is the user's choice. However, the choice should be consistent across different studies, if A1 column is the effect allele in one study but the noneffect allele in another, meta-analysis results will be incorrect.
}
  \item{A2}{
a character vector of allele 2 column names in each input file. The length and order must match the length and order of \code{files} (default = rep("A2", length(files))). Note that \code{glmm.score.meta( )} does not define A2 as the effect allele or noneffect allele: it is the user's choice. However, the choice should be consistent across different studies, if A2 column is the effect allele in one study but the noneffect allele in another, meta-analysis results will be incorrect.
}
}
\value{
a data frame containing the following:
  \item{SNP}{SNP name.}
  \item{A1}{allele 1.}
  \item{A2}{allele 2.}
  \item{N}{total sample size.}
  \item{AF}{effect allele frequency (user-defined: can be either allele 1 or allele 2).}
  \item{SCORE}{the summary score of the effect allele.}
  \item{VAR}{the variance of the summary score.}
  \item{PVAL}{meta-analysis p-value.}
}
\author{
Han Chen
}
\seealso{
\code{\link{glmm.score}}
}
\examples{
\donttest{
infile1 <- system.file("extdata", "meta1.txt", package = "GMMAT")
infile2 <- system.file("extdata", "meta2.txt", package = "GMMAT")
infile3 <- system.file("extdata", "meta3.txt", package = "GMMAT")
outfile <- "glmm.score.meta.testoutfile.txt"
glmm.score.meta(files = c(infile1, infile2, infile3), outfile = outfile,
	SNP = rep("SNP", 3), A1 = rep("A1", 3), A2 = rep("A2", 3))
unlink(outfile)
}
}
\keyword{ meta-analysis }
\keyword{ score test }
