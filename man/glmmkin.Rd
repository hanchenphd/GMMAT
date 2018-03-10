\name{glmmkin}
\alias{glmmkin}
\title{
Fit generalized linear mixed model with known relationship matrices
}
\description{
Fit a generalized linear mixed model with a random intercept. The covariance matrix of the random intercept is proportional to a known relationship matrix (e.g. kinship matrix in genetic association studies). Alternatively, it can be a variance components model with multiple random effects, and each component has a known relationship matrix.
}
\usage{
glmmkin(fixed, data = parent.frame(), kins, groups = NULL, 
	family = binomial(link = "logit"), method = "REML", 
	method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, 
	taumax = 1e5, tauregion = 10, verbose = FALSE, ...)
}
\arguments{
  \item{fixed}{
an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the fixed effects model to be fitted.
}
  \item{data}{
a data frame or list (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model.
}
  \item{kins}{
a known relationship matrix (e.g. kinship matrix in genetic association studies) or a list of known relationship matrices. When supplied with a matrix, this matrix should be a positive semi-definite matrix with dimensions equal to the sample size in \code{data}, and the order of subjects in this matrix should also match the order of subjects in \code{data}. When supplied with a list of matrices, each matrix should be a positive semi-definite matrix with dimensions equal to the sample size in \code{data}, and the order of subjects in each matrix should also match the order of subjects in \code{data}. Currently there is no ID checking feature implemented, and it is the user's responsibility to match the orders.
}
  \item{groups}{
an optional categorical variable indicating the groups used in a heteroscedastic linear mixed model (allowing residual variances in different groups to be different). This variable must be included in the names of \code{data}, and \code{family} must be "gaussian" and \code{method.optim} must be "AI" (default = NULL).
}
  \item{family}{
a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.)
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
a positive number specifying tolerance, the difference threshold for parameter estimates below which iterations should be stopped (default = 1e-5).
}
  \item{taumin}{
the lower bound of search space for the variance component parameter \eqn{\tau} (default = 1e-5), used when \code{method.optim = "Brent"}. See Details.
}
  \item{taumax}{
the upper bound of search space for the variance component parameter \eqn{\tau} (default = 1e5), used when \code{method.optim = "Brent"}. See Details.
}
  \item{tauregion}{
the number of search intervals for the REML or ML estimate of the variance component parameter \eqn{\tau} (default = 10), used when \code{method.optim = "Brent"}. See Details.
}
  \item{verbose}{
a logical switch for printing detailed information (parameter estimates in each iteration) for testing and debugging purpose (default = FALSE).
}
  \item{\dots}{
additional arguments that could be passed to \code{\link{glm}}.
}
}
\details{
Generalized linear mixed models (GLMM) are fitted using the penalized quasi-likelihood (PQL) method proposed by Breslow and Clayton (1993). Generally, fitting a GLMM is computationally expensive, and by default we use the Average Information REML algorithm (Gilmour, Thompson and Cullis, 1995; Yang et al., 2011) to fit the model. If only one relationship matrix is specified (\code{kins} is a matrix), iterations may be accelerated using the algorithm proposed by Zhou and Stephens (2012) for linear mixed models. An eigendecomposition is performed in each outer iteration and the estimate of the variance component parameter \eqn{\tau} is obtained by maximizing the profiled log restricted likelihood (or likelihood) in a search space from \code{taumin} to \code{taumax}, equally divided into \code{tauregion} intervals on the log scale, using Brent's method (1973). If \code{kins} is a list of matrices and \code{method = "Nelder-Mead"}, iterations are performed as a multi-dimensional maximization problem solved by Nelder and Mead's method (1965). It can be very slow, and we do not recommend using this method unless the likelihood function is badly behaved. Both Brent's method and Nelder and Mead's method are derivative-free. When the Average Information REML algorithm fails to converge, a warning message is given and the algorithm is default to derivative-free approaches: Brent's method if only one relationship matrix is specified, Nelder and Mead's method if more than one relationship matrix is specified.
}
\value{
  \item{theta}{a vector of length \code{1 + n.groups} for the dispersion parameter estimate (or in the case of heteroscedastic linear mixed models with \code{n.groups > 1}, the residual variance estimates: \code{theta[1:n.groups]}) and the variance component parameter estimate (\code{theta[1 + n.groups]}) if \code{kins} is a matrix; a vector of length equal to \code{length(kins) + n.groups} if \code{kins} is a list, \code{theta[1:n.groups]} is the dispersion parameter estimate (or in the case of heteroscedastic linear mixed models with \code{n.groups > 1}, the residual variance estimates) and the remaining are the variance components parameter estimates, with the same order as in \code{kins}.}
  \item{n.groups}{an integer indicating the number of distinct residual variance groups in heteroscedastic linear mixed models (for other models, \code{n.groups = 1}).}
  \item{coefficients}{a vector for the fixed effects parameter estimates (including the intercept).}
  \item{linear.predictors}{a vector of length equal to the sample size for the linear predictors.}
  \item{fitted.values}{a vector of length equal to the sample size for the fitted mean values on the original scale.}
  \item{Y}{a vector of length equal to the sample size for the final working vector.}
  \item{P}{the projection matrix with dimensions equal to the sample size. Used in score tests.}
  \item{residuals}{a vector of length equal to the sample size for the residuals on the original scale. NOT rescaled by the dispersion parameter.}
  \item{scaled.residuals}{a vector of length equal to the sample size for the scaled residuals, calculated as the original residuals divided by the dispersion parameter (in heteroscedastic linear mixed models, corresponding residual variance estimates by each group).}
  \item{cov}{covariance matrix for the fixed effects (including the intercept).}
  \item{converged}{a logical indicator for convergence.}
  \item{call}{the matched call.}
  \item{id_include}{a vector indicating which rows in \code{data} have nonmissing outcome and covariates, thus are included in the model fit.}
}
\references{
Brent, R.P. (1973) "Chapter 4: An Algorithm with Guaranteed Convergence for Finding a Zero of a Function", Algorithms for Minimization without Derivatives, Englewood Cliffs, NJ: Prentice-Hall, ISBN 0-13-022335-2.

Breslow, N.E. and Clayton, D.G. (1993) Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association 88, 9-25.

Chen, H., Wang, C., Conomos, M.P., Stilp, A.M., Li, Z., Sofer, T., Szpiro, A.A., Chen, W., Brehm, J.M., Celedón, J.C., Redline, S., Papanicolaou, G.J., Thornton, T.A., Laurie, C.C., Rice, K. and Lin, X. (2016) Control for population structure and relatedness for binary traits in genetic association studies via logistic mixed models. The American Journal of Human Genetics 98, 653-666.

Gilmour, A.R., Thompson, R. and Cullis, B.R. (1995) Average Information REML: An Efficient Algorithm for Variance Parameter Estimation in Linear Mixed Models. Biometrics 51, 1440-1450.

Nelder, J.A. and Mead, R. (1965) A simplex algorithm for function minimization. Computer Journal 7, 308-313.

Yang, J., Lee, S.H., Goddard, M.E. and Visscher, P.M. (2011) GCTA: A Tool for Genome-wide Complex Trait Analysis. The American Journal of Human Genetics 88, 76-82.

Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis for association studies. Nature Genetics 44, 821-824.
}
\author{
Han Chen, Matthew P. Conomos
}
\examples{
data(example)
attach(example)
model0 <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, 
       family = binomial(link = "logit"))
model0$theta
model0$coefficients
model0$cov
}
\keyword{ generalized linear mixed model }

