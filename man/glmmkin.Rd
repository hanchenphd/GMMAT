\name{glmmkin}
\alias{glmmkin}
\title{
Fit generalized linear mixed model with known relationship matrices
}
\description{
Fit a generalized linear mixed model with a random intercept, or a random intercept and an optional random slope of time effect for longitudinal data. The covariance matrix of the random intercept is proportional to a known relationship matrix (e.g. kinship matrix in genetic association studies). Alternatively, it can be a variance components model with multiple random effects, and each component has a known relationship matrix.
}
\usage{
glmmkin(fixed, data = parent.frame(), kins = NULL, id, random.slope = NULL, 
	groups = NULL, family = binomial(link = "logit"), method = "REML", 
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
a known positive semi-definite relationship matrix (e.g. kinship matrix in genetic association studies) or a list of known positive semi-definite relationship matrices. The rownames and colnames of these matrices must at least include all samples as specified in the \code{id} column of the data frame \code{data}. If not provided, glmmkin will switch to the generalized linear model with no random effects (default = NULL).
}
  \item{id}{
a column in the data frame \code{data}, indicating the id of samples. When there are duplicates in \code{id}, the data is assumed to be longitudinal with repeated measures.
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

For longitudinal data (with duplicated \code{id}), two types of models can be applied: random intercept only models, and random intercept and random slope models. The random intercept only model is appropriate for analyzing repeated measures with no time trends, and observations for the same individual are assumed to be exchangeable. The random intercept and random slope model is appropriate for analyzing longitudinal data with individual-specific time trends (therefore, a random slope for time effect). Typically, the time effect should be included in the model as a fixed effect covariate as well. Covariances of the random intercept and the random slope are estimated.

For multiple phenotype analysis, \code{\link{formula}} recognized by \code{\link{lm}}, such as \code{cbind(y1, y2, y3) ~ x1 + x2}, can be used in \code{fixed} as fixed effects. For each matrix in \code{kins}, variance components corresponding to each phenotype, as well as their covariance components, will be estimated. Currently, \code{family} must be "gaussian" and \code{method.optim} must be "AI".
}
\value{
  \item{theta}{a vector or a list of variance component parameter estimates. See below.

For cross-sectional data, if \code{kins} is not provided (unrelated individuals), \code{theta} is the dispersion parameter estimate from the generalized linear model; if \code{kins} is a matrix and \code{groups} is not provided, \code{theta} is a length 2 vector, with \code{theta[1]} being the dispersion parameter estimate and \code{theta[2]} being the variance component parameter estimate for \code{kins}; if \code{kins} is a list and \code{groups} is not provided, \code{theta} is a length \code{1 + length(kins)} vector, with \code{theta[1]} being the dispersion parameter estimate and \code{theta[2:(1 + length(kins))]} being the variance component parameter estimates, corresponding to the order of matrices in the list \code{kins}; if \code{kins} is a matrix and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{1 + n.groups} vector, with \code{theta[1:n.groups]} being the residual variance estimates for each group and \code{theta[1 + n.groups]} being the variance component parameter estimate for \code{kins}; if \code{kins} is a list and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{length(kins) + n.groups} vector, with \code{theta[1:n.groups]} being the residual variance estimates for each group and \code{theta[(1 + n.groups):(length(kins) + n.groups)]} being the variance component parameter estimates, corresponding to the order of matrices in the list \code{kins}. 

For longitudinal data (with duplicated \code{id}) in a random intercept only model, if \code{kins} is not provided (unrelated individuals) and \code{groups} is not provided, \code{theta} is a length 2 vector, with \code{theta[1]} being the dispersion parameter estimate and \code{theta[2]} being the variance component parameter estimate for the random individual effects; if \code{kins} is a matrix and \code{groups} is not provided, \code{theta} is a length 3 vector, with \code{theta[1]} being the dispersion parameter estimate, \code{theta[2]} being the variance component parameter estimate for the random individual effects attributable to relatedness from \code{kins}, and \code{theta[3]} being the variance component parameter estimate for the random individual effects not attributable to relatedness from \code{kins}; if \code{kins} is a list and \code{groups} is not provided, \code{theta} is a length \code{2 + length(kins)} vector, with \code{theta[1]} being the dispersion parameter estimate, \code{theta[2:(1 + length(kins))]} being the variance component parameter estimates for the random individual effects attributable to relatedness from \code{kins}, corresponding to the order of matrices in the list \code{kins}, and \code{theta[2 + length(kins)]} being the variance component parameter estimate for the random individual effects not attributable to relatedness from \code{kins}; if \code{kins} is not provided (unrelated individuals) and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{1 + n.groups} vector, with \code{theta[1:n.groups]} being the residual variance estimates for each group and \code{theta[1 + n.groups]} being the variance component parameter estimate for the random individual effects; if \code{kins} is a matrix and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{2 + n.groups} vector, with \code{theta[1:n.groups]} being the residual variance estimates for each group, \code{theta[1 + n.groups]} being the variance component parameter estimate for the random individual effects attributable to relatedness from \code{kins}, and \code{theta[2 + n.groups]} being the variance component parameter estimate for the random individual effects not attributable to relatedness from \code{kins}; if \code{kins} is a list and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{1 + length(kins) + n.groups} vector, with \code{theta[1:n.groups]} being the residual variance estimates for each group, \code{theta[(1 + n.groups):(length(kins) + n.groups)]} being the variance component parameter estimates for the random individual effects attributable to relatedness from \code{kins}, corresponding to the order of matrices in the list \code{kins}, and \code{theta[1 + length(kins) + n.groups]} being the variance component parameter estimate for the random individual effects not attributable to relatedness from \code{kins}.

For longitudinal data (with duplicated \code{id}) in a random intercept and random slope (for time effect) model, if \code{kins} is not provided (unrelated individuals) and \code{groups} is not provided, \code{theta} is a length 4 vector, with \code{theta[1]} being the dispersion parameter estimate, \code{theta[2]} being the variance component parameter estimate for the random individual effects of the intercept, \code{theta[3]} being the covariance estimate for the random individual effects of the intercept and the random individual effects of the time slope, and \code{theta[4]} being the variance component parameter estimate for the random individual effects of the time slope; if \code{kins} is a matrix and \code{groups} is not provided, \code{theta} is a length 7 vector, with \code{theta[1]} being the dispersion parameter estimate, \code{theta[2]} being the variance component parameter estimate for the random individual effects of the intercept attributable to relatedness from \code{kins}, \code{theta[3]} being the variance component parameter estimate for the random individual effects of the intercept not attributable to relatedness from \code{kins}, \code{theta[4]} being the covariance estimate for the random individual effects of the intercept and the random individual effects of the time slope attributable to relatedness from \code{kins}, \code{theta[5]} being the covariance estimate for the random individual effects of the intercept and the random individual effects of the time slope not attributable to relatedness from \code{kins}, \code{theta[6]} being the variance component parameter estimate for the random individual effects of the time slope attributable to relatedness from \code{kins}, and \code{theta[7]} being the variance component parameter estimate for the random individual effects of the time slope not attributable to relatedness from \code{kins}; if \code{kins} is a list and \code{groups} is not provided, \code{theta} is a length \code{4 + 3 * length(kins)} vector, with \code{theta[1]} being the dispersion parameter estimate, \code{theta[2:(1 + length(kins))]} being the variance component parameter estimates for the random individual effects of the intercept attributable to relatedness from \code{kins}, corresponding to the order of matrices in the list \code{kins}, \code{theta[2 + length(kins)]} being the variance component parameter estimate for the random individual effects of the intercept not attributable to relatedness from \code{kins}, \code{theta[(3 + length(kins)):(2 + 2 * length(kins))]} being the covariance estimates for the random individual effects of the intercept and the random individual effects of the time slope attributable to relatedness from \code{kins}, corresponding to the order of matrices in the list \code{kins}, \code{theta[3 + 2 * length(kins)]} being the covariance estimate for the random individual effects of the intercept and the random individual effects of the time slope not attributable to relatedness from \code{kins}, \code{theta[(4 + 2 * length(kins)):(3 + 3 * length(kins))]} being the variance component parameter estimates for the random individual effects of the time slope attributable to relatedness from \code{kins}, corresponding to the order of matrices in the list \code{kins}, \code{theta[4 + 3 * length(kins)]} being the variance component parameter estimate for the random individual effects of the time slope not attributable to relatedness from \code{kins}; if \code{kins} is not provided (unrelated individuals) and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{3 + n.groups} vector, with \code{theta[1:n.groups]} being the residual variance estimates for each group, \code{theta[1 + n.groups]} being the variance component parameter estimate for the random individual effects of the intercept, \code{theta[2 + n.groups]} being the covariance estimate for the random individual effect of the intercept and the random individual effects of the time slope, and \code{theta[3 + n.groups]} being the variance component parameter estimate for the random individual effects of the time slope; if \code{kins} is a matrix and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{6 + n.groups} vector, with \code{theta[1:n.groups]} being the residual variance estimates for each group, \code{theta[1 + n.groups]} being the variance component parameter estimate for the random individual effects of the intercept attributable to relatedness from \code{kins}, \code{theta[2 + n.groups]} being the variance component parameter estimate for the random individual effects of the intercept not attributable to relatedness from \code{kins}, \code{theta[3 + n.groups]} being the covariance estimate for the random individual effects of the intercept and the random individual effects of the time slope attributable to relatedness from \code{kins}, \code{theta[4 + n.groups]} being the covariance estimate for the random individual effects of the intercept and the random individual effects of the time slope not attributable to relatedness from \code{kins}, \code{theta[5 + n.groups]} being the variance component parameter estimate for the random individual effects of the time slope attributable to relatedness from \code{kins}, and \code{theta[6 + n.groups]} being the variance component parameter estimate for the random individual effects of the time slope not attributable to relatedness from \code{kins}; if \code{kins} is a list and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{3 + 3 * length(kins) + n.groups} vector, with \code{theta[1:n.groups]} being the residual variance estimates for each group, \code{theta[(1 + n.groups):(length(kins) + n.groups)]} being the variance component parameter estimates for the random individual effects of the intercept attributable to relatedness from \code{kins}, corresponding to the order of matrices in the list \code{kins}, \code{theta[1 + length(kins) + n.groups]} being the variance component parameter estimate for the random individual effects of the intercept not attributable to relatedness from \code{kins}, \code{theta[(2 + length(kins) + n.groups):(1 + 2 * length(kins) + n.groups)]} being the covariance estimates for the random individual effects of the intercept and the random individual effects of the time slope attributable to relatedness from \code{kins}, corresponding to the order of matrices in the list \code{kins}, \code{theta[2 + 2 * length(kins) + n.groups]} being the covariance estimate for the random individual effects of the intercept and the random individual effects of the time slope not attributable to relatedness from \code{kins}, \code{theta[(3 + 2 * length(kins) + n.groups):(2 + 3 * length(kins) + n.groups)]} being the variance component parameter estimates for the random individual effects of the time slope attributable to relatedness from \code{kins}, corresponding to the order of matrices in the list \code{kins}, and \code{theta[3 + 3 * length(kins) + n.groups]} being the variance component parameter estimate for the random individual effects of the time slope not attributable to relatedness from \code{kins}.

For multiple phenotype analysis, \code{theta} is a list of variance-covariance matrices. If \code{kins} is not provided (unrelated individuals), \code{theta} is an \code{n.pheno} by \code{n.pheno} variance-covariance matrix for the residuals of the multiple phenotypes from the linear model; if \code{kins} is a matrix and \code{groups} is not provided, \code{theta} is a length 2 list, with \code{theta[[1]]} being the variance-covariance matrix for the residuals and \code{theta[[2]]} being the variance-covariance matrix for \code{kins}; if \code{kins} is a list and \code{groups} is not provided, \code{theta} is a length \code{1 + length(kins)} list, with \code{theta[[1]]} being the variance-covariance matrix for the residuals and \code{theta[[2]]} to \code{theta[[1 + length(kins)]]} being the variance-covariance matrices, corresponding to the order of matrices in the list \code{kins}; if \code{kins} is a matrix and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{1 + n.groups} list, with \code{theta[[1]]} to \code{theta[[n.groups]]} being the variance-covariance matrices for the residuals in each group and \code{theta[[1 + n.groups]]} being the variance-covariance matrix for \code{kins}; if \code{kins} is a list and \code{groups} is provided (a heteroscedastic linear mixed model with \code{n.groups} residual variance groups), \code{theta} is a length \code{length(kins) + n.groups} list, with \code{theta[[1]]} to \code{theta[[n.groups]]} being the variance-covariance matrices for the residuals in each group and \code{theta[[1 + n.groups]]} to \code{theta[[length(kins) + n.groups]]} being the variance-covariance matrices, corresponding to the order of matrices in the list \code{kins}.}
  \item{n.pheno}{an integer indicating the number of phenotypes in multiple phenotype analysis (for single phenotype analysis, \code{n.pheno = 1}).}
  \item{n.groups}{an integer indicating the number of distinct residual variance groups in heteroscedastic linear mixed models (for other models, \code{n.groups = 1}).}
  \item{coefficients}{a vector or a matrix for the fixed effects parameter estimates (including the intercept).}
  \item{linear.predictors}{a vector or a matrix for the linear predictors.}
  \item{fitted.values}{a vector or a matrix for the fitted mean values on the original scale.}
  \item{Y}{a vector or a matrix for the final working vector.}
  \item{X}{model matrix for the fixed effects.}
  \item{P}{the projection matrix with dimensions equal to the sample size multiplied by \code{n.pheno}. Used in \code{\link{glmm.score}} and \code{\link{SMMAT}} for dense matrices.}
  \item{residuals}{a vector or a matrix for the residuals on the original scale. NOT rescaled by the dispersion parameter.}
  \item{scaled.residuals}{a vector or a matrix for the scaled residuals, calculated as the original residuals divided by the dispersion parameter (in heteroscedastic linear mixed models, corresponding residual variance estimates by each group).}
  \item{cov}{covariance matrix for the fixed effects (including the intercept).}
  \item{Sigma_i}{the inverse of the estimated covariance matrix for samples, with dimensions equal to the sample size multiplied by \code{n.pheno}. Used in \code{\link{glmm.score}} and \code{\link{SMMAT}} for sparse matrices.}
  \item{Sigma_iX}{Sigma_i multiplied by X. Used in \code{\link{glmm.score}} and \code{\link{SMMAT}} for sparse matrices.}
  \item{converged}{a logical indicator for convergence.}
  \item{call}{the matched call.}
  \item{id_include}{a vector indicating the \code{id} of rows in \code{data} with nonmissing outcome and covariates, thus are included in the model fit.}
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
model0 <- glmmkin(disease ~ age + sex, data = pheno, kins = GRM, id = "id",
       family = binomial(link = "logit"))
model0$theta
model0$coefficients
model0$cov
\donttest{
model1 <- glmmkin(y.repeated ~ sex, data = pheno2, kins = GRM, id = "id", 
       family = gaussian(link = "identity"))
model1$theta
model1$coefficients
model1$cov
model2 <- glmmkin(y.trend ~ sex + time, data = pheno2, kins = GRM, id = "id", 
       random.slope = "time", family = gaussian(link = "identity"))
model2$theta
model2$coefficients
model2$cov
}
}
\keyword{ generalized linear mixed model }

