glmmkin <- function(fixed, data = parent.frame(), kins = NULL, id, random.slope = NULL, groups = NULL, family = binomial(link = "logit"), method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE, ...) {
  call <- match.call()
  if(!is.null(kins) && !class(kins)[1] %in% c("matrix", "list")) {
    if(is.null(attr(class(kins), "package"))) stop("Error: \"kins\" must be a matrix or a list.")
    else if(attr(class(kins), "package") != "Matrix") stop("Error: if \"kins\" is a sparse matrix, it must be created using the Matrix package.")
  }
  if(!method %in% c("REML", "ML"))
    stop("Error: \"method\" must be \"REML\" or \"ML\".")
  method.optim <- try(match.arg(method.optim, c("AI", "Brent", "Nelder-Mead")))
  if(class(method.optim) == "try-error")
    stop("Error: \"method.optim\" must be \"AI\", \"Brent\" or \"Nelder-Mead\".")
  if(method.optim == "AI" && method == "ML")
    stop("Error: method \"ML\" not available for method.optim \"AI\", use method \"REML\" instead.")
  if(method.optim == "Brent" && class(kins)[1] == "list")
    stop("Error: method.optim \"Brent\" can only be applied in one-dimensional optimization, use a matrix for \"kins\".")
  if(method.optim != "AI" && ((!is.null(attr(class(kins), "package")) && attr(class(kins), "package") == "Matrix") || (class(kins)[1] == "list" && any(sapply(kins, function(xx) !is.null(attr(class(xx), "package")) && attr(class(xx), "package") == "Matrix")))))
    stop("Error: sparse matrices can only be handled by method.optim \"AI\".")
  if(class(family) != "family")
    stop("Error: \"family\" must be an object of class \"family\".")
  if(!family$family %in% c("binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"))
    stop("Error: \"family\" must be one of the following: binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, quasipoisson.")
  if(!is.null(groups)) {
    if(family$family != "gaussian") stop("Error: heteroscedastic linear mixed models are only applicable when \"family\" is gaussian.")
    if(method.optim != "AI") stop("Error: heteroscedastic linear mixed models are currently only implemented for method.optim \"AI\".")
    if(!groups %in% names(data)) stop("Error: \"groups\" must be one of the variables in the names of \"data\".")
  }
  if(!id %in% names(data)) stop("Error: \"id\" must be one of the variables in the names of \"data\".")
  if("data.frame" %in% class(data) && length(class(data)) > 1) data <- as.data.frame(data)
  if(!is.null(random.slope)) {
    if(method.optim != "AI") stop("Error: random slope for longitudinal data is currently only implemented for method.optim \"AI\".")
    if(!random.slope %in% names(data)) stop("Error: \"random.slope\" must be one of the variables in the names of \"data\".")
  }
  if(!is.null(attr(class(kins), "package")) && attr(class(kins), "package") == "Matrix")  kins <- list(kins1 = kins)
  if(method.optim != "Brent" && class(kins)[1] == "matrix") kins <- list(kins1 = kins)
  mdl <- model.frame(formula = fixed, data = data, na.action = na.omit)
  idx <- match(rownames(mdl), rownames(model.frame(formula = fixed, data = data, na.action = na.pass)))
  multi.pheno <- !is.null(ncol(mdl[,1])) && ncol(mdl[,1]) > 1
  rm(mdl)
  if(multi.pheno) {
    cat("Fitting a multiple phenotype analysis model...\n")
    if(method != "REML") stop("Error: \"method\" must be \"REML\".")
    if(method.optim != "AI") stop("Error: \"method.optim\" must be \"AI\".")
    if(family$family != "gaussian") stop("Error: \"family\" must be gaussian.")
    fit0 <- do.call("lm", list(formula = fixed, data = data, ...))
  } else {
    fit0 <- do.call("glm", list(formula = fixed, data = data, family = family, ...))
  }
  if(any(duplicated(data[idx, id]))) {
    if(multi.pheno) stop("Error: Duplicated id detected...\nNot currently supported for multiple phenotypes...\n")
    cat("Duplicated id detected...\nAssuming longitudinal data with repeated measures...\n")
    if(method.optim == "Brent") {
      if(is.null(kins)) {
        kins <- diag(length(unique(data[idx, id])))
        rownames(kins) <- colnames(kins) <- unique(data[idx, id])
      } else stop("Error: method.optim \"Brent\" can only be applied to unrelated individuals in longitudinal data analysis.")
    } else {
      if(method.optim != "AI") kins[[length(kins) + 1]] <- diag(length(unique(data[idx,id]))) 
      else if(length(kins) > 0) kins[[length(kins) + 1]] <- as(diag(length(unique(data[idx, id]))), "sparseMatrix")
      else kins <- list(kins1 = as(diag(length(unique(data[idx, id]))), "sparseMatrix"))
      rownames(kins[[length(kins)]]) <- colnames(kins[[length(kins)]]) <- unique(data[idx, id])
    }
  } else if(!is.null(random.slope)) stop("Error: no duplicated \"id\" found, \"random.slope\" must be used for longitudinal data with duplicated \"id\".")
  if(class(kins)[1] == "matrix") {
    match.idx1 <- match(data[idx, id], rownames(kins))
    match.idx2 <- match(data[idx, id], colnames(kins))
    if(any(is.na(c(match.idx1, match.idx2)))) stop("Error: kins matrix does not include all individuals in the data.")
    kins <- kins[match.idx1, match.idx2]
  } else if(class(kins)[1] == "list") {
    for(i in 1:length(kins)) {
      match.idx1 <- match(data[idx, id], rownames(kins[[i]]))
      match.idx2 <- match(data[idx, id], colnames(kins[[i]]))
      if(any(is.na(c(match.idx1, match.idx2)))) stop("Error: kins matrix ", i, " does not include all individuals in the data.")
      kins[[i]] <- kins[[i]][match.idx1, match.idx2]
    }
  } else {
    if(!is.null(groups)) stop("Error: heteroscedastic linear models for unrelated observations have not been implemented.")
    y <- fit0$model[[1]]
    n <- if(multi.pheno) nrow(y) else length(y)
    n.pheno <- if(multi.pheno) ncol(y) else NULL
    offset <- fit0$offset
    if(is.null(offset)) offset <- rep(0, n)
    X <- model.matrix(fit0)
    alpha <- fit0$coef
    if(multi.pheno) {
      mu <- fit0$fitted.values
      Y <- y - offset
      res <- fit0$residuals
      tau <- crossprod(res)/fit0$df.residual
      Sigma_i <- solve(tau) %x% Diagonal(n = n)
      rownames(Sigma_i) <- colnames(Sigma_i) <- rep(rownames(X), n.pheno)
      fit <- list(theta=list(tau), n.pheno=n.pheno, n.groups=1, coefficients=alpha, linear.predictors=mu, fitted.values=mu, Y=Y, X=X, P=NULL, residuals=res, scaled.residuals=tcrossprod(res, solve(tau)), cov=vcov(fit0), Sigma_i=Sigma_i,Sigma_iX=crossprod(Sigma_i, Diagonal(n=n.pheno) %x% X),converged=TRUE, call = call, id_include = data[idx, id])
      class(fit) <- "glmmkin.multi"
    } else {
      family <- fit0$family
      eta <- fit0$linear.predictors
      mu <- fit0$fitted.values
      mu.eta <- family$mu.eta(eta)
      Y <- eta - offset + (y - mu)/mu.eta
      sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0))*family$variance(mu))
      res <- y - mu
      tau <- summary(fit0)$dispersion
      Sigma_i <- Diagonal(x = sqrtW^2/tau)
      rownames(Sigma_i) <- colnames(Sigma_i) <- rownames(X)
      fit <- list(theta=tau, n.pheno=1, n.groups=1, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, X=X, P=NULL, residuals=res, scaled.residuals=res*as.vector(weights(fit0))/tau, cov=vcov(fit0), Sigma_i=Sigma_i,Sigma_iX=X*sqrtW^2/tau,converged=TRUE, call = call, id_include = data[idx, id])
      class(fit) <- "glmmkin"
    }
    return(fit)
  }
  group.id <- if(is.null(groups)) rep(1, length(idx)) else data[idx, groups]
  time.var <- if(is.null(random.slope)) NULL else data[idx, random.slope]
  fit <- glmmkin.fit(fit0, kins, time.var, group.id, method = method, method.optim = method.optim, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
  fit$call <- call
  fit$id_include <- data[idx, id]
  class(fit) <- if(multi.pheno) "glmmkin.multi" else "glmmkin"
  return(fit)
}

glmmkin.fit <- function(fit0, kins, time.var, group.id, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE) {
  if(method.optim == "Brent") {
    fit <- glmmkin.brent(fit0, kins, method = method, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
    if(fit$theta[2]/fit$theta[1] < 1.01 * tol) {
      warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
      fit <- glmmkin.brent(fit0, kins, method = method, tau = 0, fixtau = 1, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
    }
  } else {
    names(kins) <- paste("kins", 1:length(kins), sep="")
    if(method.optim == "AI") {
      group.unique <- unique(group.id)
      group.idx <- list()
      for(i in 1:length(group.unique)) group.idx[[i]] <- which(group.id == group.unique[i])
      covariance.idx <- fixrho.old <- fixrho.new <- NULL
      n <- nrow(fit0$model[[1]])
      n.pheno <- ncol(fit0$model[[1]])
      q <- length(kins)
      ng <- length(group.idx)
      if(!is.null(n.pheno)) {
        n.params <- n.pheno*(n.pheno+1)/2*(q+ng)
        covariance.idx <- matrix(0, n.pheno*(n.pheno-1)/2*(q+ng), 3)
        kins.old <- kins
        kins.idx <- 1
        for(ii in 1:ng) {
          diagSigma <- rep(0, n)
          diagSigma[group.idx[[ii]]] <- 1
          for(jj in 1:n.pheno) {
            for(kk in jj:n.pheno) {
              kins[[kins.idx]] <- sparseMatrix(i=jj,j=kk,x=1,symmetric=TRUE,dims=c(n.pheno,n.pheno)) %x% Diagonal(x = diagSigma)
              kins.idx <- kins.idx + 1
              if(kk > jj) covariance.idx[(ii-1)*n.pheno*(n.pheno-1)/2+(2*n.pheno-jj)*(jj-1)/2+kk-jj, ] <- c((2*n.pheno+2-jj)*(jj-1)/2+kk-jj+1, (2*n.pheno+2-jj)*(jj-1)/2+1, (2*n.pheno+2-kk)*(kk-1)/2+1) + (ii-1)*n.pheno*(n.pheno+1)/2
            }
          }
        }
        for(ii in 1:q) {
          for(jj in 1:n.pheno) {
            for(kk in jj:n.pheno) {
              kins[[kins.idx]] <- sparseMatrix(i=jj,j=kk,x=1,symmetric=TRUE,dims=c(n.pheno,n.pheno)) %x% kins.old[[ii]]
              kins.idx <- kins.idx + 1
              if(kk > jj) covariance.idx[(ng+ii-1)*n.pheno*(n.pheno-1)/2+(2*n.pheno-jj)*(jj-1)/2+kk-jj, ] <- c((2*n.pheno+2-jj)*(jj-1)/2+kk-jj+1, (2*n.pheno+2-jj)*(jj-1)/2+1, (2*n.pheno+2-kk)*(kk-1)/2+1) + (ng+ii-1)*n.pheno*(n.pheno+1)/2
            }
          }
        }
        rm(kins.old)
        gc()
        names(kins) <- paste("kins", 1:length(kins), sep="")
        fixrho.old <- rep(0, nrow(covariance.idx))
        fixtau.old <- rep(0, n.params)
        fit <- glmmkin.multi.ai(fit0 = fit0, kins = kins, covariance.idx = covariance.idx, group.idx = group.idx, maxiter = maxiter, tol = tol, verbose = verbose)
      } else {
        if(!is.null(time.var)) {
          covariance.idx <- matrix(0, q, 3)
          for(i in 1:q) kins[[q + i]] <- if(!is.null(attr(class(kins[[i]]), "package")) && attr(class(kins[[i]]), "package") == "Matrix") forceSymmetric(kins[[i]] * time.var + t(kins[[i]] * time.var)) else kins[[i]] * time.var + t(kins[[i]] * time.var)
          for(i in 1:q) {
            kins[[2*q + i]] <- if(!is.null(attr(class(kins[[i]]), "package")) && attr(class(kins[[i]]), "package") == "Matrix") forceSymmetric(t(kins[[i]] * time.var) * time.var) else t(kins[[i]] * time.var) * time.var
            covariance.idx[i, ] <- c(q + i, i, 2*q + i) + length(group.idx)
          }
          names(kins) <- paste("kins", 1:length(kins), sep="")
          fixrho.old <- rep(0, q)
        }
        fixtau.old <- rep(0, length(kins)+ng)
        fit <- glmmkin.ai(fit0 = fit0, kins = kins, covariance.idx = covariance.idx, group.idx = group.idx, maxiter = maxiter, tol = tol, verbose = verbose)
      }
      fixtau.new <- 1*(fit$theta < 1.01 * tol)
      if(!is.null(covariance.idx)) {
        fixtau.new[covariance.idx[, 1]] <- 0
        fixrho.new <- rep(0, nrow(covariance.idx))
        fixrho.idx <- apply(covariance.idx, 1, function(x) abs(fit$theta[x[1]]) > (1 - 1.01 * tol) * sqrt(fit$theta[x[2]] * fit$theta[x[3]]))
        fixrho.new[fixrho.idx] <- sign(fit$theta[covariance.idx[fixrho.idx, 1]])
      }
      while(any(fixtau.new != fixtau.old) || (!is.null(fixrho.new) && any(fixrho.new != fixrho.old))) {
        warning("Variance estimate on the boundary of the parameter space observed, refitting model...", call. = FALSE)
        fixtau.old <- fixtau.new
        if(!is.null(covariance.idx)) fixrho.old <- fixrho.new
        if(!is.null(n.pheno)) {
          fit <- glmmkin.multi.ai(fit0 = fit0, kins = kins, covariance.idx = covariance.idx, group.idx = group.idx, fixtau = fixtau.old, fixrho = fixrho.old, maxiter = maxiter, tol = tol, verbose = verbose)
        } else {
          fit <- glmmkin.ai(fit0 = fit0, kins = kins, covariance.idx = covariance.idx, group.idx = group.idx, fixtau = fixtau.old, fixrho = fixrho.old, maxiter = maxiter, tol = tol, verbose = verbose)
        }
        fixtau.new <- 1*(fit$theta < 1.01 * tol)
        if(!is.null(covariance.idx)) {
          fixtau.new[covariance.idx[, 1]] <- 0
          fixrho.new <- rep(0, nrow(covariance.idx))
          fixrho.idx <- apply(covariance.idx, 1, function(x) abs(fit$theta[x[1]]) > (1 - 1.01 * tol) * sqrt(fit$theta[x[2]] * fit$theta[x[3]]))
          fixrho.new[fixrho.idx] <- sign(fit$theta[covariance.idx[fixrho.idx, 1]])
        }
      }
      if(!is.null(n.pheno)) {
        fit$n.groups <- ng
        scaled.res <- fit$residuals
        theta <- list()
        for(ii in 1:ng) {
          theta[[ii]] <- sparseMatrix(i=rep(1:n.pheno,n.pheno:1),j=unlist(lapply(1:n.pheno,function(xx) xx:n.pheno)),x=fit$theta[(ii-1)*n.pheno*(n.pheno+1)/2+(1:(n.pheno*(n.pheno+1)/2))],symmetric=T,dims=c(n.pheno,n.pheno))
          scaled.res[group.idx[[ii]],] <- as.matrix(tcrossprod(scaled.res[group.idx[[ii]],,drop = FALSE], solve(theta[[ii]])))
          if(ng > 1) names(theta)[ii] <- group.unique[ii]
          else names(theta)[1] <- "residuals"
        }
        for(ii in 1:q) {
          theta[[ng+ii]] <- sparseMatrix(i=rep(1:n.pheno,n.pheno:1),j=unlist(lapply(1:n.pheno,function(xx) xx:n.pheno)),x=fit$theta[(ng+ii-1)*n.pheno*(n.pheno+1)/2+(1:(n.pheno*(n.pheno+1)/2))],symmetric=T,dims=c(n.pheno,n.pheno))
          names(theta)[ng+ii] <- paste0("kins", ii)
        }
        fit$theta <- theta
        fit$scaled.residuals <- scaled.res
      } else {
        if(ng > 1) names(fit$theta)[1:ng] <- group.unique
        else names(fit$theta)[1] <- "dispersion"
        if(!is.null(time.var)) {
          names(fit$theta)[ng+(1:q)] <- paste0("kins", 1:q, ".var.intercept")
          names(fit$theta)[ng+q+(1:q)] <- paste0("kins", 1:q, ".cov.intercept.slope")
          names(fit$theta)[ng+2*q+(1:q)] <- paste0("kins", 1:q, ".var.slope")
        } else {
          names(fit$theta)[ng+(1:q)] <- paste0("kins", 1:q)
        }
      }
      if(!fit$converged) {
        if(length(group.idx) != 1) stop("Error: Average Information REML not converged, cannot refit heteroscedastic linear mixed model using Brent or Nelder-Mead methods.")
        if(!is.null(time.var)) stop("Error: Average Information REML not converged, cannot refit random slope model for longitudinal data using Brent or Nelder-Mead methods.")
        if(!is.null(n.pheno)) stop("Error: Average Information REML not converged, cannot refit multiple phenotype model using Brent or Nelder-Mead methods.")
        if(length(kins) == 1) {
          warning("Average Information REML not converged, refitting model using Brent method...", call. = FALSE)
          fit <- glmmkin.brent(fit0, kins[[1]], method = method, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
          if(fit$theta[2]/fit$theta[1] < 1.01 * tol) {
            warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
            fit <- glmmkin.brent(fit0, kins[[1]], method = method, tau = 0, fixtau = 1, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
          }
        } else {
          warning("Average Information REML not converged, refitting model using Nelder-Mead method...", call. = FALSE)
          fixtau.old <- rep(0, length(kins))
          fit <- glmmkin.nm(fit0, kins, method = method, maxiter = maxiter, tol = tol, verbose = verbose)
          fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
          while(any(fixtau.new != fixtau.old)) {
            warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
            fixtau.old <- fixtau.new
            tau <- rep(1, length(kins))
            tau[which(fixtau.old == 1)] <- 0
            fit <- glmmkin.nm(fit0, kins, method = method, tau = tau, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
            fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
          }
        }
      }
    } else {
      fixtau.old <- rep(0, length(kins))
      fit <- glmmkin.nm(fit0, kins, method = method, maxiter = maxiter, tol = tol, verbose = verbose)
      fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
      while(any(fixtau.new != fixtau.old)) {
        warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
        fixtau.old <- fixtau.new
        tau <- rep(1, length(kins))
        tau[which(fixtau.old == 1)] <- 0
        fit <- glmmkin.nm(fit0, kins, method = method, tau = tau, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
        fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
      }
    }
  }
  return(fit)
}

glmmkin.ai <- function(fit0, kins, covariance.idx = NULL, group.idx, tau = rep(0, length(kins)+length(group.idx)), fixtau = rep(0, length(kins)+length(group.idx)), fixrho = NULL, maxiter = 500, tol = 1e-5, verbose = FALSE) {
  is.Matrix <- any(sapply(kins, function(xx) !is.null(attr(class(xx), "package")) && attr(class(xx), "package") == "Matrix"))
  y <- fit0$y
  n <- length(y)
  offset <- fit0$offset
  if(is.null(offset)) offset <- rep(0, n)
  family <- fit0$family
  eta <- fit0$linear.predictors
  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu)/mu.eta
  sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0))*family$variance(mu))
  X <- model.matrix(fit0)
  alpha <- fit0$coef
  if(verbose) {
    cat("Fixed-effect coefficients:\n")
    print(alpha)
  }
  if(family$family %in% c("poisson", "binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
  }
  q <- length(kins)
  ng <- length(group.idx)
  idxtau <- which(fixtau == 0)
  if(!is.null(covariance.idx)) idxtau2 <- intersect(covariance.idx[, 1], idxtau)
  q2 <- sum(fixtau == 0)
  if(q2 > 0) {
    tau[idxtau] <- rep(var(Y)/(q+ng), q2)
    if(!is.null(covariance.idx)) tau[idxtau2] <- 0
    diagSigma <- rep(0, n)
    for(i in 1:ng) diagSigma[group.idx[[i]]] <- tau[i]/sqrtW[group.idx[[i]]]^2
    Sigma <- if(is.Matrix) Diagonal(x = diagSigma) else diag(diagSigma)
    for(i in 1:q) {
      tau[i+ng] <- tau[i+ng]/mean(diag(kins[[i]]))
      Sigma <- Sigma + tau[i+ng]*kins[[i]]
    }
    #		if(is.Matrix) Sigma <- forceSymmetric(Sigma)
    Sigma_i <- chol2inv(chol(Sigma))
    rm(Sigma, diagSigma)
    gc()
    Sigma_iX <- crossprod(Sigma_i, X)
    #P <- Sigma_i - tcrossprod(tcrossprod(Sigma_iX, chol2inv(chol(crossprod(X, Sigma_iX)))), Sigma_iX)
    #rm(Sigma_i)
    #gc()
    #PY <- crossprod(P, Y)
    XSigma_iX <- crossprod(X, Sigma_iX)
    if(is.Matrix) XSigma_iX <- forceSymmetric(XSigma_iX)
    cov <- chol2inv(chol(XSigma_iX))
    Sigma_iXcov <- tcrossprod(Sigma_iX, cov)
    PY <- crossprod(Sigma_i, Y) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, Y)))
    tau0 <- tau
    for(i in 1:q2) {
      #if(idxtau[i] <= ng) tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (sum((PY/sqrtW)[group.idx[[idxtau[i]]]]^2) - sum((diag(P)/sqrtW^2)[group.idx[[idxtau[i]]]]))/n)
      if(idxtau[i] <= ng) tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (sum((PY/sqrtW)[group.idx[[idxtau[i]]]]^2) - sum(((diag(Sigma_i)-rowSums(Sigma_iX*Sigma_iXcov))/sqrtW^2)[group.idx[[idxtau[i]]]]))/n)
      else {
        #PAPY <- crossprod(P, crossprod(kins[[idxtau[i]-ng]], PY))
        #if(!is.null(covariance.idx)) tau[idxtau[i]] <- if(idxtau[i] %in% idxtau2) 0 else max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (crossprod(Y, PAPY) - sum(P*kins[[idxtau[i]-ng]]))/n)
        #else tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (crossprod(Y, PAPY) - sum(P*kins[[idxtau[i]-ng]]))/n)
        APY <- crossprod(kins[[idxtau[i]-ng]], PY)
        PAPY <- crossprod(Sigma_i, APY) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, APY)))
        if(!is.null(covariance.idx)) tau[idxtau[i]] <- if(idxtau[i] %in% idxtau2) 0 else max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (sum(Y* PAPY) - (sum(Sigma_i*kins[[idxtau[i]-ng]])-sum(Sigma_iX*crossprod(kins[[idxtau[i]-ng]],Sigma_iXcov))))/n)
        else tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (sum(Y* PAPY) - (sum(Sigma_i*kins[[idxtau[i]-ng]])-sum(Sigma_iX*crossprod(kins[[idxtau[i]-ng]],Sigma_iXcov))))/n)
      }
    }
    #rm(P)
    #gc()
  }
  for (i in seq_len(maxiter)) {
    if(verbose) cat("\nIteration ", i, ":\n")
    alpha0 <- alpha
    tau0 <- tau
    #fit <- .Call(C_fitglmm_ai, Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau, tol)
    if(is.Matrix) fit <- R_fitglmm_ai(Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau)
    else if(n < 2^15.5) fit <- .Call(C_fitglmm_ai, Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau)
    else fit <- R_fitglmm_ai_dense(Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau)
    if(q2 > 0) {
      Dtau <- as.numeric(fit$Dtau)
      tau[idxtau] <- tau0[idxtau] + Dtau
      if(is.null(covariance.idx)) {
        tau[tau < tol & tau0 < tol] <- 0
        while(any(tau < 0)) {
          Dtau <- Dtau / 2
          tau[idxtau] <- tau0[idxtau] + Dtau
          tau[tau < tol & tau0 < tol] <- 0
        }
        tau[tau < tol] <- 0
      } else {
        fixrho.idx0 <- apply(covariance.idx, 1, function(x) abs(tau0[x[1]]) > (1 - 1.01 * tol) * sqrt(tau0[x[2]] * tau0[x[3]]))
        tau[-covariance.idx[, 1]][tau[-covariance.idx[, 1]] < tol & tau0[-covariance.idx[, 1]] < tol] <- 0
        if(any(fixrho != 0)) {
          idxrho <- which(fixrho != 0)
          tau[covariance.idx[idxrho, 1]] <- suppressWarnings(fixrho[idxrho] * sqrt(tau[covariance.idx[idxrho, 2]] * tau[covariance.idx[idxrho, 3]]))
        }
        fixrho.idx <- suppressWarnings(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) > (1 - 1.01 * tol) * sqrt(tau[x[2]] * tau[x[3]])))
        tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]] <- sign(tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]]) * sqrt(tau[covariance.idx[fixrho.idx & fixrho.idx0, 2]] * tau[covariance.idx[fixrho.idx & fixrho.idx0, 3]])
        while(any(tau[-covariance.idx[, 1]] < 0) || any(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) > sqrt(tau[x[2]] * tau[x[3]])))) {
          Dtau <- Dtau / 2
          tau[idxtau] <- tau0[idxtau] + Dtau
          tau[-covariance.idx[, 1]][tau[-covariance.idx[, 1]] < tol & tau0[-covariance.idx[, 1]] < tol] <- 0
          if(any(fixrho != 0)) {
            idxrho <- which(fixrho != 0)
            tau[covariance.idx[idxrho, 1]] <- suppressWarnings(fixrho[idxrho] * sqrt(tau[covariance.idx[idxrho, 2]] * tau[covariance.idx[idxrho, 3]]))
          }
          fixrho.idx <- suppressWarnings(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) > (1 - 1.01 * tol) * sqrt(tau[x[2]] * tau[x[3]])))
          tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]] <- sign(tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]]) * sqrt(tau[covariance.idx[fixrho.idx & fixrho.idx0, 2]] * tau[covariance.idx[fixrho.idx & fixrho.idx0, 3]])
        }
        tau[-covariance.idx[, 1]][tau[-covariance.idx[, 1]] < tol] <- 0
        tau[covariance.idx[fixrho.idx, 1]] <- sign(tau[covariance.idx[fixrho.idx, 1]]) * sqrt(tau[covariance.idx[fixrho.idx, 2]] * tau[covariance.idx[fixrho.idx, 3]])
      }
    }
    cov <- as.matrix(fit$cov)
    alpha <- as.numeric(fit$alpha)
    eta <- as.numeric(fit$eta) + offset
    if(verbose) {
      cat("Variance component estimates:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    Y <- eta - offset + (y - mu)/mu.eta
    sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0))*family$variance(mu))
    if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
    if(max(abs(tau)) > tol^(-2)) {
      warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
      i <- maxiter
      break
    }
  }
  converged <- ifelse(i < maxiter, TRUE, FALSE)
  res <- y - mu
  res.var <- rep(1, n)
  for(i in 1:ng) res.var[group.idx[[i]]] <- tau[i]
  if(!is.Matrix) fit$Sigma_i <- fit$Sigma_iX <- NULL
  names(alpha) <- names(fit0$coef)
  rownames(cov) <- rownames(vcov(fit0))
  colnames(cov) <- colnames(vcov(fit0))
  return(list(theta=tau, n.pheno=1, n.groups=ng, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, X=X, P=fit$P, residuals=res, scaled.residuals=res*as.vector(weights(fit0))/res.var, cov=cov, Sigma_i=fit$Sigma_i, Sigma_iX=fit$Sigma_iX, converged=converged))
}

glmmkin.multi.ai <- function(fit0, kins, covariance.idx = NULL, group.idx, tau = rep(0, length(kins)), fixtau = rep(0, length(kins)), fixrho = NULL, maxiter = 500, tol = 1e-5, verbose = FALSE) {
  y <- fit0$model[[1]]
  n <- nrow(y)
  n.pheno <- ncol(y)
  offset <- fit0$offset
  if(is.null(offset)) offset <- rep(0, n)
  mu <- fit0$fitted.values
  Y <- y - offset
  Y2 <- as.vector(Y)
  X <- model.matrix(fit0)
  X2 <- Diagonal(n=n.pheno) %x% X
  alpha <- fit0$coef
  if(verbose) {
    cat("Fixed-effect coefficients:\n")
    print(alpha)
  }
  q <- length(kins)
  ng <- length(group.idx)*n.pheno*(n.pheno+1)/2
  idxtau <- which(fixtau == 0)
  idxtau2 <- intersect(covariance.idx[, 1], idxtau)
  q2 <- sum(fixtau == 0)
  if(q2 > 0) {
    tau[idxtau] <- rep(var(Y2)/q, q2)
    tau[idxtau2] <- 0
    Sigma <- 0*Diagonal(n = n*n.pheno)
    for(i in 1:q) {
      Sigma <- Sigma + tau[i]*kins[[i]]
    }
    Sigma <- forceSymmetric(Sigma)
    Sigma_i <- chol2inv(chol(Sigma))
    rm(Sigma)
    gc()
    Sigma_iX <- crossprod(Sigma_i, X2)
    XSigma_iX <- crossprod(X2, Sigma_iX)
    XSigma_iX <- forceSymmetric(XSigma_iX)
    cov <- chol2inv(chol(XSigma_iX))
    Sigma_iXcov <- tcrossprod(Sigma_iX, cov)
    PY <- crossprod(Sigma_i, Y2) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, Y2)))
    tau0 <- tau
    for(i in 1:q2) {
      APY <- crossprod(kins[[idxtau[i]]], PY)
      PAPY <- crossprod(Sigma_i, APY) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, APY)))
      tau[idxtau[i]] <- if(idxtau[i] %in% idxtau2) 0 else max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (sum(Y2 * PAPY) - (sum(Sigma_i*kins[[idxtau[i]]])-sum(Sigma_iX*crossprod(kins[[idxtau[i]]],Sigma_iXcov))))/n)
    }
  }
  for (i in seq_len(maxiter)) {
    if(verbose) cat("\nIteration ", i, ":\n")
    alpha0 <- alpha
    tau0 <- tau
    fit <- R_fitglmm_ai_noresidual(Y2, X2, q, kins, ng, tau, fixtau)
    if(q2 > 0) {
      Dtau <- as.numeric(fit$Dtau)
      tau[idxtau] <- tau0[idxtau] + Dtau
      if(is.null(covariance.idx)) {
        tau[tau < tol & tau0 < tol] <- 0
        while(any(tau < 0)) {
          Dtau <- Dtau / 2
          tau[idxtau] <- tau0[idxtau] + Dtau
          tau[tau < tol & tau0 < tol] <- 0
        }
        tau[tau < tol] <- 0
      } else {
        fixrho.idx0 <- apply(covariance.idx, 1, function(x) abs(tau0[x[1]]) > (1 - 1.01 * tol) * sqrt(tau0[x[2]] * tau0[x[3]]))
        tau[-covariance.idx[, 1]][tau[-covariance.idx[, 1]] < tol & tau0[-covariance.idx[, 1]] < tol] <- 0
        if(any(fixrho != 0)) {
          idxrho <- which(fixrho != 0)
          tau[covariance.idx[idxrho, 1]] <- suppressWarnings(fixrho[idxrho] * sqrt(tau[covariance.idx[idxrho, 2]] * tau[covariance.idx[idxrho, 3]]))
        }
        fixrho.idx <- suppressWarnings(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) > (1 - 1.01 * tol) * sqrt(tau[x[2]] * tau[x[3]])))
        tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]] <- sign(tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]]) * sqrt(tau[covariance.idx[fixrho.idx & fixrho.idx0, 2]] * tau[covariance.idx[fixrho.idx & fixrho.idx0, 3]])
        while(any(tau[-covariance.idx[, 1]] < 0) || any(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) > sqrt(tau[x[2]] * tau[x[3]])))) {
          Dtau <- Dtau / 2
          tau[idxtau] <- tau0[idxtau] + Dtau
          tau[-covariance.idx[, 1]][tau[-covariance.idx[, 1]] < tol & tau0[-covariance.idx[, 1]] < tol] <- 0
          if(any(fixrho != 0)) {
            idxrho <- which(fixrho != 0)
            tau[covariance.idx[idxrho, 1]] <- suppressWarnings(fixrho[idxrho] * sqrt(tau[covariance.idx[idxrho, 2]] * tau[covariance.idx[idxrho, 3]]))
          }
          fixrho.idx <- suppressWarnings(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) > (1 - 1.01 * tol) * sqrt(tau[x[2]] * tau[x[3]])))
          tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]] <- sign(tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]]) * sqrt(tau[covariance.idx[fixrho.idx & fixrho.idx0, 2]] * tau[covariance.idx[fixrho.idx & fixrho.idx0, 3]])
        }
        tau[-covariance.idx[, 1]][tau[-covariance.idx[, 1]] < tol] <- 0
        tau[covariance.idx[fixrho.idx, 1]] <- sign(tau[covariance.idx[fixrho.idx, 1]]) * sqrt(tau[covariance.idx[fixrho.idx, 2]] * tau[covariance.idx[fixrho.idx, 3]])
      }
    }
    cov <- as.matrix(fit$cov)
    alpha <- as.numeric(fit$alpha)
    eta <- as.numeric(fit$eta) + rep(offset, n.pheno)
    if(verbose) {
      cat("Variance component estimates:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
    if(max(abs(tau)) > tol^(-2)) {
      warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
      i <- maxiter
      break
    }
  }
  converged <- ifelse(i < maxiter, TRUE, FALSE)
  eta <- matrix(eta, ncol = n.pheno)
  colnames(eta) <- colnames(y)
  res <- y - eta
  alpha <- matrix(alpha, ncol = n.pheno)
  rownames(alpha) <- rownames(fit0$coef)
  colnames(alpha) <- colnames(fit0$coef)
  rownames(Sigma_i) <- colnames(Sigma_i) <- rep(rownames(X), n.pheno)
  rownames(cov) <- rownames(vcov(fit0))
  colnames(cov) <- colnames(vcov(fit0))
  return(list(theta=tau, n.pheno = n.pheno, n.groups=NULL, coefficients=alpha, linear.predictors=eta, fitted.values=eta, Y=Y, X=X, P=fit$P, residuals=res, scaled.residuals=NULL, cov=cov, Sigma_i=fit$Sigma_i, Sigma_iX=fit$Sigma_iX, converged=converged))
}

glmmkin.brent <- function(fit0, kins, method = "REML", tau = 1, fixtau = 0, maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE) {
  y <- fit0$y
  n <- length(y)
  offset <- fit0$offset
  if(is.null(offset)) offset <- rep(0, n)
  family <- fit0$family
  eta <- fit0$linear.predictors
  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu)/mu.eta
  sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0))*fit0$family$variance(mu))
  X <- model.matrix(fit0)
  alpha <- fit0$coef
  if(verbose) {
    cat("Fixed-effect coefficients:\n")
    print(alpha)
  }
  dispersion <- ifelse(family$family %in% c("poisson", "binomial"), "N", "Y")
  method <- ifelse(method=="REML", "R", "L")
  for (i in seq_len(maxiter)) {
    if(verbose) cat("\nIteration ", i, ":\n")
    alpha0 <- alpha
    tau0 <- tau
    fit <- .Call(C_fitglmm_brent, Y, X, kins, sqrtW, method, dispersion, tau, fixtau, tol, taumin, taumax, tauregion)
    tau <- as.numeric(fit$tau)
    cov <- as.matrix(fit$cov)
    alpha <- as.numeric(fit$alpha)
    eta <- as.numeric(fit$eta) + offset
    if(verbose) {
      cat("Variance component estimates:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    Y <- eta - offset + (y - mu)/mu.eta
    sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0))*family$variance(mu))
    if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
    if(family$family == "gaussian") break
  }
  converged <- ifelse(i < maxiter, TRUE, FALSE)
  res <- y - mu
  fit$eval <- as.numeric(fit$eval)
  P <- fit$U %*% (t(fit$U) * fit$eval) - fit$U %*% (fit$UtX * fit$eval) %*% (cov %*% t(fit$UtX) %*% (t(fit$U) * fit$eval))
  rm(fit)
  gc()
  if(dispersion=="N") {
    theta <- c(1, tau)
  } else {
    phi <- ifelse(method == "R", as.numeric(t(Y) %*% P %*% Y)/(n - ncol(X)), as.numeric(t(Y) %*% P %*% Y)/n)
    theta <- phi * c(1, tau)
    P <- P/phi
    cov <- phi*cov
  }
  names(alpha) <- names(fit0$coef)
  names(theta) <- c("dispersion", "kins1")
  rownames(cov) <- rownames(vcov(fit0))
  colnames(cov) <- colnames(vcov(fit0))
  return(list(theta=theta, n.pheno=1, n.groups=1, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, X=X, P=P, residuals=res, scaled.residuals=res*as.vector(weights(fit0))/theta[1], cov=cov, Sigma_i=NULL, Sigma_iX=NULL, converged=converged))
}

glmmkin.nm <- function(fit0, kins, method = "REML", tau = rep(1, length(kins)), fixtau = rep(0, length(kins)), maxiter = 500, tol = 1e-5, verbose = FALSE) {
  y <- fit0$y
  n <- length(y)
  offset <- fit0$offset
  if(is.null(offset)) offset <- rep(0, n)
  family <- fit0$family
  eta <- fit0$linear.predictors
  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu)/mu.eta
  sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0))*fit0$family$variance(mu))
  X <- model.matrix(fit0)
  alpha <- fit0$coef
  if(verbose) {
    cat("Fixed-effect coefficients:\n")
    print(alpha)
  }
  dispersion <- ifelse(family$family %in% c("poisson", "binomial"), "N", "Y")
  method <- ifelse(method=="REML", "R", "L")
  for (i in seq_len(maxiter)) {
    if(verbose) cat("\nIteration ", i, ":\n")
    alpha0 <- alpha
    tau0 <- tau
    fit <- .Call(C_fitglmm_nm, Y, X, length(kins), kins, sqrtW^2, method, dispersion, tau, fixtau, maxiter, tol)
    tau <- as.numeric(fit$tau)
    cov <- as.matrix(fit$cov)
    alpha <- as.numeric(fit$alpha)
    eta <- as.numeric(fit$eta) + offset
    if(verbose) {
      cat("Variance component estimates:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    Y <- eta - offset + (y - mu)/mu.eta
    sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0))*family$variance(mu))
    if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
    if(family$family == "gaussian") break
  }
  converged <- ifelse(i < maxiter, TRUE, FALSE)
  res <- y - mu
  fit$eval <- as.numeric(fit$eval)
  P <- fit$U %*% (t(fit$U) * fit$eval) - fit$U %*% (fit$UtX * fit$eval) %*% (cov %*% t(fit$UtX) %*% (t(fit$U) * fit$eval))
  rm(fit)
  gc()
  if(dispersion=="N") {
    theta <- c(1, tau)
  } else {
    phi <- ifelse(method == "R", as.numeric(t(Y) %*% P %*% Y)/(n - ncol(X)), as.numeric(t(Y) %*% P %*% Y)/n)
    theta <- phi * c(1, tau)
    P <- P/phi
    cov <- phi*cov
  }
  names(alpha) <- names(fit0$coef)
  names(theta) <- c("dispersion", paste0("kins", 1:length(kins)))
  rownames(cov) <- rownames(vcov(fit0))
  colnames(cov) <- colnames(vcov(fit0))
  return(list(theta=theta, n.pheno=1, n.groups=1, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, X=X, P=P, residuals=res, scaled.residuals=res*as.vector(weights(fit0))/theta[1], cov=cov, Sigma_i=NULL, Sigma_iX=NULL, converged=converged))
}

R_fitglmm_ai <- function(Y, X, q, kins, ng, group.idx, W, tau, fixtau) {
  n <- nrow(X)
  p <- ncol(X)
  q2 <- sum(fixtau == 0)
  diagSigma <- rep(0, n)
  for(i in 1:ng) diagSigma[group.idx[[i]]] <- tau[i]/W[group.idx[[i]]]
  Sigma <- Diagonal(x = diagSigma)
  for(i in 1:q) Sigma <- Sigma + tau[i+ng]*kins[[i]]
  #	Sigma <- forceSymmetric(Sigma)
  Sigma_i <- chol2inv(chol(Sigma))
  rm(Sigma)
  gc()
  Sigma_iX <- crossprod(Sigma_i, X)
  XSigma_iX <- crossprod(X, Sigma_iX)
  XSigma_iX <- forceSymmetric(XSigma_iX)
  cov <- chol2inv(chol(XSigma_iX))
  Sigma_iXcov <- tcrossprod(Sigma_iX, cov)
  alpha <- crossprod(cov, crossprod(Sigma_iX, Y))
  eta <- Y - diagSigma*(crossprod(Sigma_i,Y)-tcrossprod(Sigma_iX,t(alpha)))
  if(q2 > 0) {
    idxtau <- which(fixtau == 0)
    PY <- crossprod(Sigma_i, Y) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, Y)))
    wPY <- PY/W
    diagP <- (diag(Sigma_i) - rowSums(Sigma_iX * Sigma_iXcov))/W
    AI <- matrix(NA, q2, q2)
    score <- rep(NA, q2)
    for(i in 1:q2) {
      if(idxtau[i] <= ng) {
        score[i] <- sum(wPY[group.idx[[idxtau[i]]]] * PY[group.idx[[idxtau[i]]]] - diagP[group.idx[[idxtau[i]]]])
        for(j in 1:i) {
          AI[i,j] <- sum(wPY[group.idx[[idxtau[i]]]] * crossprod(Sigma_i[group.idx[[idxtau[j]]], group.idx[[idxtau[i]]]], wPY[group.idx[[idxtau[j]]]])) - sum(crossprod(Sigma_iXcov[group.idx[[idxtau[i]]], ], wPY[group.idx[[idxtau[i]]]]) * crossprod(Sigma_iX[group.idx[[idxtau[j]]], ], wPY[group.idx[[idxtau[j]]]]))
          if(j != i) AI[j,i] <- AI[i,j]
        }
      } else {
        APY <- crossprod(kins[[idxtau[i]-ng]], PY)
        PAPY <- crossprod(Sigma_i, APY) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, APY)))
        score[i] <- sum(Y * PAPY) - (sum(Sigma_i*kins[[idxtau[i]-ng]])-sum(Sigma_iX*crossprod(kins[[idxtau[i]-ng]],Sigma_iXcov)))
        for(j in 1:i) {
          if(idxtau[j] <= ng) AI[i,j] <- sum(wPY[group.idx[[idxtau[j]]]] * PAPY[group.idx[[idxtau[j]]]])
          else AI[i,j] <- sum(PY * crossprod(kins[[idxtau[j]-ng]], PAPY))
          if(j != i) AI[j,i] <- AI[i,j]
        }
      }
    }
    Dtau <- solve(AI, score)
    return(list(Dtau = Dtau, P = NULL, cov = cov, alpha = alpha, eta = eta, Sigma_i = Sigma_i, Sigma_iX = Sigma_iX))
  }
  return(list(Dtau = NULL, P = NULL, cov = cov, alpha = alpha, eta = eta, Sigma_i = Sigma_i, Sigma_iX = Sigma_iX))
}

R_fitglmm_ai_dense <- function(Y, X, q, kins, ng, group.idx, W, tau, fixtau) {
  n <- nrow(X)
  p <- ncol(X)
  q2 <- sum(fixtau == 0)
  diagSigma <- rep(0, n)
  for(i in 1:ng) diagSigma[group.idx[[i]]] <- tau[i]/W[group.idx[[i]]]
  #	Sigma <- Diagonal(x = diagSigma)
  Sigma <- diag(diagSigma)
  for(i in 1:q) Sigma <- Sigma + tau[i+ng]*kins[[i]]
  ##	Sigma <- forceSymmetric(Sigma)
  Sigma_i <- chol2inv(chol(Sigma))
  rm(Sigma)
  gc()
  Sigma_iX <- crossprod(Sigma_i, X)
  XSigma_iX <- crossprod(X, Sigma_iX)
  #	XSigma_iX <- forceSymmetric(XSigma_iX)
  cov <- chol2inv(chol(XSigma_iX))
  Sigma_iXcov <- tcrossprod(Sigma_iX, cov)
  P <- Sigma_i - tcrossprod(Sigma_iXcov, Sigma_iX)
  alpha <- crossprod(cov, crossprod(Sigma_iX, Y))
  eta <- Y - diagSigma*(crossprod(Sigma_i,Y)-tcrossprod(Sigma_iX,t(alpha)))
  rm(Sigma_i)
  gc()
  if(q2 > 0) {
    idxtau <- which(fixtau == 0)
    #		PY <- crossprod(Sigma_i, Y) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, Y)))
    PY <- crossprod(P, Y)
    wPY <- PY/W
    #		diagP <- (diag(Sigma_i) - rowSums(Sigma_iX * Sigma_iXcov))/W
    diagP <- diag(P)/W
    AI <- matrix(NA, q2, q2)
    score <- rep(NA, q2)
    for(i in 1:q2) {
      if(idxtau[i] <= ng) {
        score[i] <- sum(wPY[group.idx[[idxtau[i]]]] * PY[group.idx[[idxtau[i]]]] - diagP[group.idx[[idxtau[i]]]])
        for(j in 1:i) {
          #				      	AI[i,j] <- sum(wPY[group.idx[[idxtau[i]]]] * crossprod(Sigma_i[group.idx[[idxtau[j]]], group.idx[[idxtau[i]]]], wPY[group.idx[[idxtau[j]]]])) - sum(crossprod(Sigma_iXcov[group.idx[[idxtau[i]]], ], wPY[group.idx[[idxtau[i]]]]) * crossprod(Sigma_iX[group.idx[[idxtau[j]]], ], wPY[group.idx[[idxtau[j]]]]))
          AI[i,j] <- crossprod(wPY[group.idx[[idxtau[i]]]], crossprod(P[group.idx[[idxtau[j]]], group.idx[[idxtau[i]]]], wPY[group.idx[[idxtau[j]]]]))
          if(j != i) AI[j,i] <- AI[i,j]
        }
      } else {
        #				APY <- crossprod(kins[[idxtau[i]-ng]], PY)
        #	        	        PAPY <- crossprod(Sigma_i, APY) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, APY)))
        PAPY <- crossprod(P, crossprod(kins[[idxtau[i]-ng]], PY))
        #				score[i] <- sum(Y * PAPY) - (sum(Sigma_i*kins[[idxtau[i]-ng]])-sum(Sigma_iX*crossprod(kins[[idxtau[i]-ng]],Sigma_iXcov)))
        score[i] <- sum(Y * PAPY) - sum(P * kins[[idxtau[i]-ng]])
        for(j in 1:i) {
          if(idxtau[j] <= ng) AI[i,j] <- sum(wPY[group.idx[[idxtau[j]]]] * PAPY[group.idx[[idxtau[j]]]])
          else AI[i,j] <- sum(PY * crossprod(kins[[idxtau[j]-ng]], PAPY))
          if(j != i) AI[j,i] <- AI[i,j]
        }
      }
    }
    Dtau <- solve(AI, score)
    return(list(Dtau = Dtau, P = P, cov = cov, alpha = alpha, eta = eta))
  }
  return(list(Dtau = NULL, P = P, cov = cov, alpha = alpha, eta = eta))
}

R_fitglmm_ai_noresidual <- function(Y, X, q, kins, ng, tau, fixtau) {
  n <- nrow(X)
  p <- ncol(X)
  q2 <- sum(fixtau == 0)
  Sigma <- 0*Diagonal(n = n)
  for(i in 1:q) {
    Sigma <- Sigma + tau[i]*kins[[i]]
    if(i == ng) resSigma <- forceSymmetric(Sigma)
  }
  Sigma <- forceSymmetric(Sigma)
  Sigma_i <- chol2inv(chol(Sigma))
  rm(Sigma)
  gc()
  Sigma_iX <- crossprod(Sigma_i, X)
  XSigma_iX <- crossprod(X, Sigma_iX)
  XSigma_iX <- forceSymmetric(XSigma_iX)
  cov <- chol2inv(chol(XSigma_iX))
  Sigma_iXcov <- tcrossprod(Sigma_iX, cov)
  alpha <- crossprod(cov, crossprod(Sigma_iX, Y))
  eta <- Y - crossprod(resSigma, (crossprod(Sigma_i,Y)-tcrossprod(Sigma_iX,t(alpha))))
  rm(resSigma)
  gc()
  if(q2 > 0) {
    idxtau <- which(fixtau == 0)
    PY <- crossprod(Sigma_i, Y) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, Y)))
    AI <- matrix(NA, q2, q2)
    score <- rep(NA, q2)
    for(i in 1:q2) {
      APY <- crossprod(kins[[idxtau[i]]], PY)
      PAPY <- crossprod(Sigma_i, APY) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, APY)))
      score[i] <- sum(Y * PAPY) - (sum(Sigma_i*kins[[idxtau[i]]])-sum(Sigma_iX*crossprod(kins[[idxtau[i]]],Sigma_iXcov)))
      for(j in 1:i) {
        AI[i,j] <- sum(PY * crossprod(kins[[idxtau[j]]], PAPY))
        if(j != i) AI[j,i] <- AI[i,j]
      }
    }
    Dtau <- solve(AI, score)
    return(list(Dtau = Dtau, P = NULL, cov = cov, alpha = alpha, eta = eta, Sigma_i = Sigma_i, Sigma_iX = Sigma_iX))
  }
  return(list(Dtau = NULL, P = NULL, cov = cov, alpha = alpha, eta = eta, Sigma_i = Sigma_i, Sigma_iX = Sigma_iX))
}