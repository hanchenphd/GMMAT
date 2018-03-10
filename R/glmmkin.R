glmmkin <- function(fixed, data = parent.frame(), kins, groups = NULL, family = binomial(link = "logit"), method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE, ...) {
	call <- match.call()
	if(!class(kins) %in% c("matrix", "list"))
		stop("Error: \"kins\" must be a matrix or a list.")
	if(class(kins) == "list" && length(table(sapply(kins, dim))) != 1)
		stop("Error: when \"kins\" is a list, all its elements must be square matrices of the same size.")
	if(!method %in% c("REML", "ML"))
		stop("Error: \"method\" must be \"REML\" or \"ML\".")
	method.optim <- try(match.arg(method.optim, c("AI", "Brent", "Nelder-Mead")))
	if(class(method.optim) == "try-error")
		stop("Error: \"method.optim\" must be \"AI\", \"Brent\" or \"Nelder-Mead\".")
	if(method.optim == "AI" && method == "ML")
		stop("Error: method \"ML\" not available for method.optim \"AI\", use method \"REML\" instead.")
	if(method.optim == "Brent" && class(kins) == "list")
		stop("Error: method.optim \"Brent\" can only be applied in one-dimensional optimization, use a matrix for \"kins\".")
	if(class(family) != "family")
		stop("Error: \"family\" must be an object of class \"family\".")
	if(!family$family %in% c("binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"))
		stop("Error: \"family\" must be one of the following: binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, quasipoisson.")
	if(!is.null(groups)) {
		if(family$family != "gaussian") stop("Error: heteroscedastic linear mixed models are only applicable when \"family\" is gaussian.")
		if(method.optim != "AI") stop("Error: heteroscedastic linear mixed models are currently only implemented for method.optim \"AI\".")
		if(!groups %in% names(data)) stop("Error: \"groups\" must be one of the variables in the names of \"data\".")
	}
	if(method.optim != "Brent" && class(kins) == "matrix") kins <- list(kins1 = kins)
	fit0 <- glm(formula = fixed, data = data, family = family, ...)
	idx <- match(rownames(model.frame(formula = fixed, data = data, na.action = na.omit)), rownames(model.frame(formula = fixed, data = data, na.action = na.pass)))
	if(class(kins) == "matrix") kins <- kins[idx, idx]
	else {
	        for(i in 1:length(kins)) kins[[i]] <- kins[[i]][idx, idx]
	}
	group.id <- if(is.null(groups)) rep(1, length(idx)) else data[idx, groups]
	fit <- glmmkin.fit(fit0, kins, group.id, method = method, method.optim = method.optim, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
	fit$call <- call
	fit$id_include <- idx
	class(fit) <- "glmmkin"
	return(fit)
}

glmmkin.fit <- function(fit0, kins, group.id, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE) {
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
			fixtau.old <- rep(0, length(kins)+length(group.idx))
			fit <- glmmkin.ai(fit0, kins, group.idx, maxiter = maxiter, tol = tol, verbose = verbose)
			fixtau.new <- 1*(fit$theta < 1.01 * tol)
			while(any(fixtau.new != fixtau.old)) {
				warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
				fixtau.old <- fixtau.new
				fit <- glmmkin.ai(fit0, kins, group.idx, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
				fixtau.new <- 1*(fit$theta < 1.01 * tol)
			}
			if(!fit$converged) {
				if(length(group.idx) != 1) stop("Error: Average Information REML not converged, cannot refit heteroscedastic linear mixed model using Brent or Nelder-Mead methods.")
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

glmmkin.ai <- function(fit0, kins, group.idx, tau = rep(0, length(kins)+length(group.idx)), fixtau = rep(0, length(kins)+length(group.idx)), maxiter = 500, tol = 1e-5, verbose = FALSE) {
	y <- fit0$y
	n <- length(y)
	offset <- fit0$offset
	if(is.null(offset)) offset <- rep(0, n)
	family <- fit0$family
	eta <- fit0$linear.predictors
	mu <- fit0$fitted.values
	mu.eta <- family$mu.eta(eta)
	Y <- eta - offset + (y - mu)/mu.eta
	sqrtW <- mu.eta/sqrt(family$variance(mu))
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
	q2 <- sum(fixtau == 0)
	if(q2 > 0) {
	        tau[fixtau == 0] <- rep(var(Y)/(q+ng), q2)
		diagSigma <- rep(0, n)
		for(i in 1:ng) diagSigma[group.idx[[i]]] <- tau[i]/sqrtW[group.idx[[i]]]^2
		Sigma <- diag(diagSigma)
		for(i in 1:q) Sigma <- Sigma + tau[i+ng]*kins[[i]]
		Sigma_i <- chol2inv(chol(Sigma))
		rm(Sigma, diagSigma)
		gc()
		Sigma_iX <- crossprod(Sigma_i, X)
		P <- Sigma_i - tcrossprod(tcrossprod(Sigma_iX, chol2inv(chol(crossprod(X, Sigma_iX)))), Sigma_iX)
		rm(Sigma_i)
		gc()
		PY <- crossprod(P, Y)
		tau0 <- tau
		for(i in 1:q2) {
		        if(idxtau[i] <= ng) tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (sum((PY/sqrtW)[group.idx[[idxtau[i]]]]^2) - sum((diag(P)/sqrtW^2)[group.idx[[idxtau[i]]]]))/n)
			else {
	        	        PAPY <- crossprod(P, crossprod(kins[[idxtau[i]-ng]], PY))
				tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (crossprod(Y, PAPY) - sum(P*kins[[idxtau[i]-ng]]))/n)
			}
		}
		rm(P)
		gc()
	}
	for (i in seq_len(maxiter)) {
		if(verbose) cat("\nIteration ", i, ":\n")
		alpha0 <- alpha
		tau0 <- tau
		fit <- .Call(C_fitglmm_ai, Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau, tol)
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
		sqrtW <- mu.eta/sqrt(family$variance(mu))
		if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
		if(max(tau) > tol^(-2)) {
			warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
			i <- maxiter
			break
		}
	}
	converged <- ifelse(i < maxiter, TRUE, FALSE)
	res <- y - mu
	res.var <- rep(1, n)
	for(i in 1:ng) res.var[group.idx[[i]]] <- tau[i]
	return(list(theta=tau, n.groups=ng, coefficients=alpha,
	linear.predictors=eta, fitted.values=mu, Y=Y, P=fit$P, residuals=res,
	scaled.residuals=res/res.var, cov=cov, converged=converged))
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
	sqrtW <- mu.eta/sqrt(fit0$family$variance(mu))
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
		sqrtW <- mu.eta/sqrt(family$variance(mu))
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
	return(list(theta=theta, n.groups=1, coefficients=alpha,
	linear.predictors=eta, fitted.values=mu, Y=Y, P=P, residuals=res,
	scaled.residuals=res/theta[1], cov=cov, converged=converged))
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
	sqrtW <- mu.eta/sqrt(fit0$family$variance(mu))
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
		sqrtW <- mu.eta/sqrt(family$variance(mu))
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
	return(list(theta=theta, n.groups=1, coefficients=alpha,
	linear.predictors=eta, fitted.values=mu, Y=Y, P=P, residuals=res,
	scaled.residuals=res/theta[1], cov=cov, converged=converged))
}
