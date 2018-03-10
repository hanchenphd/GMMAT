glmm.wald <- function(fixed, data = parent.frame(), kins, groups = NULL, family = binomial(link = "logit"), infile, snps, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, center = T, select = NULL, missing.method = "impute2mean", infile.nrow = NULL, infile.nrow.skip = 0, infile.sep = "\t", infile.na = "NA", snp.col = 1, infile.ncol.skip = 1, infile.ncol.print = 1, infile.header.print = "SNP", verbose = FALSE, ...) {
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
		stop("Error: family must be one of the following: binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, quasipoisson.")
        if(!is.null(groups)) {
                if(family$family != "gaussian") stop("Error: heteroscedastic linear mixed models are only applicable when \"family\" is gaussian.")
                if(method.optim != "AI") stop("Error: heteroscedastic linear mixed models are currently only implemented for method.optim \"AI\".")
                if(!groups %in% names(data)) stop("Error: \"groups\" must be one of the variables in the names of \"data\".")
        }
        if(method.optim != "Brent" && class(kins) == "matrix") kins <- list(kins1 = kins)
	miss.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
	if(class(miss.method) == "try-error") stop("Error: missing.method should be one of the following: impute2mean, omit!")
	miss.method <- substr(miss.method, 1, 1)
	nn <- ifelse(class(kins) == "matrix", nrow(kins), nrow(kins[[1]]))
	if(is.null(select)) select <- 1:nn
	select2 <- select[select > 0]
	if(length(select2) != nn | any(sort(select2) != 1:nn)) stop("Error: select is a vector of orders, individuals not in kins should be coded 0!")
	idx <- match(rownames(model.frame(formula = fixed, data = data, na.action = na.omit)), rownames(model.frame(formula = fixed, data = data, na.action = na.pass)))
	nn <- length(idx)
	data <- data[idx, ]
        group.id <- if(is.null(groups)) rep(1, nn) else data[, groups]
	select2 <- match(idx, select)
	select[-select2] <- 0
	select[select2] <- 1:nn
	if(class(kins) == "matrix") kins <- kins[idx, idx]
	else {
	        for(i in 1:length(kins)) kins[[i]] <- kins[[i]][idx, idx]
	}
	is.plinkfiles <- all(file.exists(paste(infile, c("bim", "bed", "fam"), sep=".")))
	is.gds <- grepl("\\.gds$", infile)
	if(is.plinkfiles) {
		bimfile <- paste(infile, "bim", sep=".")
		bedfile <- paste(infile, "bed", sep=".")
		famfile <- paste(infile, "fam", sep=".")
		if(length(select) != as.integer(system(paste("wc -l", famfile, "| awk '{print $1}'"), intern = T))) stop("Error: number of individuals in plink fam file incorrect!")
		snpinfo <- matrix(NA, length(snps), 6)
	} else if(is.gds) { # GDS genotype file
	        gds <- SeqArray::seqOpen(infile)
		variant.idx <- SeqArray::seqGetData(gds, "variant.id")
		variant.id <- SeqArray::seqGetData(gds, "annotation/id")
		snpinfo <- matrix(NA, length(snps), 5)
	} else { # text genotype files
		if(is.null(infile.nrow)) {
			if(grepl("\\.gz$", infile)) infile.nrow <- as.integer(system(paste("zcat", infile, "| wc -l | awk '{print $1}'"), intern = T))
			else if(grepl("\\.bz2$", infile)) infile.nrow <- as.integer(system(paste("bzcat", infile, "| wc -l | awk '{print $1}'"), intern = T))
			else infile.nrow <- as.integer(system(paste("wc -l", infile, "| awk '{print $1}'"), intern = T))
		}
		if(!is.numeric(infile.nrow) | infile.nrow < 0)
			stop("Error: number of rows of the input file is incorrect!")
		if(!is.numeric(infile.nrow.skip) | infile.nrow.skip < 0)
		       	stop("Error: number of skipped rows of the input file is incorrect!")
		if(!is.numeric(infile.ncol.skip) | infile.ncol.skip <= 0)
			stop("Error: number of skipped cols of the input file is incorrect!")
		if(length(infile.ncol.print) != length(infile.header.print))
			stop("Error: number of cols selected to print does not match number of header names!")
		if(!is.numeric(infile.ncol.print) | (!snp.col %in% infile.ncol.print))
			stop("Error: snp.col is not in infile.ncol.print!")
		if(any(!is.numeric(infile.ncol.print)) | any(infile.ncol.print < 0) | any(infile.ncol.print > infile.ncol.skip))
			stop("Error: cols selected to print have incorrect indices!")
		if(any(infile.ncol.print != sort(infile.ncol.print)))
			stop("Error: col indices must be sorted increasingly in infile.ncol.print!")
		snpinfo <- matrix(NA, length(snps), length(infile.header.print))
	}
	N <- AF <- BETA <- SE <- PVAL <- converged <- rep(NA, length(snps))
	for(ii in 1:length(snps)) {
	        snp <- snps[ii]
		if(verbose) cat("\nAnalyze SNP ", ii, ": ", snp, "\n")
		if(is.plinkfiles) {
			if(center) {
				readfile <- .Call(C_glmm_wald_bed, nn, snp, bimfile, bedfile, 'c', miss.method, select)
			} else {
				readfile <- .Call(C_glmm_wald_bed, nn, snp, bimfile, bedfile, 'n', miss.method, select)
			}
		} else if(is.gds) { # GDS genotype file
			if(!snp %in% variant.id) readfile <- list(skip=2)
			else {
				SeqArray::seqSetFilter(gds, variant.id = variant.idx[variant.id == snp], verbose = FALSE)
				alleles <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
				readfile <- list(snpinfo = c(snp, SeqArray::seqGetData(gds, "chromosome"), SeqArray::seqGetData(gds, "position"), alleles[[1]][1], paste(alleles[[1]][-1], collapse=",")))
				geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
				geno <- as.numeric(geno)[select>0][order(select[select>0])]
				readfile$N <- sum(!is.na(geno))
				readfile$AF <- mean(geno, na.rm = TRUE)/2
				readfile$skip <- ifelse(max(geno, na.rm = TRUE)-min(geno, na.rm = TRUE)<tol, 1, 0)
				if(sum(is.na(geno))>0 & missing.method == "impute2mean") geno[is.na(geno)] <- readfile$AF * 2
				if(center) geno <- geno - readfile$AF * 2
				readfile$G <- geno
			}
		} else { # text genotype files
			if(center) {
				readfile <- .Call(C_glmm_wald_text, nn, snp, infile, tol, 'c', miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, snp.col, select)
			} else {
				readfile <- .Call(C_glmm_wald_text, nn, snp, infile, tol, 'n', miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, snp.col, select)
			}
		}
		if(readfile$skip == 2) { # snp not found in the file
			if(is.plinkfiles) {
				snpinfo[ii, 2] <- snp
			} else if(is.gds) { # GDS genotype file
			        snpinfo[ii, 1] <- snp
			} else { # text genotype files
				snpinfo[ii, which(infile.ncol.print == snp.col)] <- snp
			}
		} else {
			snpinfo[ii, ] <- readfile$snpinfo
			N[ii] <- readfile$N
			AF[ii] <- readfile$AF
			if(readfile$skip != 1) { # snp
				data$SNP__ <- as.numeric(readfile$G)
				data$SNP__[data$SNP__ < (-999)] <- NA
				fit0 <- glm(formula = as.formula(paste(deparse(fixed), "SNP__", sep=" + ")), data = data, family = family, ...)
				idx <- match(rownames(model.frame(formula = as.formula(paste(deparse(fixed), "SNP__", sep=" + ")), data = data, na.action = na.omit)), rownames(model.frame(formula = as.formula(paste(deparse(fixed), "SNP__", sep=" + ")), data = data, na.action = na.pass)))
				tmpkins <- kins
				if(class(tmpkins) == "matrix") tmpkins <- tmpkins[idx, idx]
				else {
				        for(i in 1:length(tmpkins)) tmpkins[[i]] <- tmpkins[[i]][idx, idx]
				}
				fit <- try(glmmkin.fit(fit0, tmpkins, group.id, method = method, method.optim = method.optim, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose))
				if(class(fit) != "try-error") {
					BETA[ii] <- fit$coefficients[length(fit$coefficients)]
					SE[ii] <- sqrt(diag(fit$cov)[length(fit$coefficients)])
					PVAL[ii] <- pchisq((BETA[ii]/SE[ii])^2, 1, lower.tail=F)
					converged[ii] <- fit$converged
				}
			}
		}
	}
	res <- data.frame(snpinfo, N, AF, BETA, SE, PVAL, converged)
	if(is.plinkfiles) {
		names(res)[1:6] <- c("CHR", "SNP", "cM", "POS", "A1", "A2")
	} else if(is.gds) { # GDS genotype file
	        names(res)[1:5] <- c("SNP", "CHR", "POS", "REF", "ALT")
		SeqArray::seqClose(gds)
	} else { # text genotype files
		names(res)[1:length(infile.header.print)] <- infile.header.print
	}
	return(res)
}
