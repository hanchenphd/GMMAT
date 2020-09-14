glmm.wald <- function(fixed, data = parent.frame(), kins = NULL, id, random.slope = NULL, groups = NULL, family = binomial(link = "logit"), infile, snps, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, center = T, select = NULL, missing.method = "impute2mean", infile.nrow = NULL, infile.nrow.skip = 0, infile.sep = "\t", infile.na = "NA", snp.col = 1, infile.ncol.skip = 1, infile.ncol.print = 1, infile.header.print = "SNP", is.dosage = FALSE, verbose = FALSE, ...) {
	is.Windows <- Sys.info()["sysname"] == "Windows"
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
	miss.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
	if(class(miss.method) == "try-error") stop("Error: missing.method should be one of the following: impute2mean, omit!")
	miss.method <- substr(miss.method, 1, 1)
	is.plinkfiles <- all(file.exists(paste(infile, c("bim", "bed", "fam"), sep=".")))
	is.gds <- grepl("\\.gds$", infile)
	if(is.plinkfiles) {
		bimfile <- paste(infile, "bim", sep=".")
		bedfile <- paste(infile, "bed", sep=".")
		famfile <- paste(infile, "fam", sep=".")
		sample.id <- read.table(famfile, as.is=T)[,2]
		snpinfo <- matrix(NA, length(snps), 6)
	} else if(is.gds) { # GDS genotype file
	        gds <- SeqArray::seqOpen(infile)
		sample.id <- SeqArray::seqGetData(gds, "sample.id")
		variant.idx <- SeqArray::seqGetData(gds, "variant.id")
		variant.id <- SeqArray::seqGetData(gds, "annotation/id")
		snpinfo <- matrix(NA, length(snps), 5)
	} else { # text genotype files
		if(is.null(infile.nrow)) {
                        if(Sys.info()["sysname"] != "Linux") {
				infile.nrow <- length(readLines(infile))
                        } else {
                                if(grepl("\\.gz$", infile)) infile.nrow <- suppressWarnings(as.integer(system(paste("gzip -dc", infile, "| wc -l | gawk '{print $1}'"), intern = T)))
                                else if(grepl("\\.bz2$", infile)) infile.nrow <- suppressWarnings(as.integer(system(paste("bzip2 -dc", infile, "| wc -l | gawk '{print $1}'"), intern = T)))
                                else infile.nrow <- suppressWarnings(as.integer(system(paste("wc -l", infile, "| gawk '{print $1}'"), intern = T)))
                        }
                        if(any(is.na(infile.nrow))) infile.nrow <- length(readLines(infile))
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
		if(is.null(select)) {
			warning("Argument select is unspecified... Assuming the order of individuals in infile matches unique id in data...")
			select <- 1:length(unique(data[, id]))
		}
	}
	if(is.null(select)) {
		if(any(is.na(match(unique(data[, id]), sample.id)))) warning("Check your data... Some individuals in data are missing in sample.id of infile!")
		select <- match(sample.id, unique(data[, id]))
                select[is.na(select)] <- 0
		if(all(select == 0)) stop("Error: ID in data does not match sample.id in infile!")
	}
	if((is.plinkfiles || is.gds) && length(select) != length(sample.id)) stop("Error: number of individuals in select does not match infile!")
	select2 <- select[select > 0]
	if(any(duplicated(select2)) || max(select2) > length(unique(data[, id]))) stop("Error: select is a vector of orders, individuals not in the phenotype data should be coded 0!")
	data <- data[data[, id] %in% unique(data[, id])[select2], ]
	select[select > 0] <- rank(select[select > 0])
	idx <- match(rownames(model.frame(formula = fixed, data = data, na.action = na.omit)), rownames(model.frame(formula = fixed, data = data, na.action = na.pass)))
        group.id <- if(is.null(groups)) rep(1, length(idx)) else data[idx, groups]
        time.var <- if(is.null(random.slope)) NULL else data[idx, random.slope]
	nn <- length(unique(data[idx, id]))
	select2 <- match(unique(data[idx, id]), unique(data[, id]))
	select[!select %in% select2] <- 0
	select[select > 0] <- rank(select[select > 0])
	data <- data[idx, ]
	data.idx <- match(data[, id], unique(data[, id]))
        if(any(duplicated(data[, id]))) {
                cat("Duplicated id detected...\nAssuming longitudinal data with repeated measures...\n")
                if(method.optim == "Brent") {
                        if(is.null(kins)) {
                                kins <- diag(length(unique(data[, id])))
                                rownames(kins) <- colnames(kins) <- unique(data[, id])
                        } else stop("Error: method.optim \"Brent\" can only be applied to unrelated individuals in longitudinal data analysis.")
                } else {
                        if(method.optim != "AI") kins[[length(kins) + 1]] <- diag(length(unique(data[,id])))
                        else if(length(kins) > 0) kins[[length(kins) + 1]] <- as(diag(length(unique(data[, id]))), "sparseMatrix")
                        else kins <- list(kins1 = as(diag(length(unique(data[, id]))), "sparseMatrix"))
                        rownames(kins[[length(kins)]]) <- colnames(kins[[length(kins)]]) <- unique(data[, id])
                }
        } else if(!is.null(random.slope)) stop("Error: \"random.slope\" must be used for longitudinal data with duplicated \"id\".")
	if(class(kins)[1] == "matrix") {
		match.idx1 <- match(data[, id], rownames(kins))
		match.idx2 <- match(data[, id], colnames(kins))
		if(any(is.na(c(match.idx1, match.idx2)))) stop("Error: kins matrix does not include all individuals in the data.")
		kins <- kins[match.idx1, match.idx2]
	} else if(class(kins)[1] == "list") {
	        for(i in 1:length(kins)) {
		      	match.idx1 <- match(data[, id], rownames(kins[[i]]))
			match.idx2 <- match(data[, id], colnames(kins[[i]]))
			if(any(is.na(c(match.idx1, match.idx2)))) stop("Error: kins matrix ", i, " does not include all individuals in the data.")
			kins[[i]] <- kins[[i]][match.idx1, match.idx2]
		}
        } else {
                if(!is.null(groups)) stop("Error: heteroscedastic linear models for unrelated observations have not been implemented.")
	}
	N <- AF <- BETA <- SE <- PVAL <- converged <- rep(NA, length(snps))
	if(verbose) {
		if(is.Windows) pb <- winProgressBar(min = 0, max = length(snps))
		else {
	                cat("Progress of Wald test:\n")
			pb <- txtProgressBar(min = 0, max = length(snps), style = 3)
			cat("\n")
		}
	}
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
				geno <- if(is.dosage) SeqVarTools::imputedDosage(gds, use.names = FALSE) else SeqVarTools::altDosage(gds, use.names = FALSE)
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
				data$SNP__ <- as.numeric(readfile$G)[data.idx]
				data$SNP__[data$SNP__ < (-999)] <- NA
				fit0 <- do.call("glm", list(formula = as.formula(paste(paste(deparse(fixed), collapse = ""), "SNP__", sep=" + ")), data = data, family = family, ...))
				if(is.null(kins)) {
					coef <- summary(fit0)$coef
                                        BETA[ii] <- coef[nrow(coef), 1]
                                        SE[ii] <- coef[nrow(coef), 2]
                                        PVAL[ii] <- coef[nrow(coef), 4]
					converged[ii] <- TRUE
					next
				}
				idx <- match(rownames(model.frame(formula = as.formula(paste(paste(deparse(fixed), collapse = ""), "SNP__", sep=" + ")), data = data, na.action = na.omit)), rownames(model.frame(formula = as.formula(paste(paste(deparse(fixed), collapse = ""), "SNP__", sep=" + ")), data = data, na.action = na.pass)))
				tmpkins <- kins
				if(class(tmpkins)[1] == "matrix") tmpkins <- tmpkins[idx, idx]
				else {
				        for(i in 1:length(tmpkins)) tmpkins[[i]] <- tmpkins[[i]][idx, idx]
				}
				if(!is.null(time.var)) time.var <- time.var[idx]
				group.id <- group.id[idx]
				fit <- try(glmmkin.fit(fit0, tmpkins, time.var, group.id, method = method, method.optim = method.optim, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose))
				if(class(fit) != "try-error") {
					BETA[ii] <- fit$coefficients[length(fit$coefficients)]
					SE[ii] <- sqrt(diag(fit$cov)[length(fit$coefficients)])
					PVAL[ii] <- pchisq((BETA[ii]/SE[ii])^2, 1, lower.tail=F)
					converged[ii] <- fit$converged
				}
			}
		}
		if(verbose) {
		        if(is.Windows) setWinProgressBar(pb, ii, title=paste0("Progress of Wald test: ",round(ii/length(snps)*100),"%"))
			else setTxtProgressBar(pb, ii)
		}
	}
	if(verbose) close(pb)
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
