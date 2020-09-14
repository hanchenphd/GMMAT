SMMAT <- function(null.obj, geno.file, group.file, group.file.sep = "\t", meta.file.prefix = NULL, MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(1, 25), miss.cutoff = 1, missing.method = "impute2mean", method = "davies", tests = "E", rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1), use.minor.allele = FALSE, auto.flip = FALSE, Garbage.Collection = FALSE, is.dosage = FALSE, ncores = 1, verbose = FALSE)
{
    is.Windows <- Sys.info()["sysname"] == "Windows"
    if(is.Windows && ncores > 1) {
        warning("The package doMC is not available on Windows... Switching to single thread...")
        ncores <- 1
    }
    if(!class(null.obj) %in% c("glmmkin", "glmmkin.multi")) stop("Error: null.obj must be a class glmmkin or glmmkin.multi object!")
    n.pheno <- null.obj$n.pheno
    missing.method <- try(match.arg(missing.method, c("impute2mean", "impute2zero")))
    if(class(missing.method) == "try-error") stop("Error: \"missing.method\" must be \"impute2mean\" or \"impute2zero\".")
    if(any(!tests %in% c("B", "S", "O", "E"))) stop("Error: \"tests\" should only include \"B\" for the burden test, \"S\" for SKAT, \"O\" for SKAT-O or \"E\" for the efficient hybrid test of the burden test and SKAT.")
    Burden <- "B" %in% tests
    SKAT <- "S" %in% tests
    SKATO <- "O" %in% tests
    SMMAT <- "E" %in% tests
    if(any(duplicated(null.obj$id_include))) {
        J <- Matrix(sapply(unique(null.obj$id_include), function(x) 1*(null.obj$id_include==x)), sparse = TRUE)
        residuals <- as.numeric(as.matrix(crossprod(J, null.obj$scaled.residuals)))
        if(!is.null(null.obj$P)) null.obj$P <- as.matrix(crossprod(J, crossprod(null.obj$P, J)))
	else {
	    null.obj$Sigma_iX <- crossprod(J, null.obj$Sigma_iX)
	    null.obj$Sigma_i <- forceSymmetric(crossprod(J,crossprod(null.obj$Sigma_i,J)))
	    null.obj$Sigma_i <- Matrix(null.obj$Sigma_i, sparse = TRUE)
	}
        rm(J)
    } else residuals <- null.obj$scaled.residuals
    n <- length(unique(null.obj$id_include))
    if(!grepl("\\.gds$", geno.file)) stop("Error: currently only .gds format is supported in geno.file!")
    gds <- SeqArray::seqOpen(geno.file)
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    if(any(is.na(match(null.obj$id_include, sample.id)))) warning("Check your data... Some individuals in null.obj$id_include are missing in sample.id of geno.file!")
    sample.id <- sample.id[sample.id %in% null.obj$id_include]
    if(length(sample.id) == 0) stop("Error: null.obj$id_include does not match sample.id in geno.file!")
    match.id <- match(sample.id, unique(null.obj$id_include))
    if(class(null.obj) == "glmmkin.multi") {
        residuals <- residuals[match.id, , drop = FALSE]
    	match.id <- rep(match.id, n.pheno) + rep((0:(n.pheno-1))*n, each = length(match.id))
    } else {
        residuals <- residuals[match.id]
    }
    if(!is.null(null.obj$P)) null.obj$P <- null.obj$P[match.id, match.id]
    else {
    	null.obj$Sigma_iX <- null.obj$Sigma_iX[match.id, , drop = FALSE]
	null.obj$Sigma_i <- null.obj$Sigma_i[match.id, match.id]
    }
    variant.idx <- SeqArray::seqGetData(gds, "variant.id")
    #variant.id <- SeqArray::seqGetData(gds, "annotation/id")
    #variant.id <- paste(SeqArray::seqGetData(gds, "chromosome"), SeqArray::seqGetData(gds, "position"), SeqArray::seqGetData(gds, "allele"), sep=":")
    chr <- SeqArray::seqGetData(gds, "chromosome")
    pos <- SeqArray::seqGetData(gds, "position")
    alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
    ref <- unlist(lapply(alleles.list, function(x) x[1]))
    alt <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
    rm(alleles.list); gc()
    SeqArray::seqClose(gds)
    variant.id <- paste(chr, pos, ref, alt, sep = ":")
    rm(chr, pos, ref, alt); gc()
    group.info <- try(read.table(group.file, header = FALSE, col.names = c("group", "chr", "pos", "ref", "alt", "weight"), colClasses = c("character","character","integer","character","character","numeric"), sep = group.file.sep), silent = TRUE)
    if (class(group.info) == "try-error") {
        stop("Error: cannot read group.file!")
    }
    variant.id1 <- paste(group.info$chr, group.info$pos, group.info$ref, group.info$alt, sep = ":")
    group.info <- group.info[!duplicated(paste(group.info$group, variant.id1, sep = ":")), ]
    variant.idx1 <- variant.idx[match(variant.id1, variant.id)]
    group.info$variant.idx <- variant.idx1
    group.info$flip <- 0
    if(auto.flip) {
        cat("Automatic allele flipping enabled...\nVariants matching alt/ref but not ref/alt alleles will also be included, with flipped effects\n")
    	variant.id2 <- paste(group.info$chr, group.info$pos, group.info$alt, group.info$ref, sep = ":")
    	variant.idx2 <- variant.idx[match(variant.id2, variant.id)]
    	if(any(!is.na(variant.idx1) & !is.na(variant.idx2))) {
            tmp.dups <- which(!is.na(variant.idx1) & !is.na(variant.idx2))
	    cat("The following ambiguous variants were found:\n")
	    cat("chr:", chr[tmp.dups], "\n")
	    cat("pos:", pos[tmp.dups], "\n")
	    cat("ref:", ref[tmp.dups], "\n")
	    cat("alt:", alt[tmp.dups], "\n")
	    cat("Warning: both variants with alleles ref/alt and alt/ref were present at the same position and coding should be double checked!\nFor these variants, only those with alleles ref/alt were used in the analysis...\n")
	    variant.idx2[tmp.dups] <- NA
	    rm(tmp.dups)
    	}
    	group.info$flip <- 1*(!is.na(variant.idx2))
    	group.info$variant.idx[!is.na(variant.idx2)] <- variant.idx2[!is.na(variant.idx2)]
	rm(variant.id2, variant.idx2)
    }
    rm(variant.id, variant.id1, variant.idx1); gc()
    #group.info$variant.idx <- variant.idx[match(group.info$variant, variant.id)]
    group.info <- subset(group.info, !is.na(variant.idx))
    groups <- unique(group.info$group)
    n.groups.all <- length(groups)
    group.info$group.idx <- as.numeric(factor(group.info$group))
    group.info <- group.info[order(group.info$group.idx, group.info$variant.idx), ]
    group.idx.end <- findInterval(1:n.groups.all, group.info$group.idx)
    group.idx.start <- c(1, group.idx.end[-n.groups.all] + 1)
    ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
    if(ncores > 1) {
    	doMC::registerDoMC(cores = ncores)
    	n.groups.percore <- (n.groups.all-1) %/% ncores + 1
    	n.groups.percore_1 <- n.groups.percore * ncores - n.groups.all
	b <- NULL
    	out <- foreach(b=1:ncores, .combine=rbind, .multicombine = TRUE, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
    	    idx <- if(b <= n.groups.percore_1) ((b-1)*(n.groups.percore-1)+1):(b*(n.groups.percore-1)) else (n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1-1)*n.groups.percore+1):(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1)*n.groups.percore)
	    n.groups <- length(idx)
	    if(verbose) {
	        if(b==1) cat("Progress of SMMAT:\n")
		pb <- txtProgressBar(min = 0, max = n.groups, style = 3)
		cat("\n")
	    }
    	    gds <- SeqArray::seqOpen(geno.file)
    	    SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
	    n.variants <- rep(0,n.groups)
    	    miss.min <- rep(NA,n.groups)
    	    miss.mean <- rep(NA, n.groups)
    	    miss.max <- rep(NA, n.groups)
    	    freq.min <- rep(NA, n.groups)
    	    freq.mean <- rep(NA, n.groups)
    	    freq.max <- rep(NA, n.groups)
	    if(Burden | SKATO | SMMAT) {
	        Burden.score <- rep(NA, n.groups)
    	    	Burden.var <- rep(NA, n.groups)
    	    	Burden.pval <- rep(NA, n.groups)
	    }
    	    if(SKAT | SKATO) SKAT.pval <- rep(NA, n.groups)
	    if(SKATO) {
	        SKATO.pval <- rep(NA, n.groups)
		SKATO.minp <- rep(NA, n.groups)
		SKATO.minp.rho <- rep(NA, n.groups)
	    }
    	    if(SMMAT) SMMAT.pval <- rep(NA, n.groups)
	    if(!is.null(meta.file.prefix)) {
	    	if(class(null.obj) == "glmmkin.multi") stop("Error: meta-analysis not supported yet for multiple phenotypes.")
	        if(.Platform$endian!="little") stop("Error: platform must be little endian.")
		meta.file.score <- paste0(meta.file.prefix, ".score.", b)
		meta.file.var <- paste0(meta.file.prefix, ".var.", b)
		write.table(t(c("group", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE", "VAR", "PVAL")), meta.file.score, quote = FALSE, row.names = FALSE, col.names = FALSE)
		meta.file.var.handle <- file(meta.file.var, "wb")
	    }
    	    for(i in 1:n.groups) {
	        if(verbose && i %% ceiling(n.groups/100) == 0) setTxtProgressBar(pb, i)
    	    	tmp.idx <- group.idx.start[idx[i]]:group.idx.end[idx[i]]
		tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
	    	SeqArray::seqSetFilter(gds, variant.id = tmp.group.info$variant.idx, verbose = FALSE)
                geno <- if(is.dosage) SeqVarTools::imputedDosage(gds, use.names = FALSE) else SeqVarTools::altDosage(gds, use.names = FALSE)
                miss <- colMeans(is.na(geno))
                freq <- colMeans(geno, na.rm = TRUE)/2
	    	include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
	    	n.p <- sum(include)
	    	if(n.p == 0) next
		tmp.group.info <- tmp.group.info[include, , drop = FALSE]
	    	miss <- miss[include]
	    	freq <- freq[include]
                geno <- geno[, include, drop = FALSE]
		N <- nrow(geno) - colSums(is.na(geno))
		if(sum(tmp.group.info$flip) > 0) {
		    freq[tmp.group.info$flip==1] <- 1 - freq[tmp.group.info$flip==1]
		    geno[, tmp.group.info$flip==1] <- 2 - geno[, tmp.group.info$flip==1]
		}
	    	if(max(miss)>0) {
	    	    miss.idx <- which(is.na(geno))
	    	    geno[miss.idx] <- if(missing.method=="impute2mean") 2*freq[ceiling(miss.idx/nrow(geno))] else 0
	    	}
	    	U <- as.vector(crossprod(geno, residuals))
		if(class(null.obj) == "glmmkin.multi") geno <- Diagonal(n = n.pheno) %x% geno
		if(!is.null(null.obj$P)) V <- crossprod(geno, crossprod(null.obj$P, geno))
		else {
		    GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
		    V <- crossprod(geno, crossprod(null.obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
		}
	    	V <- as.matrix(V)
		if(!is.null(meta.file.prefix)) {
		    VAR <- diag(V)
		    PVAL <- ifelse(VAR>0, pchisq(U^2/VAR, df=1, lower.tail=FALSE), NA)
		    write.table(cbind(tmp.group.info[,c("group","chr","pos","ref","alt")], N, miss, freq, U, VAR, PVAL), meta.file.score, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
		    writeBin(V[lower.tri(V, diag = TRUE)], meta.file.var.handle, size = 4)
		}
		if(use.minor.allele) {
		    tmp.group.info$weight[freq > 0.5] <- -tmp.group.info$weight[freq > 0.5]
		    freq[freq > 0.5] <- 1 - freq[freq > 0.5]
		}
		weights <- rep(tmp.group.info$weight * MAF.weights.beta.fun(freq, MAF.weights.beta[1], MAF.weights.beta[2]), n.pheno)
	    	n.variants[i] <- n.p
	    	miss.min[i] <- min(miss)
	    	miss.mean[i] <- mean(miss)
	    	miss.max[i] <- max(miss)
	    	freq.min[i] <- min(freq)
	    	freq.mean[i] <- mean(freq)
	    	freq.max[i] <- max(freq)
		U <- U*weights
		V <- t(V*weights)*weights
		if(max(V)-min(V) < sqrt(.Machine$double.eps)) {
	    	    burden.score <- sum(U)
    	    	    burden.var <- sum(V)
    	    	    burden.pval <- pchisq(burden.score^2/burden.var, df=1, lower.tail=FALSE)
	    	    if(Burden | SKATO | SMMAT) {
		       	Burden.score[i] <- burden.score
	    		Burden.var[i] <- burden.var
			Burden.pval[i] <- burden.pval
	    	    }
    	    	    if(SKAT | SKATO) SKAT.pval[i] <- burden.pval
	    	    if(SKATO) {
	                SKATO.pval[i] <- burden.pval
			SKATO.minp[i] <- burden.pval
			SKATO.minp.rho[i] <- 1
	    	    }
    	    	    if(SMMAT) SMMAT.pval[i] <- burden.pval
		} else {
		    if(SKATO) {
		        re <- .skato_pval(U = U, V = V, rho = rho, method = method)
			Burden.score[i] <- re$Burden.score
			Burden.var[i] <- re$Burden.var
			Burden.pval[i] <- re$Burden.pval
			SKAT.pval[i] <- re$SKAT.pval
			SKATO.pval[i] <-re$p
			SKATO.minp[i] <- re$minp
			SKATO.minp.rho[i] <- re$minp.rho
		    } else {
		        if(SKAT) SKAT.pval[i] <- .quad_pval(U = U, V = V, method = method)
		    	if(Burden | SMMAT) {
	    	    	    Burden.score[i] <- sum(U)
    	    	    	    Burden.var[i] <- sum(V)
    	    	    	    Burden.pval[i] <- pchisq(Burden.score[i]^2/Burden.var[i], df=1, lower.tail=FALSE)
		    	}
		    }	    
		    if(SMMAT) {
	    	        V.rowSums <- rowSums(V)
    	    		U <- U - V.rowSums * Burden.score[i] / Burden.var[i]
    	    		V <- V - tcrossprod(V.rowSums) / Burden.var[i]
    	    		if(mean(abs(V)) < sqrt(.Machine$double.eps)) SMMAT.pval[i] <- Burden.pval[i]
	    		else SMMAT.pval[i] <- tryCatch(pchisq(-2*log(Burden.pval[i])-2*log(.quad_pval(U = U, V = V, method = method)), df = 4, lower.tail = FALSE), error = function(e) { Burden.pval[i] })
		    }
		}
	    	rm(geno)
	    	if(Garbage.Collection) gc()
    	    }
    	    SeqArray::seqClose(gds)
	    if(verbose) {
	        setTxtProgressBar(pb, n.groups)
		close(pb)
	    }
	    if(!is.null(meta.file.prefix)) close(meta.file.var.handle)
    	    tmp.out <- data.frame(group=unique(group.info$group)[idx], n.variants=n.variants, miss.min=miss.min, miss.mean=miss.mean, miss.max=miss.max, freq.min=freq.min, freq.mean=freq.mean, freq.max=freq.max)
	    if(Burden | SKATO | SMMAT) {
	        tmp.out$B.score <- Burden.score
		tmp.out$B.var <- Burden.var
		tmp.out$B.pval <- Burden.pval
	    }
    	    if(SKAT | SKATO) tmp.out$S.pval <- SKAT.pval
	    if(SKATO) {
	        tmp.out$O.pval <- SKATO.pval
		tmp.out$O.minp <- SKATO.minp
		tmp.out$O.minp.rho <- SKATO.minp.rho
	    }
    	    if(SMMAT) tmp.out$E.pval <- SMMAT.pval
	    tmp.out
	}
    } else { # use a single core
	n.groups <- n.groups.all
	if(verbose) {
	    if(is.Windows) pb <- winProgressBar(min = 0, max = n.groups)
	    else {
	        cat("Progress of SMMAT:\n")
	    	pb <- txtProgressBar(min = 0, max = n.groups, style = 3)
		cat("\n")
	    }
	}
    	gds <- SeqArray::seqOpen(geno.file)
    	SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
	n.variants <- rep(0,n.groups)
    	miss.min <- rep(NA,n.groups)
    	miss.mean <- rep(NA, n.groups)
    	miss.max <- rep(NA, n.groups)
    	freq.min <- rep(NA, n.groups)
    	freq.mean <- rep(NA, n.groups)
    	freq.max <- rep(NA, n.groups)
	if(Burden | SKATO | SMMAT) {
	    Burden.score <- rep(NA, n.groups)
    	    Burden.var <- rep(NA, n.groups)
    	    Burden.pval <- rep(NA, n.groups)
	}
    	if(SKAT | SKATO) SKAT.pval <- rep(NA, n.groups)
	if(SKATO) {
	    SKATO.pval <- rep(NA, n.groups)
	    SKATO.minp <- rep(NA, n.groups)
	    SKATO.minp.rho <- rep(NA, n.groups)
	}
    	if(SMMAT) SMMAT.pval <- rep(NA, n.groups)
    	if(!is.null(meta.file.prefix)) {
	    if(class(null.obj) == "glmmkin.multi") stop("Error: meta-analysis not supported yet for multiple phenotypes.")
            if(.Platform$endian!="little") stop("Error: platform must be little endian.")
	    meta.file.score <- paste0(meta.file.prefix, ".score.1")
	    meta.file.var <- paste0(meta.file.prefix, ".var.1")
	    write.table(t(c("group", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE", "VAR", "PVAL")), meta.file.score, quote = FALSE, row.names = FALSE, col.names = FALSE)
	    meta.file.var.handle <- file(meta.file.var, "wb")
	}
    	for(i in 1:n.groups) {
	    if(verbose && i %% ceiling(n.groups/100) == 0) {
		if(is.Windows) setWinProgressBar(pb, i, title=paste0("Progress of SMMAT: ",round(i/n.groups*100),"%"))
		else setTxtProgressBar(pb, i)
	    }
    	    tmp.idx <- group.idx.start[i]:group.idx.end[i]
	    tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
    	    SeqArray::seqSetFilter(gds, variant.id = tmp.group.info$variant.idx, verbose = FALSE)
            geno <- if(is.dosage) SeqVarTools::imputedDosage(gds, use.names = FALSE) else SeqVarTools::altDosage(gds, use.names = FALSE)
            miss <- colMeans(is.na(geno))
            freq <- colMeans(geno, na.rm = TRUE)/2
	    include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
	    n.p <- sum(include)
	    if(n.p == 0) next
	    tmp.group.info <- tmp.group.info[include, , drop = FALSE]
	    miss <- miss[include]
	    freq <- freq[include]
            geno <- geno[, include, drop = FALSE]
	    N <- nrow(geno) - colSums(is.na(geno))
	    if(sum(tmp.group.info$flip) > 0) {
		freq[tmp.group.info$flip==1] <- 1 - freq[tmp.group.info$flip==1]
		geno[, tmp.group.info$flip==1] <- 2 - geno[, tmp.group.info$flip==1]
	    }
	    if(max(miss)>0) {
	    	miss.idx <- which(is.na(geno))
	    	geno[miss.idx] <- if(missing.method=="impute2mean") 2*freq[ceiling(miss.idx/nrow(geno))] else 0
	    }
	    U <- as.vector(crossprod(geno, residuals))
	    if(class(null.obj) == "glmmkin.multi") geno <- Diagonal(n = n.pheno) %x% geno
    	    if(!is.null(null.obj$P)) V <- crossprod(geno, crossprod(null.obj$P, geno))
	    else {
		GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
		V <- crossprod(geno, crossprod(null.obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
	    }
	    V <- as.matrix(V)
	    if(!is.null(meta.file.prefix)) {
		VAR <- diag(V)
		PVAL <- ifelse(VAR>0, pchisq(U^2/VAR, df=1, lower.tail=FALSE), NA)
		write.table(cbind(tmp.group.info[,c("group","chr","pos","ref","alt")], N, miss, freq, U, VAR, PVAL), meta.file.score, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
		writeBin(V[lower.tri(V, diag = TRUE)], meta.file.var.handle, size = 4)
	    }
	    if(use.minor.allele) {
		tmp.group.info$weight[freq > 0.5] <- -tmp.group.info$weight[freq > 0.5]
		freq[freq > 0.5] <- 1 - freq[freq > 0.5]
	    }
	    weights <- rep(tmp.group.info$weight * MAF.weights.beta.fun(freq, MAF.weights.beta[1], MAF.weights.beta[2]), n.pheno)
	    n.variants[i] <- n.p
	    miss.min[i] <- min(miss)
	    miss.mean[i] <- mean(miss)
	    miss.max[i] <- max(miss)
	    freq.min[i] <- min(freq)
	    freq.mean[i] <- mean(freq)
	    freq.max[i] <- max(freq)
	    U <- U*weights
	    V <- t(V*weights)*weights
	    if(max(V)-min(V) < sqrt(.Machine$double.eps)) {
	        burden.score <- sum(U)
    	        burden.var <- sum(V)
    	        burden.pval <- pchisq(burden.score^2/burden.var, df=1, lower.tail=FALSE)
	        if(Burden | SKATO | SMMAT) {
		    Burden.score[i] <- burden.score
	    	    Burden.var[i] <- burden.var
		    Burden.pval[i] <- burden.pval
	    	}
    	    	if(SKAT | SKATO) SKAT.pval[i] <- burden.pval
	    	if(SKATO) {
	            SKATO.pval[i] <- burden.pval
		    SKATO.minp[i] <- burden.pval
		    SKATO.minp.rho[i] <- 1
	    	}
    	    	if(SMMAT) SMMAT.pval[i] <- burden.pval
	    } else {
	        if(SKATO) {
		    re <- .skato_pval(U = U, V = V, rho = rho, method = method)
		    Burden.score[i] <- re$Burden.score
		    Burden.var[i] <- re$Burden.var
		    Burden.pval[i] <- re$Burden.pval
		    SKAT.pval[i] <- re$SKAT.pval
		    SKATO.pval[i] <-re$p
		    SKATO.minp[i] <- re$minp
		    SKATO.minp.rho[i] <- re$minp.rho
		} else {
		    if(SKAT) SKAT.pval[i] <- .quad_pval(U = U, V = V, method = method)
		    if(Burden | SMMAT) {
	    	        Burden.score[i] <- sum(U)
    	    	        Burden.var[i] <- sum(V)
    	    	        Burden.pval[i] <- pchisq(Burden.score[i]^2/Burden.var[i], df=1, lower.tail=FALSE)
		    }
		}	    
		if(SMMAT) {
	    	    V.rowSums <- rowSums(V)
    	    	    U <- U - V.rowSums * Burden.score[i] / Burden.var[i]
    	    	    V <- V - tcrossprod(V.rowSums) / Burden.var[i]
    	    	    if(mean(abs(V)) < sqrt(.Machine$double.eps)) SMMAT.pval[i] <- Burden.pval[i]
	    	    else SMMAT.pval[i] <- tryCatch(pchisq(-2*log(Burden.pval[i])-2*log(.quad_pval(U = U, V = V, method = method)), df = 4, lower.tail = FALSE), error = function(e) { Burden.pval[i] })
		}
	    }
	    rm(geno)
	    if(Garbage.Collection) gc()
    	}
    	SeqArray::seqClose(gds)
	if(verbose) {
	    if(is.Windows) setWinProgressBar(pb, n.groups, title="Progress of SMMAT: 100%")
	    else setTxtProgressBar(pb, n.groups)
	    close(pb)
	}
    	if(!is.null(meta.file.prefix)) close(meta.file.var.handle)
    	out <- data.frame(group=unique(group.info$group), n.variants=n.variants, miss.min=miss.min, miss.mean=miss.mean, miss.max=miss.max, freq.min=freq.min, freq.mean=freq.mean, freq.max=freq.max)
	if(Burden | SKATO | SMMAT) {
	    out$B.score <- Burden.score
	    out$B.var <- Burden.var
	    out$B.pval <- Burden.pval
	}
    	if(SKAT | SKATO) out$S.pval <- SKAT.pval
	if(SKATO) {
	    out$O.pval <- SKATO.pval
	    out$O.minp <- SKATO.minp
	    out$O.minp.rho <- SKATO.minp.rho
	}
    	if(SMMAT) out$E.pval <- SMMAT.pval
    }
    return(out[match(groups, out$group),])
}

SMMAT.prep <- function(null.obj, geno.file, group.file, group.file.sep = "\t", auto.flip = FALSE)
{
    if(!class(null.obj) %in% c("glmmkin", "glmmkin.multi")) stop("Error: null.obj must be a class glmmkin or glmmkin.multi object!")
    n.pheno <- null.obj$n.pheno
    if(any(duplicated(null.obj$id_include))) {
        J <- Matrix(sapply(unique(null.obj$id_include), function(x) 1*(null.obj$id_include==x)), sparse = TRUE)
        residuals <- as.numeric(as.matrix(crossprod(J, null.obj$scaled.residuals)))
        if(!is.null(null.obj$P)) null.obj$P <- as.matrix(crossprod(J, crossprod(null.obj$P, J)))
	else {
	    null.obj$Sigma_iX <- crossprod(J, null.obj$Sigma_iX)
	    null.obj$Sigma_i <- forceSymmetric(crossprod(J,crossprod(null.obj$Sigma_i,J)))
	    null.obj$Sigma_i <- Matrix(null.obj$Sigma_i, sparse = TRUE)
	}
        rm(J)
    } else residuals <- null.obj$scaled.residuals
    n <- length(unique(null.obj$id_include))
    if(!grepl("\\.gds$", geno.file)) stop("Error: currently only .gds format is supported in geno.file!")
    gds <- SeqArray::seqOpen(geno.file)
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    if(any(is.na(match(null.obj$id_include, sample.id)))) warning("Check your data... Some individuals in null.obj$id_include are missing in sample.id of geno.file!")
    sample.id <- sample.id[sample.id %in% null.obj$id_include]
    if(length(sample.id) == 0) stop("Error: null.obj$id_include does not match sample.id in geno.file!")
    match.id <- match(sample.id, unique(null.obj$id_include))
    if(class(null.obj) == "glmmkin.multi") {
        residuals <- residuals[match.id, , drop = FALSE]
    	match.id <- rep(match.id, n.pheno) + rep((0:(n.pheno-1))*n, each = length(match.id))
    } else {
        residuals <- residuals[match.id]
    }
    if(!is.null(null.obj$P)) null.obj$P <- null.obj$P[match.id, match.id]
    else {
    	null.obj$Sigma_iX <- null.obj$Sigma_iX[match.id, , drop = FALSE]
	null.obj$Sigma_i <- null.obj$Sigma_i[match.id, match.id]
    }
    variant.idx <- SeqArray::seqGetData(gds, "variant.id")
    #variant.id <- SeqArray::seqGetData(gds, "annotation/id")
    #variant.id <- paste(SeqArray::seqGetData(gds, "chromosome"), SeqArray::seqGetData(gds, "position"), SeqArray::seqGetData(gds, "allele"), sep=":")
    chr <- SeqArray::seqGetData(gds, "chromosome")
    pos <- SeqArray::seqGetData(gds, "position")
    alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
    ref <- unlist(lapply(alleles.list, function(x) x[1]))
    alt <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
    rm(alleles.list); gc()
    SeqArray::seqClose(gds)
    variant.id <- paste(chr, pos, ref, alt, sep = ":")
    rm(chr, pos, ref, alt); gc()
    group.info <- try(read.table(group.file, header = FALSE, col.names = c("group", "chr", "pos", "ref", "alt", "weight"), colClasses = c("character","character","integer","character","character","numeric"), sep = group.file.sep), silent = TRUE)
    if (class(group.info) == "try-error") {
        stop("Error: cannot read group.file!")
    }
    variant.id1 <- paste(group.info$chr, group.info$pos, group.info$ref, group.info$alt, sep = ":")
    group.info <- group.info[!duplicated(paste(group.info$group, variant.id1, sep = ":")), ]
    variant.idx1 <- variant.idx[match(variant.id1, variant.id)]
    group.info$variant.idx <- variant.idx1
    group.info$flip <- 0
    if(auto.flip) {
        cat("Automatic allele flipping enabled...\nVariants matching alt/ref but not ref/alt alleles will also be included, with flipped effects\n")
    	variant.id2 <- paste(group.info$chr, group.info$pos, group.info$alt, group.info$ref, sep = ":")
    	variant.idx2 <- variant.idx[match(variant.id2, variant.id)]
    	if(any(!is.na(variant.idx1) & !is.na(variant.idx2))) {
            tmp.dups <- which(!is.na(variant.idx1) & !is.na(variant.idx2))
	    cat("The following ambiguous variants were found:\n")
	    cat("chr:", chr[tmp.dups], "\n")
	    cat("pos:", pos[tmp.dups], "\n")
	    cat("ref:", ref[tmp.dups], "\n")
	    cat("alt:", alt[tmp.dups], "\n")
	    cat("Warning: both variants with alleles ref/alt and alt/ref were present at the same position and coding should be double checked!\nFor these variants, only those with alleles ref/alt were used in the analysis...\n")
	    variant.idx2[tmp.dups] <- NA
	    rm(tmp.dups)
    	}
    	group.info$flip <- 1*(!is.na(variant.idx2))
    	group.info$variant.idx[!is.na(variant.idx2)] <- variant.idx2[!is.na(variant.idx2)]
	rm(variant.id2, variant.idx2)
    }
    rm(variant.id, variant.id1, variant.idx1); gc()
    #group.info$variant.idx <- variant.idx[match(group.info$variant, variant.id)]
    group.info <- subset(group.info, !is.na(variant.idx))
    groups <- unique(group.info$group)
    n.groups.all <- length(groups)
    group.info$group.idx <- as.numeric(factor(group.info$group))
    group.info <- group.info[order(group.info$group.idx, group.info$variant.idx), ]
    group.idx.end <- findInterval(1:n.groups.all, group.info$group.idx)
    group.idx.start <- c(1, group.idx.end[-n.groups.all] + 1)
    out <- list(null.obj = null.obj, geno.file = geno.file, group.file = group.file, group.file.sep = group.file.sep, auto.flip = auto.flip, residuals = residuals, sample.id = sample.id, group.info = group.info, groups = groups, group.idx.start = group.idx.start, group.idx.end = group.idx.end)
    class(out) <- "SMMAT.prep"
    return(out)
}

SMMAT.lowmem <- function(SMMAT.prep.obj, meta.file.prefix = NULL, MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(1, 25), miss.cutoff = 1, missing.method = "impute2mean", method = "davies", tests = "E", rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1), use.minor.allele = FALSE, Garbage.Collection = FALSE, is.dosage = FALSE, ncores = 1, verbose = FALSE)
{
    if(class(SMMAT.prep.obj) != "SMMAT.prep") stop("Error: SMMAT.prep.obj must be a class SMMAT.prep object!")
    is.Windows <- Sys.info()["sysname"] == "Windows"
    if(is.Windows && ncores > 1) {
        warning("The package doMC is not available on Windows... Switching to single thread...")
        ncores <- 1
    }
    null.obj <- SMMAT.prep.obj$null.obj
    geno.file <- SMMAT.prep.obj$geno.file
    residuals <- SMMAT.prep.obj$residuals
    sample.id <- SMMAT.prep.obj$sample.id
    group.info <- SMMAT.prep.obj$group.info
    groups <- SMMAT.prep.obj$groups
    n.groups.all <- length(groups)
    group.idx.start <- SMMAT.prep.obj$group.idx.start
    group.idx.end <- SMMAT.prep.obj$group.idx.end
    rm(SMMAT.prep.obj); gc()
    if(!class(null.obj) %in% c("glmmkin", "glmmkin.multi")) stop("Error: null.obj must be a class glmmkin or glmmkin.multi object!")
    n.pheno <- null.obj$n.pheno
    missing.method <- try(match.arg(missing.method, c("impute2mean", "impute2zero")))
    if(class(missing.method) == "try-error") stop("Error: \"missing.method\" must be \"impute2mean\" or \"impute2zero\".")
    if(any(!tests %in% c("B", "S", "O", "E"))) stop("Error: \"tests\" should only include \"B\" for the burden test, \"S\" for SKAT, \"O\" for SKAT-O or \"E\" for the efficient hybrid test of the burden test and SKAT.")
    Burden <- "B" %in% tests
    SKAT <- "S" %in% tests
    SKATO <- "O" %in% tests
    SMMAT <- "E" %in% tests
    if(!grepl("\\.gds$", geno.file)) stop("Error: currently only .gds format is supported in geno.file!")
    ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
    if(ncores > 1) {
    	doMC::registerDoMC(cores = ncores)
    	n.groups.percore <- (n.groups.all-1) %/% ncores + 1
    	n.groups.percore_1 <- n.groups.percore * ncores - n.groups.all
	b <- NULL
    	out <- foreach(b=1:ncores, .combine=rbind, .multicombine = TRUE, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
    	    idx <- if(b <= n.groups.percore_1) ((b-1)*(n.groups.percore-1)+1):(b*(n.groups.percore-1)) else (n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1-1)*n.groups.percore+1):(n.groups.percore_1*(n.groups.percore-1)+(b-n.groups.percore_1)*n.groups.percore)
	    n.groups <- length(idx)
	    if(verbose) {
	        if(b==1) cat("Progress of SMMAT:\n")
		pb <- txtProgressBar(min = 0, max = n.groups, style = 3)
		cat("\n")
	    }
    	    gds <- SeqArray::seqOpen(geno.file)
    	    SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
	    n.variants <- rep(0,n.groups)
    	    miss.min <- rep(NA,n.groups)
    	    miss.mean <- rep(NA, n.groups)
    	    miss.max <- rep(NA, n.groups)
    	    freq.min <- rep(NA, n.groups)
    	    freq.mean <- rep(NA, n.groups)
    	    freq.max <- rep(NA, n.groups)
	    if(Burden | SKATO | SMMAT) {
	        Burden.score <- rep(NA, n.groups)
    	    	Burden.var <- rep(NA, n.groups)
    	    	Burden.pval <- rep(NA, n.groups)
	    }
    	    if(SKAT | SKATO) SKAT.pval <- rep(NA, n.groups)
	    if(SKATO) {
	        SKATO.pval <- rep(NA, n.groups)
		SKATO.minp <- rep(NA, n.groups)
		SKATO.minp.rho <- rep(NA, n.groups)
	    }
    	    if(SMMAT) SMMAT.pval <- rep(NA, n.groups)
	    if(!is.null(meta.file.prefix)) {
	    	if(class(null.obj) == "glmmkin.multi") stop("Error: meta-analysis not supported yet for multiple phenotypes.")
	        if(.Platform$endian!="little") stop("Error: platform must be little endian.")
		meta.file.score <- paste0(meta.file.prefix, ".score.", b)
		meta.file.var <- paste0(meta.file.prefix, ".var.", b)
		write.table(t(c("group", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE", "VAR", "PVAL")), meta.file.score, quote = FALSE, row.names = FALSE, col.names = FALSE)
		meta.file.var.handle <- file(meta.file.var, "wb")
	    }
    	    for(i in 1:n.groups) {
	        if(verbose && i %% ceiling(n.groups/100) == 0) setTxtProgressBar(pb, i)
    	    	tmp.idx <- group.idx.start[idx[i]]:group.idx.end[idx[i]]
		tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
	    	SeqArray::seqSetFilter(gds, variant.id = tmp.group.info$variant.idx, verbose = FALSE)
                geno <- if(is.dosage) SeqVarTools::imputedDosage(gds, use.names = FALSE) else SeqVarTools::altDosage(gds, use.names = FALSE)
                miss <- colMeans(is.na(geno))
                freq <- colMeans(geno, na.rm = TRUE)/2
	    	include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
	    	n.p <- sum(include)
	    	if(n.p == 0) next
		tmp.group.info <- tmp.group.info[include, , drop = FALSE]
	    	miss <- miss[include]
	    	freq <- freq[include]
                geno <- geno[, include, drop = FALSE]
		N <- nrow(geno) - colSums(is.na(geno))
		if(sum(tmp.group.info$flip) > 0) {
		    freq[tmp.group.info$flip==1] <- 1 - freq[tmp.group.info$flip==1]
		    geno[, tmp.group.info$flip==1] <- 2 - geno[, tmp.group.info$flip==1]
		}
	    	if(max(miss)>0) {
	    	    miss.idx <- which(is.na(geno))
	    	    geno[miss.idx] <- if(missing.method=="impute2mean") 2*freq[ceiling(miss.idx/nrow(geno))] else 0
	    	}
	    	U <- as.vector(crossprod(geno, residuals))
		if(class(null.obj) == "glmmkin.multi") geno <- Diagonal(n = n.pheno) %x% geno
		if(!is.null(null.obj$P)) V <- crossprod(geno, crossprod(null.obj$P, geno))
		else {
		    GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
		    V <- crossprod(geno, crossprod(null.obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
		}
	    	V <- as.matrix(V)
		if(!is.null(meta.file.prefix)) {
		    VAR <- diag(V)
		    PVAL <- ifelse(VAR>0, pchisq(U^2/VAR, df=1, lower.tail=FALSE), NA)
		    write.table(cbind(tmp.group.info[,c("group","chr","pos","ref","alt")], N, miss, freq, U, VAR, PVAL), meta.file.score, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
		    writeBin(V[lower.tri(V, diag = TRUE)], meta.file.var.handle, size = 4)
		}
		if(use.minor.allele) {
		    tmp.group.info$weight[freq > 0.5] <- -tmp.group.info$weight[freq > 0.5]
		    freq[freq > 0.5] <- 1 - freq[freq > 0.5]
		}
		weights <- rep(tmp.group.info$weight * MAF.weights.beta.fun(freq, MAF.weights.beta[1], MAF.weights.beta[2]), n.pheno)
	    	n.variants[i] <- n.p
	    	miss.min[i] <- min(miss)
	    	miss.mean[i] <- mean(miss)
	    	miss.max[i] <- max(miss)
	    	freq.min[i] <- min(freq)
	    	freq.mean[i] <- mean(freq)
	    	freq.max[i] <- max(freq)
		U <- U*weights
		V <- t(V*weights)*weights
		if(max(V)-min(V) < sqrt(.Machine$double.eps)) {
	    	    burden.score <- sum(U)
    	    	    burden.var <- sum(V)
    	    	    burden.pval <- pchisq(burden.score^2/burden.var, df=1, lower.tail=FALSE)
	    	    if(Burden | SKATO | SMMAT) {
		       	Burden.score[i] <- burden.score
	    		Burden.var[i] <- burden.var
			Burden.pval[i] <- burden.pval
	    	    }
    	    	    if(SKAT | SKATO) SKAT.pval[i] <- burden.pval
	    	    if(SKATO) {
	                SKATO.pval[i] <- burden.pval
			SKATO.minp[i] <- burden.pval
			SKATO.minp.rho[i] <- 1
	    	    }
    	    	    if(SMMAT) SMMAT.pval[i] <- burden.pval
		} else {
		    if(SKATO) {
		        re <- .skato_pval(U = U, V = V, rho = rho, method = method)
			Burden.score[i] <- re$Burden.score
			Burden.var[i] <- re$Burden.var
			Burden.pval[i] <- re$Burden.pval
			SKAT.pval[i] <- re$SKAT.pval
			SKATO.pval[i] <-re$p
			SKATO.minp[i] <- re$minp
			SKATO.minp.rho[i] <- re$minp.rho
		    } else {
		        if(SKAT) SKAT.pval[i] <- .quad_pval(U = U, V = V, method = method)
		    	if(Burden | SMMAT) {
	    	    	    Burden.score[i] <- sum(U)
    	    	    	    Burden.var[i] <- sum(V)
    	    	    	    Burden.pval[i] <- pchisq(Burden.score[i]^2/Burden.var[i], df=1, lower.tail=FALSE)
		    	}
		    }	    
		    if(SMMAT) {
	    	        V.rowSums <- rowSums(V)
    	    		U <- U - V.rowSums * Burden.score[i] / Burden.var[i]
    	    		V <- V - tcrossprod(V.rowSums) / Burden.var[i]
    	    		if(mean(abs(V)) < sqrt(.Machine$double.eps)) SMMAT.pval[i] <- Burden.pval[i]
	    		else SMMAT.pval[i] <- tryCatch(pchisq(-2*log(Burden.pval[i])-2*log(.quad_pval(U = U, V = V, method = method)), df = 4, lower.tail = FALSE), error = function(e) { Burden.pval[i] })
		    }
		}
	    	rm(geno)
	    	if(Garbage.Collection) gc()
    	    }
    	    SeqArray::seqClose(gds)
	    if(verbose) {
	        setTxtProgressBar(pb, n.groups)
		close(pb)
	    }
	    if(!is.null(meta.file.prefix)) close(meta.file.var.handle)
    	    tmp.out <- data.frame(group=unique(group.info$group)[idx], n.variants=n.variants, miss.min=miss.min, miss.mean=miss.mean, miss.max=miss.max, freq.min=freq.min, freq.mean=freq.mean, freq.max=freq.max)
	    if(Burden | SKATO | SMMAT) {
	        tmp.out$B.score <- Burden.score
		tmp.out$B.var <- Burden.var
		tmp.out$B.pval <- Burden.pval
	    }
    	    if(SKAT | SKATO) tmp.out$S.pval <- SKAT.pval
	    if(SKATO) {
	        tmp.out$O.pval <- SKATO.pval
		tmp.out$O.minp <- SKATO.minp
		tmp.out$O.minp.rho <- SKATO.minp.rho
	    }
    	    if(SMMAT) tmp.out$E.pval <- SMMAT.pval
	    tmp.out
	}
    } else { # use a single core
	n.groups <- n.groups.all
	if(verbose) {
	    if(is.Windows) pb <- winProgressBar(min = 0, max = n.groups)
	    else {
	        cat("Progress of SMMAT:\n")
	    	pb <- txtProgressBar(min = 0, max = n.groups, style = 3)
		cat("\n")
	    }
	}
    	gds <- SeqArray::seqOpen(geno.file)
    	SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
	n.variants <- rep(0,n.groups)
    	miss.min <- rep(NA,n.groups)
    	miss.mean <- rep(NA, n.groups)
    	miss.max <- rep(NA, n.groups)
    	freq.min <- rep(NA, n.groups)
    	freq.mean <- rep(NA, n.groups)
    	freq.max <- rep(NA, n.groups)
	if(Burden | SKATO | SMMAT) {
	    Burden.score <- rep(NA, n.groups)
    	    Burden.var <- rep(NA, n.groups)
    	    Burden.pval <- rep(NA, n.groups)
	}
    	if(SKAT | SKATO) SKAT.pval <- rep(NA, n.groups)
	if(SKATO) {
	    SKATO.pval <- rep(NA, n.groups)
	    SKATO.minp <- rep(NA, n.groups)
	    SKATO.minp.rho <- rep(NA, n.groups)
	}
    	if(SMMAT) SMMAT.pval <- rep(NA, n.groups)
    	if(!is.null(meta.file.prefix)) {
	    if(class(null.obj) == "glmmkin.multi") stop("Error: meta-analysis not supported yet for multiple phenotypes.")
            if(.Platform$endian!="little") stop("Error: platform must be little endian.")
	    meta.file.score <- paste0(meta.file.prefix, ".score.1")
	    meta.file.var <- paste0(meta.file.prefix, ".var.1")
	    write.table(t(c("group", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE", "VAR", "PVAL")), meta.file.score, quote = FALSE, row.names = FALSE, col.names = FALSE)
	    meta.file.var.handle <- file(meta.file.var, "wb")
	}
    	for(i in 1:n.groups) {
	    if(verbose && i %% ceiling(n.groups/100) == 0) {
		if(is.Windows) setWinProgressBar(pb, i, title=paste0("Progress of SMMAT: ",round(i/n.groups*100),"%"))
		else setTxtProgressBar(pb, i)
	    }
    	    tmp.idx <- group.idx.start[i]:group.idx.end[i]
	    tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
    	    SeqArray::seqSetFilter(gds, variant.id = tmp.group.info$variant.idx, verbose = FALSE)
            geno <- if(is.dosage) SeqVarTools::imputedDosage(gds, use.names = FALSE) else SeqVarTools::altDosage(gds, use.names = FALSE)
            miss <- colMeans(is.na(geno))
            freq <- colMeans(geno, na.rm = TRUE)/2
	    include <- (miss <= miss.cutoff & ((freq >= MAF.range[1] & freq <= MAF.range[2]) | (freq >= 1-MAF.range[2] & freq <= 1-MAF.range[1])))
	    n.p <- sum(include)
	    if(n.p == 0) next
	    tmp.group.info <- tmp.group.info[include, , drop = FALSE]
	    miss <- miss[include]
	    freq <- freq[include]
            geno <- geno[, include, drop = FALSE]
	    N <- nrow(geno) - colSums(is.na(geno))
	    if(sum(tmp.group.info$flip) > 0) {
		freq[tmp.group.info$flip==1] <- 1 - freq[tmp.group.info$flip==1]
		geno[, tmp.group.info$flip==1] <- 2 - geno[, tmp.group.info$flip==1]
	    }
	    if(max(miss)>0) {
	    	miss.idx <- which(is.na(geno))
	    	geno[miss.idx] <- if(missing.method=="impute2mean") 2*freq[ceiling(miss.idx/nrow(geno))] else 0
	    }
	    U <- as.vector(crossprod(geno, residuals))
	    if(class(null.obj) == "glmmkin.multi") geno <- Diagonal(n = n.pheno) %x% geno
    	    if(!is.null(null.obj$P)) V <- crossprod(geno, crossprod(null.obj$P, geno))
	    else {
		GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
		V <- crossprod(geno, crossprod(null.obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
	    }
	    V <- as.matrix(V)
	    if(!is.null(meta.file.prefix)) {
		VAR <- diag(V)
		PVAL <- ifelse(VAR>0, pchisq(U^2/VAR, df=1, lower.tail=FALSE), NA)
		write.table(cbind(tmp.group.info[,c("group","chr","pos","ref","alt")], N, miss, freq, U, VAR, PVAL), meta.file.score, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
		writeBin(V[lower.tri(V, diag = TRUE)], meta.file.var.handle, size = 4)
	    }
	    if(use.minor.allele) {
		tmp.group.info$weight[freq > 0.5] <- -tmp.group.info$weight[freq > 0.5]
		freq[freq > 0.5] <- 1 - freq[freq > 0.5]
	    }
	    weights <- rep(tmp.group.info$weight * MAF.weights.beta.fun(freq, MAF.weights.beta[1], MAF.weights.beta[2]), n.pheno)
	    n.variants[i] <- n.p
	    miss.min[i] <- min(miss)
	    miss.mean[i] <- mean(miss)
	    miss.max[i] <- max(miss)
	    freq.min[i] <- min(freq)
	    freq.mean[i] <- mean(freq)
	    freq.max[i] <- max(freq)
	    U <- U*weights
	    V <- t(V*weights)*weights
	    if(max(V)-min(V) < sqrt(.Machine$double.eps)) {
	        burden.score <- sum(U)
    	        burden.var <- sum(V)
    	        burden.pval <- pchisq(burden.score^2/burden.var, df=1, lower.tail=FALSE)
	        if(Burden | SKATO | SMMAT) {
		    Burden.score[i] <- burden.score
	    	    Burden.var[i] <- burden.var
		    Burden.pval[i] <- burden.pval
	    	}
    	    	if(SKAT | SKATO) SKAT.pval[i] <- burden.pval
	    	if(SKATO) {
	            SKATO.pval[i] <- burden.pval
		    SKATO.minp[i] <- burden.pval
		    SKATO.minp.rho[i] <- 1
	    	}
    	    	if(SMMAT) SMMAT.pval[i] <- burden.pval
	    } else {
	        if(SKATO) {
		    re <- .skato_pval(U = U, V = V, rho = rho, method = method)
		    Burden.score[i] <- re$Burden.score
		    Burden.var[i] <- re$Burden.var
		    Burden.pval[i] <- re$Burden.pval
		    SKAT.pval[i] <- re$SKAT.pval
		    SKATO.pval[i] <-re$p
		    SKATO.minp[i] <- re$minp
		    SKATO.minp.rho[i] <- re$minp.rho
		} else {
		    if(SKAT) SKAT.pval[i] <- .quad_pval(U = U, V = V, method = method)
		    if(Burden | SMMAT) {
	    	        Burden.score[i] <- sum(U)
    	    	        Burden.var[i] <- sum(V)
    	    	        Burden.pval[i] <- pchisq(Burden.score[i]^2/Burden.var[i], df=1, lower.tail=FALSE)
		    }
		}	    
		if(SMMAT) {
	    	    V.rowSums <- rowSums(V)
    	    	    U <- U - V.rowSums * Burden.score[i] / Burden.var[i]
    	    	    V <- V - tcrossprod(V.rowSums) / Burden.var[i]
    	    	    if(mean(abs(V)) < sqrt(.Machine$double.eps)) SMMAT.pval[i] <- Burden.pval[i]
	    	    else SMMAT.pval[i] <- tryCatch(pchisq(-2*log(Burden.pval[i])-2*log(.quad_pval(U = U, V = V, method = method)), df = 4, lower.tail = FALSE), error = function(e) { Burden.pval[i] })
		}
	    }
	    rm(geno)
	    if(Garbage.Collection) gc()
    	}
    	SeqArray::seqClose(gds)
	if(verbose) {
	    if(is.Windows) setWinProgressBar(pb, n.groups, title="Progress of SMMAT: 100%")
	    else setTxtProgressBar(pb, n.groups)
	    close(pb)
	}
    	if(!is.null(meta.file.prefix)) close(meta.file.var.handle)
    	out <- data.frame(group=unique(group.info$group), n.variants=n.variants, miss.min=miss.min, miss.mean=miss.mean, miss.max=miss.max, freq.min=freq.min, freq.mean=freq.mean, freq.max=freq.max)
	if(Burden | SKATO | SMMAT) {
	    out$B.score <- Burden.score
	    out$B.var <- Burden.var
	    out$B.pval <- Burden.pval
	}
    	if(SKAT | SKATO) out$S.pval <- SKAT.pval
	if(SKATO) {
	    out$O.pval <- SKATO.pval
	    out$O.minp <- SKATO.minp
	    out$O.minp.rho <- SKATO.minp.rho
	}
    	if(SMMAT) out$E.pval <- SMMAT.pval
    }
    return(out[match(groups, out$group),])
}

SMMAT.meta <- function(meta.files.prefix, n.files = rep(1, length(meta.files.prefix)), cohort.group.idx = NULL, group.file, group.file.sep = "\t", MAF.range = c(1e-7, 0.5), MAF.weights.beta = c(1, 25), miss.cutoff = 1, method = "davies", tests = "E", rho = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1), use.minor.allele = FALSE, verbose = FALSE)
{
    is.Windows <- Sys.info()["sysname"] == "Windows"
    if(.Platform$endian!="little") stop("Error: platform must be little endian.")
    n.cohort <- length(meta.files.prefix)
    if(length(n.files) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and n.files do not match.")
    if(!is.null(cohort.group.idx)) {
        if(length(cohort.group.idx) != n.cohort) stop("Error: numbers of cohorts specified in meta.files.prefix and cohort.group.idx do not match.")
	cohort.group.idx <- as.numeric(factor(cohort.group.idx))
	n.cohort.groups <- length(unique(cohort.group.idx))
    }
    if(any(!tests %in% c("B", "S", "O", "E"))) stop("Error: \"tests\" should only include \"B\" for the burden test, \"S\" for SKAT, \"O\" for SKAT-O or \"E\" for the efficient hybrid test of the burden test and SKAT.")
    Burden <- "B" %in% tests
    SKAT <- "S" %in% tests
    SKATO <- "O" %in% tests
    SMMAT <- "E" %in% tests
    group.info <- try(read.table(group.file, header = FALSE, col.names = c("group", "chr", "pos", "ref", "alt", "weight"), colClasses = c("character","character","integer","character","character","numeric"), sep = group.file.sep), silent = TRUE)
    if (class(group.info) == "try-error") {
        stop("Error: cannot read group.file!")
    }
    variant.id <- paste(group.info$group, group.info$chr, group.info$pos, group.info$ref, group.info$alt, sep = ":")
    group.info <- group.info[!duplicated(variant.id), ]
    variant.id <- variant.id[!duplicated(variant.id)]
    groups <- unique(group.info$group)
    n.groups <- length(groups)
    group.info$group.idx <- as.numeric(factor(group.info$group))
    sorted.order <- order(group.info$group.idx, group.info$chr, group.info$pos)
    group.info <- group.info[sorted.order, ]
    variant.id <- variant.id[sorted.order]
    rm(sorted.order)
    group.idx.end <- findInterval(1:n.groups, group.info$group.idx)
    group.idx.start <- c(1, group.idx.end[-n.groups] + 1)
    group.idx.ends <- group.idx.starts <- scores <- cons <- vector("list", n.cohort)
    if(verbose) {
    	if(is.Windows) pb <- winProgressBar(min = 0, max = n.cohort)
	else {
            cat("Progress in reading individual cohort results:\n")
	    pb <- txtProgressBar(min = 0, max = n.cohort, style = 3)
	}
    }
    for(i in 1:n.cohort) {
        tmp.scores <- NULL
	for(j in 1:n.files[i]) {
	    tmp <- try(read.table(paste0(meta.files.prefix[i], ".score.", j), header = TRUE, as.is = TRUE))
    	    if (class(tmp) == "try-error") {
                stop(paste0("Error: cannot read ", meta.files.prefix[i], ".score.", j, "!"))
    	    }
	    tmp <- tmp[,c("group", "chr", "pos", "ref", "alt", "N", "missrate", "altfreq", "SCORE")]
	    tmp$idx <- match(paste(tmp$group, tmp$chr, tmp$pos, tmp$ref, tmp$alt, sep = ":"), variant.id)
	    if(any(is.na(tmp$idx))) {
            	tmp.dups <- which(is.na(tmp$idx))
	    	cat("In", paste0(meta.files.prefix[i], ".score.", j), ", the following variants were not present in group.file:\n")
		cat("group:", tmp$group[tmp.dups], "\n")
	    	cat("chr:", tmp$chr[tmp.dups], "\n")
	    	cat("pos:", tmp$pos[tmp.dups], "\n")
	    	cat("ref:", tmp$ref[tmp.dups], "\n")
	    	cat("alt:", tmp$alt[tmp.dups], "\n")
	    	stop("Error: meta files possibly not generated using this group.file!")
    	    }
	    tmp$file <- j
	    tmp.scores <- rbind(tmp.scores, tmp)
	    rm(tmp)
	}
    	if(any(sort(tmp.scores$idx)!=tmp.scores$idx)) {
	    cat("In some", meta.files.prefix[i], "score files, the order of group and variants is not the same as in the group-sorted group.file.\n")
	    stop("Error: meta files possibly not generated using this group.file!")
    	}
    	group.idx.ends[[i]] <- findInterval(1:n.groups, group.info$group.idx[tmp.scores$idx])
    	group.idx.starts[[i]] <- c(1, group.idx.ends[[i]][-n.groups] + 1)
	scores[[i]] <- tmp.scores
	rm(tmp.scores)
	cons[[i]] <- file(paste0(meta.files.prefix[i], ".var.1"), "rb")
	if(verbose) {
	    if(is.Windows) setWinProgressBar(pb, i, title=paste0("Progress of reading: ",round(i/n.cohort*100),"%"))
	    else setTxtProgressBar(pb, i)
	}
    }
    if(verbose) close(pb)
    n.variants <- rep(0,n.groups)
    if(Burden | SKATO | SMMAT) {
	Burden.score <- rep(NA, n.groups)
    	Burden.var <- rep(NA, n.groups)
    	Burden.pval <- rep(NA, n.groups)
    }
    if(SKAT | SKATO) SKAT.pval <- rep(NA, n.groups)
    if(SKATO) {
	SKATO.pval <- rep(NA, n.groups)
	SKATO.minp <- rep(NA, n.groups)
	SKATO.minp.rho <- rep(NA, n.groups)
    }
    if(SMMAT) SMMAT.pval <- rep(NA, n.groups)
    current.lines <- current.cons <- rep(1, n.cohort)
    if(verbose) {
	if(is.Windows) pb <- winProgressBar(min = 0, max = n.groups)
	else {
	    cat("Progress of SMMAT meta-analysis:\n")
	    pb <- txtProgressBar(min = 0, max = n.groups, style = 3)
	}
    }
    for(i in 1:n.groups) {
	if(verbose && i %% ceiling(n.groups/100) == 0) {
	    if(is.Windows) setWinProgressBar(pb, i, title=paste0("Progress of meta-analysis: ",round(i/n.groups*100),"%"))
	    else setTxtProgressBar(pb, i)
	}
        tmp.idx <- group.idx.start[i]:group.idx.end[i]
	tmp.group.info <- group.info[tmp.idx, , drop = FALSE]
	U.list <- V.list <- vector("list", n.cohort)
	variant.indices <- tmp.N <- tmp.Nmiss <- tmp.AC <- c()
	for(j in 1:n.cohort) {
	    if(group.idx.starts[[j]][i] <= group.idx.ends[[j]][i]) {
	        tmp.n.p <- group.idx.ends[[j]][i]-group.idx.starts[[j]][i]+1
	        U.list[[j]] <- scores[[j]][current.lines[j]:(current.lines[j]+tmp.n.p-1), , drop = FALSE]
		if(all(U.list[[j]]$file==current.cons[j]+1)) {
		    close(cons[[j]])
		    current.cons[j] <- current.cons[j]+1
		    cons[[j]] <- file(paste0(meta.files.prefix[j], ".var.", current.cons[j]), "rb")
		} else if(any(U.list[[j]]$file!=current.cons[j])) {
	    	    cat("For test group", i, ", cohort", j, "has incorrect indices for binary var files.\n")
	    	    stop("Error: check your individual score files!")
		}
		tmp.V <- matrix(0, tmp.n.p, tmp.n.p)
		tmp.V[lower.tri(tmp.V, diag = TRUE)] <- readBin(cons[[j]], what = "numeric", n = (1+tmp.n.p)*tmp.n.p/2, size = 4)
		tmp.V[upper.tri(tmp.V)] <- t(tmp.V)[upper.tri(tmp.V)]
	    	V.list[[j]] <- tmp.V
		rm(tmp.V)
		current.lines[j] <- current.lines[j]+tmp.n.p
		variant.indices <- c(variant.indices, U.list[[j]]$idx)
		tmp.N <- c(tmp.N, U.list[[j]]$N)
		tmp.Nmiss <- c(tmp.Nmiss, U.list[[j]]$N * U.list[[j]]$missrate/(1-U.list[[j]]$missrate))
		tmp.AC <- c(tmp.AC, 2*U.list[[j]]$N*U.list[[j]]$altfreq)
	    }
	}
	if(length(variant.indices) == 0) next
	tmp.variant.indices <- variant.indices
	variant.indices <- sort(unique(variant.indices))
	N <- sapply(variant.indices, function(x) sum(tmp.N[tmp.variant.indices==x]))
	Nmiss <- sapply(variant.indices, function(x) sum(tmp.Nmiss[tmp.variant.indices==x]))
	AF <- sapply(variant.indices, function(x) sum(tmp.AC[tmp.variant.indices==x]))/2/N
	include <- (Nmiss/(N+Nmiss) <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
	rm(tmp.N, tmp.Nmiss, tmp.AC, tmp.variant.indices, N, Nmiss)
	if(sum(include) == 0) next
	variant.indices <- variant.indices[include]
	n.p <- length(variant.indices)
	n.variants[i] <- n.p
	U <- if(!is.null(cohort.group.idx)) rep(0, n.cohort.groups*n.p) else rep(0, n.p)
	V <- if(!is.null(cohort.group.idx)) matrix(0, n.cohort.groups*n.p, n.cohort.groups*n.p) else matrix(0, n.p, n.p)
	for(j in 1:n.cohort) {
	    if(!is.null(U.list[[j]]) & !is.null(V.list[[j]])) {
		IDX <- match(U.list[[j]]$idx, variant.indices)
		if(sum(!is.na(IDX)) == 0) next
		IDX2 <- which(!is.na(IDX))
		IDX <- IDX[IDX2]
		if(!is.null(cohort.group.idx)) IDX <- IDX+n.p*(cohort.group.idx[j]-1)
		U[IDX] <- U[IDX]+U.list[[j]]$SCORE[IDX2]
		V[IDX, IDX] <- V[IDX,IDX]+V.list[[j]][IDX2,IDX2]
	    }
	}
    	tmp.weight <- tmp.group.info$weight[match(variant.indices, tmp.idx)]
	if(use.minor.allele) tmp.weight[AF[include] > 0.5] <- -tmp.weight[AF[include] > 0.5]
	tmp.weight <- tmp.weight * MAF.weights.beta.fun(AF[include], MAF.weights.beta[1], MAF.weights.beta[2])
	if(!is.null(cohort.group.idx)) tmp.weight <- rep(tmp.weight, n.cohort.groups)
	U <- U*tmp.weight
	V <- t(V*tmp.weight)*tmp.weight
	if(max(V)-min(V) < sqrt(.Machine$double.eps)) {
	    burden.score <- sum(U)
    	    burden.var <- sum(V)
    	    burden.pval <- pchisq(burden.score^2/burden.var, df=1, lower.tail=FALSE)
	    if(Burden | SKATO | SMMAT) {
		Burden.score[i] <- burden.score
	    	Burden.var[i] <- burden.var
		Burden.pval[i] <- burden.pval
	    }
    	    if(SKAT | SKATO) SKAT.pval[i] <- burden.pval
	    if(SKATO) {
	        SKATO.pval[i] <- burden.pval
		SKATO.minp[i] <- burden.pval
		SKATO.minp.rho[i] <- 1
	    }
    	    if(SMMAT) SMMAT.pval[i] <- burden.pval
	} else {
	    if(SKATO) {
		re <- .skato_pval(U = U, V = V, rho = rho, method = method)
		Burden.score[i] <- re$Burden.score
		Burden.var[i] <- re$Burden.var
		Burden.pval[i] <- re$Burden.pval
		SKAT.pval[i] <- re$SKAT.pval
		SKATO.pval[i] <-re$p
		SKATO.minp[i] <- re$minp
		SKATO.minp.rho[i] <- re$minp.rho
	    } else {
		if(SKAT) SKAT.pval[i] <- .quad_pval(U = U, V = V, method = method)
		if(Burden | SMMAT) {
	    	    Burden.score[i] <- sum(U)
    	    	    Burden.var[i] <- sum(V)
    	    	    Burden.pval[i] <- pchisq(Burden.score[i]^2/Burden.var[i], df=1, lower.tail=FALSE)
		}
	    }	    
	    if(SMMAT) {
	    	V.rowSums <- rowSums(V)
    	    	U <- U - V.rowSums * Burden.score[i] / Burden.var[i]
    	    	V <- V - tcrossprod(V.rowSums) / Burden.var[i]
    	    	if(mean(abs(V)) < sqrt(.Machine$double.eps)) SMMAT.pval[i] <- Burden.pval[i]
	    	else SMMAT.pval[i] <- tryCatch(pchisq(-2*log(Burden.pval[i])-2*log(.quad_pval(U = U, V = V, method = method)), df = 4, lower.tail = FALSE), error = function(e) { Burden.pval[i] })
	    }
	}
    }
    if(verbose) {
        if(is.Windows) setWinProgressBar(pb, n.groups, title="Progress of meta-analysis: 100%")
    	else setTxtProgressBar(pb, n.groups)
	close(pb)
    }
    for(i in 1:n.cohort) close(cons[[i]])
    out <- data.frame(group=unique(group.info$group), n.variants=n.variants)
    if(Burden | SKATO | SMMAT) {
	out$B.score <- Burden.score
	out$B.var <- Burden.var
	out$B.pval <- Burden.pval
    }
    if(SKAT | SKATO) out$S.pval <- SKAT.pval
    if(SKATO) {
	out$O.pval <- SKATO.pval
	out$O.minp <- SKATO.minp
	out$O.minp.rho <- SKATO.minp.rho
    }
    if(SMMAT) out$E.pval <- SMMAT.pval
    return(out[match(groups, out$group),])
}

.skato_pval <- function(U, V, rho, method = "davies") {
    n.r <- length(rho)
    n.p <- length(U)
    lambdas <- vector("list", n.r)
    pval <- qval <- rep(NA, n.r)
    Q <- (1-rho)*sum(U^2)+rho*sum(U)^2
    Burden.score <- Burden.var <- Burden.pval <- SKAT.pval <- NA
    for(i in 1:n.r) {
	if(rho[i]==1) {
	    Burden.score <- sum(U)
	    Burden.var <- sum(V)
	    Burden.pval <- pchisq(Burden.score^2/Burden.var, df=1, lower.tail=FALSE)
	    lambdas[[i]] <- Burden.var
	    pval[i] <- Burden.pval
	    next
	}
	if(rho[i]!=0) {
	    R.M <- matrix(rho[i], n.p, n.p)
	    diag(R.M) <- 1
	    R.M.chol <- t(chol(R.M, pivot = TRUE))
	    V.temp <- crossprod(R.M.chol, crossprod(V, R.M.chol))
	} else V.temp <- V
	lambda <- eigen(V.temp, only.values = TRUE, symmetric=TRUE)$values
    	lambdas[[i]] <- lambda[lambda > 0]
	pval[i] <- .Q_pval(Q[i], lambdas[[i]], method = method)
	if(rho[i]==0) SKAT.pval <- pval[i]
    }
    minp <- min(pval)
    for(i in 1:n.r) {
	df <- sum(lambdas[[i]]^2)^2/sum(lambdas[[i]]^4)
	qval[i] <- (qchisq(minp, df, lower.tail = FALSE)-df)/sqrt(2*df)*sqrt(2*sum(lambdas[[i]]^2))+sum(lambdas[[i]])
    }
    ZMZ <- tcrossprod(rowSums(V))/sum(V)
    V.temp <- V - ZMZ
    lambda <- eigen(V.temp, only.values = TRUE, symmetric = TRUE)$values
    lambda <- lambda[lambda > 0]
    muq <- sum(lambda)
    varq <- sum(lambda^2) * 2 + sum(ZMZ * V.temp) * 4
    df <- sum(lambda^2)^2/sum(lambda^4)
    tau <- rho * sum(V) + sum(V %*% V)/sum(V) * (1 - rho)
    re <- tryCatch({
        integrate(function(x){
    	    t1 <- tau %x% t(x)
    	    re<-pchisq((apply((qval - t1)/(1-rho),2,min) - muq)/sqrt(varq)*sqrt(2*df) + df, df=df) * dchisq(x,df=1)
    	    return(re)
	}, lower = 0, upper = 40, subdivisions = 2000, abs.tol = 10^-25)
    }, error=function(e) NA)
    return(list(p = min(1-re[[1]], minp*n.r), minp = minp, minp.rho = rho[which.min(pval)],
    Burden.score=Burden.score, Burden.var=Burden.var, Burden.pval=Burden.pval,
    SKAT.pval=SKAT.pval))
}

.Q_pval <- function(Q, lambda, method = "davies") {
    if(method == "davies") {
        tmp <- suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-6))
        pval <- tmp$Qq
	if((tmp$ifault > 0) | (pval <= 1e-5) | (pval >= 1)) method <- "kuonen"
    }
    if(method == "kuonen") {
    	pval <- .pKuonen(x = Q, lambda = lambda)
	if(is.na(pval)) method <- "liu"
    }
    if(method == "liu") pval <- CompQuadForm::liu(q = Q, lambda = lambda)
    return(pval)
}

.quad_pval <- function(U, V, method = "davies") {
    Q <- sum(U^2)
    lambda <- eigen(V, only.values = TRUE, symmetric=TRUE)$values
    lambda <- lambda[lambda > 0]
    pval <- .Q_pval(Q, lambda, method = method)
    return(pval)
}

.pKuonen<-function (x, lambda, delta = rep(0, length(lambda)), df = rep(1, length(lambda)))
{
    delta <- delta[lambda != 0]
    df <- df[lambda != 0]
    lambda <- lambda[lambda != 0]
    if(length(lambda) != length(delta)) stop("Error: inconsistent length in lambda and delta!")
    if(length(lambda) != length(df)) stop("Error: inconsistent length in lambda and df!")
    if (length(lambda) == 1) {
        pchisq(x/lambda, df = df, ncp = delta, lower.tail = FALSE)
    }
    d <- max(lambda)
    lambda <- lambda/d
    x <- x/d
    k0 <- function(zeta) {
        -sum(df * log(1 - 2 * zeta * lambda))/2 + sum((delta * lambda *
            zeta)/(1 - 2 * zeta * lambda))
    }
    kprime0 <- function(zeta) {
        sapply(zeta, function(zz) {
            sum(((delta + df) * lambda)/(1 - 2 * zz * lambda) + 2 * (delta *
                zz * lambda^2)/(1 - 2 * zz * lambda)^2)
        })
    }
    kpprime0 <- function(zeta) {
        sum((2 * (2 * delta + df) * lambda^2)/(1 - 2 * zeta * lambda)^2 + 8 *
            delta * zeta * lambda^3/(1 - 2 * zeta * lambda)^3)
    }
    if (any(lambda < 0)) {
        lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
    }
    else if (x > sum((df+delta)*lambda)) {
        lmin <- -0.01
    }
    else {
        lmin <- -length(lambda)*max(df+delta)/(2 * x)
    }
    lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
    hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin,
        upper = lmax, tol = 1e-08)$root
    w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v <- hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-04)
        NA
    else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}

MAF.weights.beta.fun <- function(freq, beta1, beta2) {
    freq[freq > 0.5] <- 1 - freq[freq > 0.5]
    ifelse(freq <= 0, 0, dbeta(freq, beta1, beta2))
}
