glmm.score <- function(obj, infile, outfile, BGEN.samplefile = NULL, center = T, select = NULL, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, missing.method = "impute2mean", nperbatch = 100, tol = 1e-5, infile.nrow = NULL, infile.nrow.skip = 0, infile.sep = "\t", infile.na = "NA", infile.ncol.skip = 1, infile.ncol.print = 1, infile.header.print = "SNP", is.dosage = FALSE, ncores = 1, verbose = FALSE) {
  is.Windows <- Sys.info()["sysname"] == "Windows"
  if(!class(obj) %in% c("glmmkin", "glmmkin.multi")) stop("Error: obj must be a class glmmkin or glmmkin.multi object!")
  n.pheno <- obj$n.pheno
  if(any(duplicated(obj$id_include))) {
    J <- Matrix(sapply(unique(obj$id_include), function(x) 1*(obj$id_include==x)), sparse = TRUE)
    res <- as.vector(as.matrix(crossprod(J, obj$scaled.residuals)))
    if(!is.null(obj$P)) obj$P <- as.matrix(crossprod(J, crossprod(obj$P, J)))
    else {
      obj$Sigma_iX <- crossprod(J, obj$Sigma_iX)
      obj$Sigma_i <- forceSymmetric(crossprod(J,crossprod(obj$Sigma_i,J)))
      obj$Sigma_i <- Matrix(obj$Sigma_i, sparse = TRUE)
    }
    rm(J)
  } else res <- obj$scaled.residuals
  n <- length(unique(obj$id_include))
  miss.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
  if(class(miss.method) == "try-error") stop("Error: missing.method should be one of the following: impute2mean, omit!")
  miss.method <- substr(miss.method, 1, 1)
  if(all(file.exists(paste(infile, c("bim", "bed", "fam"), sep=".")))) {
    if(ncores != 1) stop("Error: parallel computing currently not implemented for PLINK binary format genotypes.")
    if(class(obj) == "glmmkin.multi") stop("Error: multiple phenotypes currently not implemented for PLINK binary format genotypes.")
    bimfile <- paste(infile, "bim", sep=".")
    bedfile <- paste(infile, "bed", sep=".")
    famfile <- paste(infile, "fam", sep=".")
    sample.id <- read.table(famfile, as.is=T)[,2]
    if(is.null(select)) {
      if(any(is.na(match(unique(obj$id_include), sample.id)))) warning("Check your data... Some id_include in obj are missing in sample.id of infile!")
      select <- match(sample.id, unique(obj$id_include))
      select[is.na(select)] <- 0
      if(all(select == 0)) stop("Error: id_include in obj does not match sample.id in infile!")
    }
    if(length(select) != length(sample.id)) stop("Error: number of individuals in select does not match infile!")
    select2 <- select[select > 0]
    if(any(duplicated(select2)) || max(select2) > length(unique(obj$id_include))) stop("Error: select is a vector of orders, individuals not in obj should be coded 0!")
    res <- res[select2]
    if(!is.null(obj$P)) obj$P <- obj$P[select2, select2]
    else {
      obj$Sigma_iX <- Matrix(obj$Sigma_iX[select2, , drop = FALSE], sparse = TRUE)
      obj$Sigma_i <- Matrix(obj$Sigma_i[select2, select2], sparse = TRUE)
      obj$cov <- Matrix(obj$cov, sparse = TRUE)
    }
    select[select > 0] <- 1:sum(select > 0)
    rm(select2)
    if(center) {
      if(!is.null(obj$P)) time <- .Call(C_glmm_score_bed, res, obj$P, bimfile, bedfile, outfile, 'c', MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select)
      else time <- .Call(C_glmm_score_bed_sp, res, obj$Sigma_i, obj$Sigma_iX, obj$cov, bimfile, bedfile, outfile, 'c', MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select)
    } else {
      if(!is.null(obj$P)) time <- .Call(C_glmm_score_bed, res, obj$P, bimfile, bedfile, outfile, 'n', MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select)
      else time <- .Call(C_glmm_score_bed_sp, res, obj$Sigma_i, obj$Sigma_iX, obj$cov, bimfile, bedfile, outfile, 'n', MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select)
    }
    #print(sprintf("Computational time: %.2f seconds", time))
    return(invisible(time))
  } else if(grepl("\\.gds$", infile)) { # GDS genotype file
    if(is.Windows && ncores > 1) {
      warning("The package doMC is not available on Windows... Switching to single thread...")
      ncores <- 1
    }
    ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
    gds <- SeqArray::seqOpen(infile)
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    if(is.null(select)) {
      if(any(is.na(match(unique(obj$id_include), sample.id)))) warning("Check your data... Some id_include in obj are missing in sample.id of infile!")
      select <- match(sample.id, unique(obj$id_include))
      select[is.na(select)] <- 0
      if(all(select == 0)) stop("Error: id_include in obj does not match sample.id in infile!")
    }
    if(length(select) != length(sample.id)) stop("Error: number of individuals in select does not match infile!")
    select2 <- select[select > 0]
    if(any(duplicated(select2)) || max(select2) > length(unique(obj$id_include))) stop("Error: select is a vector of orders, individuals not in obj should be coded 0!")
    if(class(obj) == "glmmkin.multi") {
      res <- res[select2, , drop = FALSE]
      select2 <- rep(select2, n.pheno) + rep((0:(n.pheno-1))*n, each = length(select2))
    } else {
      res <- res[select2]
    }
    if(!is.null(obj$P)) obj$P <- obj$P[select2, select2]
    else {
      obj$Sigma_iX <- obj$Sigma_iX[select2, , drop = FALSE]
      obj$Sigma_i <- obj$Sigma_i[select2, select2]
    }
    rm(select2)
    variant.idx.all <- SeqArray::seqGetData(gds, "variant.id")
    SeqArray::seqClose(gds)
    p.all <- length(variant.idx.all)
    if(ncores > 1) {
      doMC::registerDoMC(cores = ncores)
      p.percore <- (p.all-1) %/% ncores + 1
      n.p.percore_1 <- p.percore * ncores - p.all
      foreach(b = 1:ncores, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
        variant.idx <- if(b <= n.p.percore_1) variant.idx.all[((b-1)*(p.percore-1)+1):(b*(p.percore-1))] else variant.idx.all[(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1-1)*p.percore+1):(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1)*p.percore)]
        p <- length(variant.idx)
	if(verbose) {
	  if(b==1) cat("Progress of score test:\n")
	  pb <- txtProgressBar(min = 0, max = p, style = 3)
	  cat("\n")
	}
        gds <- SeqArray::seqOpen(infile)
        SeqArray::seqSetFilter(gds, sample.id = sample.id[select > 0], verbose = FALSE)
        rm(sample.id); rm(select)
        nbatch.flush <- (p-1) %/% 100000 + 1
        ii <- 0
        for(i in 1:nbatch.flush) {
	  if(verbose) setTxtProgressBar(pb, (i-1)*100000)
          gc()
          tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*100000+1):p] else variant.idx[((i-1)*100000+1):(i*100000)]
          SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
          MISSRATE <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
          AF <- 1 - SeqVarTools::alleleFrequency(gds)
          include <- (MISSRATE <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
          if(sum(include) == 0) next
          ii <- ii + 1
          tmp.variant.idx <- tmp.variant.idx[include]
          tmp.p <- length(tmp.variant.idx)
          SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
          SNP <- SeqArray::seqGetData(gds, "annotation/id")
          SNP[SNP == ""] <- NA
          out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
          rm(SNP)
          alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
          out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
          out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
          out$MISSRATE <- MISSRATE[include]
          out$AF <- AF[include]
          rm(alleles.list, include)
          tmp.out <- lapply(1:((tmp.p-1) %/% nperbatch + 1), function(j) {
            tmp2.variant.idx <- if(j == (tmp.p-1) %/% nperbatch + 1) tmp.variant.idx[((j-1)*nperbatch+1):tmp.p] else tmp.variant.idx[((j-1)*nperbatch+1):(j*nperbatch)]
            SeqArray::seqSetFilter(gds, variant.id = tmp2.variant.idx, verbose = FALSE)
            geno <- if(is.dosage) SeqVarTools::imputedDosage(gds, use.names = FALSE) else SeqVarTools::altDosage(gds, use.names = FALSE)
            N <- nrow(geno) - colSums(is.na(geno))
            if(center) geno <- scale(geno, scale = FALSE)
            miss.idx <- which(is.na(geno))
            if(length(miss.idx)>0) {
              geno[miss.idx] <- if(!center & missing.method == "impute2mean") colMeans(geno, na.rm = TRUE)[ceiling(miss.idx/nrow(geno))] else 0
            }
            if(class(obj) == "glmmkin.multi") {
              SCORE <- crossprod(res, geno)
              geno <- Diagonal(n = n.pheno) %x% geno
            } else {
              SCORE <- as.vector(crossprod(geno, res))
            }
            if(!is.null(obj$P)) VAR <- crossprod(geno, crossprod(obj$P, geno))
            else {
              GSigma_iX <- crossprod(geno, obj$Sigma_iX)
              VAR <- crossprod(geno, crossprod(obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, obj$cov))
            }
            if(class(obj) == "glmmkin.multi") {
              VAR_PVAL <- sapply(1:length(N), function(xx) {
                VAR2 <- as.matrix(VAR[(0:(n.pheno-1))*length(N)+xx,(0:(n.pheno-1))*length(N)+xx])
                PVAL <- try(pchisq(crossprod(SCORE[,xx], solve(VAR2, SCORE[,xx])), df = n.pheno, lower.tail = FALSE))
                if(class(PVAL)[1] == "try-error") PVAL <- NA
                c(VAR2[lower.tri(VAR2, diag = TRUE)], PVAL)
              })
              return(rbind(N, SCORE, VAR_PVAL))
            } else {
              VAR <- diag(VAR)
              PVAL <- ifelse(VAR>0, pchisq(SCORE^2/VAR, df=1, lower.tail=FALSE), NA)
              return(rbind(N, SCORE, VAR, PVAL))
            }
          })
          if(class(obj) == "glmmkin.multi") {
            tmp.out <- matrix(unlist(tmp.out), ncol = 2+n.pheno+n.pheno*(n.pheno+1)/2, byrow = TRUE, dimnames = list(NULL, c("N", paste0("SCORE",1:n.pheno), paste0("VAR",rep(1:n.pheno,n.pheno:1),unlist(lapply(1:n.pheno,function(xx) xx:n.pheno))), "PVAL")))
            out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N"], out[,c("MISSRATE","AF")], tmp.out[,c(paste0("SCORE",1:n.pheno), paste0("VAR",rep(1:n.pheno,n.pheno:1),unlist(lapply(1:n.pheno,function(xx) xx:n.pheno))), "PVAL")])
          } else {
            tmp.out <- matrix(unlist(tmp.out), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("N", "SCORE", "VAR", "PVAL")))
            out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N"], out[,c("MISSRATE","AF")], tmp.out[,c("SCORE","VAR","PVAL")])
          }
          names(out)[6] <- "N"
          rm(tmp.out)
          if(b == 1) {
            write.table(out, outfile, quote=FALSE, row.names=FALSE, col.names=(ii == 1), sep="\t", append=(ii > 1), na="NA")
          } else {
            write.table(out, paste0(outfile, "_tmp.", b), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=(ii > 1), na="NA")
          }
          rm(out)
        }
        SeqArray::seqClose(gds)
	if(verbose) {
	  setTxtProgressBar(pb, p)
	  close(pb)
	}
      }
      for(b in 2:ncores) {
        system(paste0("cat ", outfile, "_tmp.", b, " >> ", outfile))
        unlink(paste0(outfile, "_tmp.", b))
      }
    } else {
      variant.idx <- variant.idx.all
      rm(variant.idx.all)
      p <- length(variant.idx)
      if(verbose) {
	if(is.Windows) pb <- winProgressBar(min = 0, max = p)
	else {
          cat("Progress of score test:\n")
	  pb <- txtProgressBar(min = 0, max = p, style = 3)
	  cat("\n")
	}
      }
      gds <- SeqArray::seqOpen(infile)
      SeqArray::seqSetFilter(gds, sample.id = sample.id[select > 0], verbose = FALSE)
      rm(sample.id); rm(select)
      nbatch.flush <- (p-1) %/% 100000 + 1
      ii <- 0
      for(i in 1:nbatch.flush) {
	if(verbose) {
	  if(is.Windows) setWinProgressBar(pb, (i-1)*100000, title=paste0("Progress of score test: ",round((i-1)*100000/p*100),"%"))
	  else setTxtProgressBar(pb, (i-1)*100000)
	}
        gc()
        tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*100000+1):p] else variant.idx[((i-1)*100000+1):(i*100000)]
        tmp.p <- length(tmp.variant.idx)
        SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
        MISSRATE <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
        AF <- 1 - SeqVarTools::alleleFrequency(gds)
        include <- (MISSRATE <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
        if(sum(include) == 0) next
        ii <- ii + 1
        tmp.variant.idx <- tmp.variant.idx[include]
        tmp.p <- length(tmp.variant.idx)
        SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
        SNP <- SeqArray::seqGetData(gds, "annotation/id")
        SNP[SNP == ""] <- NA
        out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
        rm(SNP)
        alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
        out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
        out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
        out$MISSRATE <- MISSRATE[include]
        out$AF <- AF[include]
        rm(alleles.list, include)
        tmp.out <- lapply(1:((tmp.p-1) %/% nperbatch + 1), function(j) {
          tmp2.variant.idx <- if(j == (tmp.p-1) %/% nperbatch + 1) tmp.variant.idx[((j-1)*nperbatch+1):tmp.p] else tmp.variant.idx[((j-1)*nperbatch+1):(j*nperbatch)]
          SeqArray::seqSetFilter(gds, variant.id = tmp2.variant.idx, verbose = FALSE)
          geno <- if(is.dosage) SeqVarTools::imputedDosage(gds, use.names = FALSE) else SeqVarTools::altDosage(gds, use.names = FALSE)
          N <- nrow(geno) - colSums(is.na(geno))
          if(center) geno <- scale(geno, scale = FALSE)
          miss.idx <- which(is.na(geno))
          if(length(miss.idx)>0) {
            geno[miss.idx] <- if(!center & missing.method == "impute2mean") colMeans(geno, na.rm = TRUE)[ceiling(miss.idx/nrow(geno))] else 0
          }
          if(class(obj) == "glmmkin.multi") {
            SCORE <- crossprod(res, geno)
            geno <- Diagonal(n = n.pheno) %x% geno
          } else {
            SCORE <- as.vector(crossprod(geno, res))
          }
          if(!is.null(obj$P)) VAR <- crossprod(geno, crossprod(obj$P, geno))
          else {
            GSigma_iX <- crossprod(geno, obj$Sigma_iX)
            VAR <- crossprod(geno, crossprod(obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, obj$cov))
          }
          if(class(obj) == "glmmkin.multi") {
            VAR_PVAL <- sapply(1:length(N), function(xx) {
              VAR2 <- as.matrix(VAR[(0:(n.pheno-1))*length(N)+xx,(0:(n.pheno-1))*length(N)+xx])
              PVAL <- try(pchisq(crossprod(SCORE[,xx], solve(VAR2, SCORE[,xx])), df = n.pheno, lower.tail = FALSE))
              if(class(PVAL)[1] == "try-error") PVAL <- NA
              c(VAR2[lower.tri(VAR2, diag = TRUE)], PVAL)
            })
            return(rbind(N, SCORE, VAR_PVAL))
          } else {
            VAR <- diag(VAR)
            PVAL <- ifelse(VAR>0, pchisq(SCORE^2/VAR, df=1, lower.tail=FALSE), NA)
            return(rbind(N, SCORE, VAR, PVAL))
          }
        })
        if(class(obj) == "glmmkin.multi") {
          tmp.out <- matrix(unlist(tmp.out), ncol = 2+n.pheno+n.pheno*(n.pheno+1)/2, byrow = TRUE, dimnames = list(NULL, c("N", paste0("SCORE",1:n.pheno), paste0("VAR",rep(1:n.pheno,n.pheno:1),unlist(lapply(1:n.pheno,function(xx) xx:n.pheno))), "PVAL")))
          out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N"], out[,c("MISSRATE","AF")], tmp.out[,c(paste0("SCORE",1:n.pheno), paste0("VAR",rep(1:n.pheno,n.pheno:1),unlist(lapply(1:n.pheno,function(xx) xx:n.pheno))), "PVAL")])
        } else {
          tmp.out <- matrix(unlist(tmp.out), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("N", "SCORE", "VAR", "PVAL")))
          out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N"], out[,c("MISSRATE","AF")], tmp.out[,c("SCORE","VAR","PVAL")])
        }
        names(out)[6] <- "N"
        rm(tmp.out)
        write.table(out, outfile, quote=FALSE, row.names=FALSE, col.names=(ii == 1), sep="\t", append=(ii > 1), na="NA")
        rm(out)
      }
      SeqArray::seqClose(gds)
      if(verbose) {
	if(is.Windows) setWinProgressBar(pb, p, title="Progress of score test: 100%")
	else setTxtProgressBar(pb, p)
	close(pb)
      }
    }
    return(invisible(NULL))
  } else if (grepl("\\.bgen$", infile)){
    if(class(obj) == "glmmkin.multi") stop("Error: multiple phenotypes currently not implemented for BGEN format genotypes.")
    if(is.Windows && ncores > 1) {
      warning("The package doMC is not available on Windows... Switching to single thread...")
      ncores <- 1
    }
    ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
    bgenInfo <- .Call(C_bgenHeader, infile)
    if(!is.null(select)){
      if(length(select) != bgenInfo$N) stop("Error: number of individuals in select does not match infile!")
    } else if (!is.null(BGEN.samplefile)) {
      sampleFile <- read.table(file = BGEN.samplefile, header = TRUE, sep = " ")
      if ((nrow(sampleFile)-1) != bgenInfo$N){
        stop(paste0("Error: Number of sample identifiers in sample file (", nrow(sampleFile)-1, ") does not match number of samples in BGEN file (", bgenInfo$N,")."))
      }
      sample.id <- sampleFile[-1,1]
      if(any(is.na(match(unique(obj$id_include), sample.id)))) warning("Check your data... Some id_include in obj are missing in BGEN.samplefile!")
      select <- match(sample.id, unique(obj$id_include))
      select[is.na(select)] <- 0
      if(all(select == 0)) stop("Error: id_include in obj does not match sample.id in BGEN.samplefile!")
      
    } else {
      if (is.null(BGEN.samplefile) && bgenInfo$SampleIdFlag == 0) {
        stop("Error: BGEN file does not contain sample identifiers. Use the BGEN.samplefile or select parameter.")
      } 
      if(any(is.na(match(unique(obj$id_include), bgenInfo$SampleIds)))) warning("Check your data... Some id_include in obj are missing in sample.id of infile!")
      select <- match(bgenInfo$SampleIds, unique(obj$id_include))
      select[is.na(select)] <- 0
      if(all(select == 0)) stop("Error: id_include in obj does not match sample.id in infile!")
    }
    select2 <- select[select > 0]
    if(any(duplicated(select2)) || max(select2) > length(unique(obj$id_include))) stop("Error: select is a vector of orders, individuals not in obj should be coded 0!")
   
    res <- res[select2]
    if(!is.null(obj$P)) obj$P <- obj$P[select2, select2]
    else {
      obj$Sigma_iX <- Matrix(obj$Sigma_iX[select2, , drop = FALSE], sparse = TRUE)
      obj$Sigma_i <- Matrix(obj$Sigma_i[select2, select2], sparse = TRUE)
      obj$cov <- Matrix(obj$cov, sparse = TRUE)
    }
    select[select > 0] <- 1:sum(select > 0)
    rm(select2)
    center2 <- ifelse(center, 'c', 'n')
    
    if (ncores > 1){
      if (ncores > bgenInfo$M) {
        ncores = bgenInfo$M
        warning(paste0("The number of threads specified is greater than the number of variants. Using ", ncores, " cores instead."))
      }
      
      threadInfo <- .Call(C_getVariantPos, infile, bgenInfo$offset, bgenInfo$M, bgenInfo$N, bgenInfo$CompressionFlag, bgenInfo$LayoutFlag, ncores)

      doMC::registerDoMC(cores = ncores)
      
      foreach(i=1:ncores) %dopar% {
        if (!is.null(obj$P)) {
          if (bgenInfo$LayoutFlag == 2) {
            .Call(C_glmm_score_bgen13, as.numeric(res), obj$P, infile, paste0(outfile, "_tmp.", i), center2, MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select, threadInfo$begin[i], threadInfo$end[i], threadInfo$pos[i], bgenInfo$N, bgenInfo$CompressionFlag, 1)
          } else {
            .Call(C_glmm_score_bgen11, as.numeric(res), obj$P, infile, paste0(outfile, "_tmp.", i), center2, MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select, threadInfo$begin[i], threadInfo$end[i], threadInfo$pos[i], bgenInfo$N, bgenInfo$CompressionFlag, 1)
          }
        } else {
          if (bgenInfo$LayoutFlag == 2) {
            .Call(C_glmm_score_bgen13_sp, as.numeric(res), obj$Sigma_i, obj$Sigma_iX, obj$cov, infile, paste0(outfile, "_tmp.", i), center2, MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select, threadInfo$begin[i], threadInfo$end[i], threadInfo$pos[i], bgenInfo$N, bgenInfo$CompressionFlag, 1)
          } else {
            .Call(C_glmm_score_bgen11_sp, as.numeric(res), obj$Sigma_i, obj$Sigma_iX, obj$cov, infile, paste0(outfile, "_tmp.", i), center2, MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select, threadInfo$begin[i], threadInfo$end[i], threadInfo$pos[i], bgenInfo$N, bgenInfo$CompressionFlag, 1)
          }
        }
      }
      outTmp <- file(outfile, "w")
      writeLines("SNP\tRSID\tCHR\tPOS\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL", outTmp)
      for(i in 1:ncores){
        inTmp <- readLines(paste0(outfile, "_tmp.", i))
        writeLines(inTmp, outTmp)
        unlink(paste0(outfile, "_tmp.", i))
      }
      close(outTmp)
    } else {
      if(!is.null(obj$P)) {
        if (bgenInfo$LayoutFlag == 2) {
          .Call(C_glmm_score_bgen13, as.numeric(res), obj$P, infile, outfile, center2, MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select, 0, bgenInfo$M, bgenInfo$offset + 4, bgenInfo$N, bgenInfo$CompressionFlag, 0)
        } else{
          .Call(C_glmm_score_bgen11, as.numeric(res), obj$P, infile, outfile, center2, MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select, 0, bgenInfo$M, bgenInfo$offset + 4, bgenInfo$N, bgenInfo$CompressionFlag, 0)
        }
      } else {
        if (bgenInfo$LayoutFlag == 2) {
          .Call(C_glmm_score_bgen13_sp, as.numeric(res), obj$Sigma_i, obj$Sigma_iX, obj$cov, infile, outfile, center2, MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select, 0, bgenInfo$M, bgenInfo$offset + 4, bgenInfo$N, bgenInfo$CompressionFlag, 0)
        } else{
          .Call(C_glmm_score_bgen11_sp, as.numeric(res), obj$Sigma_i, obj$Sigma_iX, obj$cov, infile, outfile, center2, MAF.range[1], MAF.range[2], miss.cutoff, miss.method, nperbatch, select, 0, bgenInfo$M, bgenInfo$offset + 4, bgenInfo$N, bgenInfo$CompressionFlag, 0)
        }
      }
    }
    
    return(invisible(NULL))
    
  } else { # text genotype files
    if(ncores != 1) stop("Error: parallel computing currently not implemented for plain text format genotypes.")
    if(class(obj) == "glmmkin.multi") stop("Error: multiple phenotypes currently not implemented for plain text format genotypes.")
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
    if(!is.numeric(infile.ncol.skip) | infile.ncol.skip < 0)
      stop("Error: number of skipped cols of the input file is incorrect!")
    if(length(infile.ncol.print) != length(infile.header.print))
      stop("Error: number of cols selected to print does not match number of header names!")
    if(is.null(infile.ncol.print))
      infile.ncol.print <- 0
    if(is.null(infile.header.print))
      infile.header.print <- infile.na
    if(any(!is.numeric(infile.ncol.print)) | any(infile.ncol.print < 0) | any(infile.ncol.print > infile.ncol.skip))
      stop("Error: cols selected to print have incorrect indices!")
    if(any(infile.ncol.print != sort(infile.ncol.print)))
      stop("Error: col indices must be sorted increasingly in infile.ncol.print!")
    if(is.null(select)) {
      warning("Argument select is unspecified... Assuming the order of individuals in infile matches unique id_include in obj...")
      select <- 1:length(unique(obj$id_include))
    }
    select2 <- select[select > 0]
    if(any(duplicated(select2)) || max(select2) > length(unique(obj$id_include))) stop("Error: select is a vector of orders, individuals not in obj should be coded 0!")
    res <- res[select2]
    if(!is.null(obj$P)) obj$P <- obj$P[select2, select2]
    else {
      obj$Sigma_iX <- Matrix(obj$Sigma_iX[select2, , drop = FALSE], sparse = TRUE)
      obj$Sigma_i <- Matrix(obj$Sigma_i[select2, select2], sparse = TRUE)
      obj$cov <- Matrix(obj$cov, sparse = TRUE)
    }
    select[select > 0] <- 1:sum(select > 0)
    rm(select2)
    if(center) {
      if(!is.null(obj$P)) time <- .Call(C_glmm_score_text, res, obj$P, infile, outfile, tol, 'c', MAF.range[1], MAF.range[2], miss.cutoff, miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, infile.header.print, nperbatch, select)
      else time <- .Call(C_glmm_score_text_sp, res, obj$Sigma_i, obj$Sigma_iX, obj$cov, infile, outfile, tol, 'c', MAF.range[1], MAF.range[2], miss.cutoff, miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, infile.header.print, nperbatch, select)
    } else {
      if(!is.null(obj$P)) time <- .Call(C_glmm_score_text, res, obj$P, infile, outfile, tol, 'n', MAF.range[1], MAF.range[2], miss.cutoff, miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, infile.header.print, nperbatch, select)
      else time <- .Call(C_glmm_score_text_sp, res, obj$Sigma_i, obj$Sigma_iX, obj$cov, infile, outfile, tol, 'n', MAF.range[1], MAF.range[2], miss.cutoff, miss.method, infile.nrow, infile.nrow.skip, infile.sep, infile.na, infile.ncol.skip, infile.ncol.print, infile.header.print, nperbatch, select)
    }
    #print(sprintf("Computational time: %.2f seconds", time))
    return(invisible(time))
  }
}

glmm.score.meta <- function(files, outfile, SNP = rep("SNP", length(files)), A1 = rep("A1", length(files)), A2 = rep("A2", length(files))) {
  k <- length(files)
  if(length(SNP) != k) stop("Error: \"SNP\" must have the same length as \"files\"!")
  if(length(A1) != k) stop("Error: \"A1\" must have the same length as \"files\"!")
  if(length(A2) != k) stop("Error: \"A2\" must have the same length as \"files\"!")
  master <- read.table(files[1], header=T, as.is=T)[, c(SNP[1], A1[1], A2[1], "N", "AF", "SCORE", "VAR", "PVAL")]
  names(master)[1:3] <- c("SNP", "A1", "A2")
  master <- master[!is.na(master$SCORE) & !is.na(master$VAR) & !is.na(master$PVAL), ]
  flag <- rep(0, nrow(master))
  if(k > 1) {
    for(i in 2:k) {
      tmp <- read.table(files[i], header=T, as.is=T)[, c(SNP[i], A1[i], A2[i], "N", "AF", "SCORE", "VAR", "PVAL")]
      names(tmp)[1:3] <- c("SNP", "A1", "A2")
      tmp <- tmp[!is.na(tmp$SCORE) & !is.na(tmp$VAR) & !is.na(tmp$PVAL), ]
      idx <- tmp$SNP %in% master$SNP
      if(sum(!idx) > 0) {
        flag <- c(flag, rep(0, sum(!idx)))
        master <- rbind(master, tmp[!idx, ])
      }
      idx2 <- match(tmp$SNP[idx], master$SNP)
      noflip <- master$A1[idx2] == tmp$A1[idx] & master$A2[idx2] == tmp$A2[idx]
      flip <- master$A1[idx2] == tmp$A2[idx] & master$A2[idx2] == tmp$A1[idx]
      flag[idx2] <- flag[idx2] + as.numeric(!noflip & !flip)
      master$AF[idx2][noflip] <- (master$AF[idx2][noflip]*master$N[idx2][noflip] + tmp$AF[idx][noflip]*tmp$N[idx][noflip])/(master$N[idx2][noflip] + tmp$N[idx][noflip])
      master$AF[idx2][flip] <- (master$AF[idx2][flip]*master$N[idx2][flip] + (1 - tmp$AF[idx][flip])*tmp$N[idx][flip])/(master$N[idx2][flip] + tmp$N[idx][flip])
      master$N[idx2] <- master$N[idx2] + tmp$N[idx]
      master$SCORE[idx2][noflip] <- master$SCORE[idx2][noflip] + tmp$SCORE[idx][noflip]
      master$SCORE[idx2][flip] <- master$SCORE[idx2][flip] - tmp$SCORE[idx][flip]
      master$VAR[idx2] <- master$VAR[idx2] + tmp$VAR[idx]
    }
    if(any(flag > 0)) {
      cat("The following SNPs have been removed due to inconsistent alleles across studies:\n")
      print(master$SNP[flag > 0])
      master <- subset(master, flag == 0)
    }
    master$PVAL <- pchisq(master$SCORE^2/master$VAR, 1, lower.tail=F)
  }
  write.table(master, outfile, sep="\t", row.names=F, col.names=T, quote=F)
  invisible(master)
}
