library(MASS)
library(edgeR)
library(DESeq2)
library(sSeq)
library(EBSeq)

paired_estimate_params <- function(rawdata, Patient, Tissue) {
	y <- DGEList(counts = rawdata)
	#y <- calcNormFactors(y)
	data.frame(Sample = colnames(y), Patient, Tissue)
	design <- model.matrix(~Patient + Tissue)
	rownames(design) <- colnames(y)

	dispCoxReidInterpolateTagwise(y$counts, design = design, offset = getOffset(y), dispersion = .1, trend = FALSE, AveLogCPM = NULL, min.row.sum = 5, prior.df = 0, span = 0.3, grid.npts = 15, grid.range = c(-8, 8)) -> dispsCR

	sample_data <- data.frame(Patient, Tissue)
	sample_data$libsize <- log(colSums(y$counts))
	libsize <- sample_data$libsize
	nofit <- 1000000
	fc <- matrix(nrow = dim(y$counts)[1], ncol = length(levels(Patient)) + length(levels(Tissue)) - 1)
	for (i in 1:dim(y$counts)[1]) {
		f <- negative.binomial(link = "log", theta = 1 / dispsCR[i])
		tryCatch({ glm(y$counts[i,] ~ Patient + Tissue + 0, offset = libsize, family = f) -> fit },
				 warning = function(w) { assign('nofit', c(nofit, i), parent.env(environment())) })
		fc[i,] <- fit$coefficients
	}
	y <- DGEList(counts = rawdata[-nofit,])
	list(y = y, fc = fc[-nofit,], dispsCR = dispsCR[-nofit], sample_data = sample_data, nofit = nofit)
}

paired_DE_call <- function(rawdata, Patient, Tissue) {

	dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = data.frame(Patient, Tissue), design = ~Patient + Tissue)
	dds <- estimateSizeFactors(dds)
	dds <- estimateDispersions(dds)
	dds <- nbinomWaldTest(dds)
	res <- results(dds)
	pval <- res$pval
	padj <- res$padj
	res <- cbind(pval, padj)
	ds2 <- as.matrix(res)
	rm(res, pval, padj)

	#DESeq#
	DESeq_cds <- newCountDataSet(rawdata, data.frame(Patient, Tissue))
	DESeq_cds <- estimateSizeFactors(DESeq_cds)
	DESeq_cds <- estimateDispersions(DESeq_cds, method = "pooled-CR", modelFrame = data.frame(Patient, Tissue), modelFormula = count ~ Tissue + Patient)
	fit1 <- fitNbinomGLMs(DESeq_cds, count ~ Tissue + Patient)
	fit0 <- fitNbinomGLMs(DESeq_cds, count ~ Patient)
	pval <- nbinomGLMTest(fit1, fit0)
	padj <- p.adjust(pval, method = "BH")
	res <- cbind(pval, padj)
	ds <- as.matrix(res)
	rm(res, pval, padj)

	#edgeR#
	design <- model.matrix(~Patient + Tissue)
	edgeR_cds <- DGEList(rawdata)
	edgeR_cds <- calcNormFactors(edgeR_cds)
	edgeR_cds <- estimateGLMCommonDisp(edgeR_cds, design)
	edgeR_cds <- estimateGLMTrendedDisp(edgeR_cds, design)
	edgeR_cds <- estimateGLMTagwiseDisp(edgeR_cds, design)
	fit <- glmFit(edgeR_cds, design)
	res <- glmLRT(fit)$table
	pval <- res$PValue
	padj <- p.adjust(pval, method = "BH")
	res <- cbind(pval, padj)
	er <- as.matrix(res)
	rm(res, pval, padj)

	#sSeq#

	as.character(Tissue) -> sSeq_condition
	temp <- order(sSeq_condition)
	sSeq_condition <- sSeq_condition[temp]
	sSeq_Patient <- Patient[temp]
	sSeq_rawdata <- rawdata[, temp]

	res <- nbTestSH(sSeq_rawdata, sSeq_condition, condA = unique(sSeq_condition)[1], condB = unique(sSeq_condition)[2], coLevels = data.frame(sSeq_Patient), pairedDesign = TRUE, pairedDesign.dispMethod = "pooled")
	pval <- res$pval
	padj <- p.adjust(pval, method = "BH")
	res <- cbind(pval, padj)
	ss <- as.matrix(res)
	rm(res, pval, padj)


	packages <- c("ds2", "er", "ss")

	de <- rep(TRUE, dim(rawdata)[1])
	for (i in packages) {
		temp <- length(which(get(i)[, "padj"] < 0.05))
		print(paste(i, ": number of DE called", temp))
		de <- de & get(i)[, "padj"] < 0.05
	}
	print(paste("intesection :", length(which(de))))
	de
}


of_estimate_params <- function(rawdata, condition) {
	y <- DGEList(counts = rawdata)
	y <- calcNormFactors(y)

	design <- model.matrix(~condition)
	rownames(design) <- colnames(y)

	dispCoxReidInterpolateTagwise(y$counts, design = design, offset = getOffset(y), dispersion = .1, trend = FALSE, AveLogCPM = NULL, min.row.sum = 5, prior.df = 0, span = 0.3, grid.npts = 15, grid.range = c(-8, 8)) -> dispsCR

	sample_data <- data.frame(condition)
	#sample_data$libsize = log(colSums(y$counts)) - log(nrow(y$counts))
	sample_data$libsize <- log(colSums(y$counts))
	libsize <- sample_data$libsize
	nofit <- 1000000
	fc <- matrix(nrow = dim(y$counts)[1], ncol = 2)
	for (i in 1:dim(y$counts)[1]) {
		f <- negative.binomial(link = "log", theta = 1 / dispsCR[i])
		tryCatch({ glm(y$counts[i,] ~ condition + 0, offset = libsize, family = f) -> fit },
				 warning = function(w) { assign('nofit', c(nofit, i), parent.env(environment())) })
		fc[i,] <- fit$coefficients
	}
	y <- DGEList(counts = rawdata[-nofit,])
	list(y = y, fc = fc[-nofit,], dispsCR = dispsCR[-nofit], sample_data = sample_data, nofit = nofit)
}


zf_estimate_params <- function(rawdata) {
	y <- DGEList(counts = rawdata)
	y <- calcNormFactors(y)

	design <- matrix(1, dim(rawdata)[2], 1)

	#rownames(design) <- colnames(y)

	dispCoxReidInterpolateTagwise(y$counts, design = design, offset = getOffset(y), dispersion = .1, trend = FALSE, AveLogCPM = NULL, min.row.sum = 5, prior.df = 0, span = 0.3, grid.npts = 15, grid.range = c(-8, 8)) -> dispsCR
	#span and prior.df are irrelevant with trend=FALSE

	list(y = y, dispsCR = dispsCR, libsize = log(colSums(y$counts)))
}


of_DE_call <- function(rawdata, condition) {
	#DESeq2#
	dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = data.frame(condition), design = ~condition)
	dds <- estimateSizeFactors(dds)
	dds <- estimateDispersions(dds)
	dds <- nbinomWaldTest(dds)
	res <- results(dds)
	pval <- res$pval
	padj <- res$padj
	res <- cbind(pval, padj)
	ds2 <- as.matrix(res)
	rm(res, pval, padj)

	# #DESeq#
	# DESeq_cds <- newCountDataSet(rawdata, condition)
	# DESeq_cds <- estimateSizeFactors(DESeq_cds)
	# DESeq_cds <- estimateDispersions(DESeq_cds)
	# pval <- nbinomTest(DESeq_cds, unique(condition)[1], unique(condition)[2], pvals_only = TRUE)
	# padj <- p.adjust(pval, method = "BH")
	# res <- cbind(pval, padj)
	# ds <- as.matrix(res)
	# rm(res, pval, padj)

	#edgeR#
	edgeR_cds <- DGEList(rawdata, group = condition)
	edgeR_cds <- calcNormFactors(edgeR_cds)
	edgeR_cds <- estimateCommonDisp(edgeR_cds)
	edgeR_cds <- estimateTagwiseDisp(edgeR_cds)
	res <- exactTest(edgeR_cds, pair = c(unique(condition)[1], unique(condition)[2]))$table
	pval <- res$PValue
	padj <- p.adjust(pval, method = "BH")
	res <- cbind(pval, padj)
	er <- as.matrix(res)
	rm(res, pval, padj)

	#sSeq#
	as.character(condition) -> sSeq_condition
	res <- nbTestSH(rawdata, sSeq_condition, condA = unique(sSeq_condition)[1], condB = unique(sSeq_condition)[2])
	pval <- res$pval
	padj <- p.adjust(pval, method = "BH")
	res <- cbind(pval, padj)
	ss <- as.matrix(res)
	rm(res, pval, padj)

	#NPEBSeq#
	#G1data <- rawdata[,which(condition==levels(condition)[1])]
	#G2data <- rawdata[,which(condition==levels(condition)[2])]
	#maxid1<-which.max(colSums(G1data))
	#maxid2<-which.max(colSums(G2data))
	#Q1<-compu_prior(G1data[,maxid1],maxiter=100,grid.length=1000)
	#Q2<-compu_prior(G2data[,maxid2],maxiter=3000,grid.length=1000)
	#resg<-NPEBSeq_biordf(G1data,G2data,Q1,Q2)

	#EBSeq

	Sizes <- MedianNorm(rawdata)
	EBOut <- EBTest(Data = rawdata, Conditions = condition, sizeFactors = Sizes, maxround = 5)
	data.frame(pval = 1 - GetPP(EBOut)) -> temp0
	temp1 <- rawdata
	merge(temp1, temp0, all.x = TRUE, by.x = 0, by.y = 0) -> temp2
	pval <- temp2[, "pval"]
	names(pval) <- temp2[, "Row.names"]
	pval <- pval[rownames(rawdata)]
	padj <- pval
	res <- cbind(pval, padj)
	eb <- as.matrix(res)
	rm(res, pval, padj)


	#AMAP.Seq#
	#mydata = RNASeq.Data(rawdata, size=Norm.GMedian(rawdata), group = sSeq_condition)
	#decom.est=MGN.EM(mydata,nK=3,p0=NULL,d0=0,iter.max=10,nK0=3)
	#res=test.AMAP(mydata, MGN=decom.est$MGN,FC=1.0)
	#pval = res$prob
	#padj = res$fdr
	#res = cbind(pval, padj)
	#am <- as.matrix(res)
	#rm(res, pval, padj)

	#packages = c("ds2", "ds","er","ss","eb")
	packages <- c("ds2", "er", "ss", "eb")

	de <- rep(TRUE, dim(rawdata)[1])
	for (i in packages) {
		temp <- length(which(get(i)[, "padj"] < 0.05))
		print(paste(i, ": number of DE called", temp))
		de <- de & get(i)[, "padj"] < 0.05
	}
	print(paste("intersection :", length(which(de))))
	de[is.na(de)] <- FALSE
	de
}

print_params <- function(params, de) {
	y <- params$y
	cat("number of genes:\t")
	cat(dim(y$counts)[1])
	cat("\n")
	temp <- cpm(y, log = TRUE, prior.count = 1)
	cat("cpm")
	cat("\t")
	cat(signif(median(temp), digits = 3))
	cat(" (")
	cat(signif(quantile(temp, 0.25), digits = 3))
	cat(" - ")
	cat(signif(quantile(temp, 0.75), digits = 3))
	cat(")\n")


	temp <- params$dispsCR
	cat("dispersion\t")
	cat(signif(median(temp), digits = 3))
	cat(" (")
	cat(signif(quantile(temp, 0.25), digits = 3))
	cat(" - ")
	cat(signif(quantile(temp, 0.75), digits = 3))
	cat(")\n")

	if (!is.null(params$fc)) {
		fc <- params$fc
		if (dim(fc)[2] == 2) {
			temp <- log2(exp(abs(fc[, 1] - fc[, 2])))
		} else {
			temp <- log2(exp(abs(fc[, dim(fc)[2]])))
		}
		temp <- temp[de]
		cat("fc\t")
		cat(signif(median(temp), digits = 3))
		cat(" (")
		cat(signif(quantile(temp, 0.25), digits = 3))
		cat(" - ")
		cat(signif(quantile(temp, 0.75), digits = 3))
		cat(")\n")
	}

	cat("libsize\t")
	if (is.null(params$sample_data)) {
		libsize <- params$libsize
	} else {
		libsize <- params$sample_data$libsize
	}
	libsize <- log10(exp(libsize))
	cat(signif(mean(libsize), digits = 3))
	cat(" +/- ")
	cat(signif(var(libsize), digits = 3))
	cat("\n")
}


#input fitted parameters from preliminary (paired) data, output NB based on fitted parameters
#disp - vector of estimated dispersions
#libsize - library sizes
#fc_patient - matrix of fold change related to patient factor
#fc - vector, fold change related to factor of interest
#n - number of samples
paired_simdata <- function(disps, libsizes, fc_patient, fc, n, de, randomseed = 1) {
	set.seed(randomseed)

	#set fc truth data
	fc[!de] <- 0

	#generate library sizes for n samples based around the mean of the library sizes in the preliminary data
	sim_libsizes <- runif(n * 2, min = min(libsizes), max = max(libsizes))

	#for each n and gene, generate a fc related to the sample
	patient_min <- apply(fc_patient, 1, min)
	patient_max <- apply(fc_patient, 1, max)

	#generate patient coefficients related to the n samples for each gene
	sim_pfc <- t(sapply(1:length(patient_min), function(i) {
		runif(n, min = patient_min[i], max = patient_max[i])
	}))

	m <- matrix(nrow = length(disps), ncol = n * 2)

	for (i in 1:length(disps)) {
		for (j in 1:(n * 2)) {
			m[i, j] <- rnbinom(1, mu = exp(sim_libsizes[j] +
											  sim_pfc[i, ceiling(j / 2)] +
											  ifelse(j %% 2 == 0, fc[i], 0)), size = 1 / disps[i])
		}
	}
	#m[m >= 100000] <- 100000
	Patient <- sort(c(1:n, 1:n))
	Tissue <- rep(c("N", "T"), n)
	colnames(m) <- paste0("s", Patient, Tissue)
	rownames(m) <- paste0("g", 1:length(disps))
	m
}


#input fitted parameters from preliminary (paired) data, output NB based on fitted parameters
#disp - vector of estimated dispersions
#libsize - library sizes
#fc_patient - matrix of fold change related to patient factor
#fc - vector, fold change related to factor of interest
#n - number of samples
of_simdata <- function(disps, libsizes, fc, n, de, randomseed = 1) {
	set.seed(randomseed)
	mean_expr <- (fc[, 1] + fc[, 2]) / 2
	#set fc truth data
	#fc[!de,] <- c(mean_expr[!de],mean_expr[!de])

	#generate library sizes for n samples based around the mean of the library sizes in the preliminary data
	sim_libsizes <- runif(n * 2, min = min(libsizes), max = max(libsizes))

	m <- matrix(nrow = length(disps), ncol = n * 2)

	for (i in 1:length(disps)) {
		for (j in 1:(n * 2)) {
			m[i, j] <- rnbinom(1, mu = exp(sim_libsizes[j] + ifelse(j <= n, fc[i, 1], fc[i, 2])), size = 1 / disps[i])
		}
	}
	#m[m >= 10000] <- 10000
	label <- paste0(c(rep("A", n), rep("B", n)), rep(1:n, 2))
	colnames(m) <- label
	rownames(m) <- paste0("g", 1:length(disps))
	m
}


printerror <- function(e, program) {
	print(e)
	print(program)
	print(paste(n, mean(sample_data$libsize), randomseed))
}


of_eval <- function(rawdata, condition, program) {
	pval_list <- list()

	if (program == "DESeq2") {
		#DESeq2#
		tryCatch({
			dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = data.frame(condition), design = ~condition)
			dds <- estimateSizeFactors(dds)
			dds <- estimateDispersions(dds)
			dds <- nbinomWaldTest(dds)
			res <- results(dds)
			pval <- res$pval
			padj <- res$padj
			res <- cbind(pval, padj)
			pval_list[["ds2"]] <- as.matrix(res)
			rm(dds, res, pval, padj)
			gc()
		}, error = function(e) { printerror(e, "DESeq2") })
	}
	#     else if (program == "DESeq") {
	# 	#DESeq#
	# 	tryCatch({
	# 		DESeq_cds <- newCountDataSet(rawdata, condition)
	# 		DESeq_cds <- estimateSizeFactors(DESeq_cds)
	# 		DESeq_cds <- estimateDispersions(DESeq_cds)
	# 		pval <- nbinomTest(DESeq_cds, "A", "B", pvals_only = TRUE)
	# 		padj <- p.adjust(pval, method = "BH")
	# 		res <- cbind(pval, padj)
	# 		pval_list[["ds"]] <- as.matrix(res)
	# 		rm(DESeq_cds, res, pval, padj)
	# 		gc()
	# 	}, error = function(e) { printerror(e, "DESeq") })
	# }
	else if (program == "edgeR") {
		#edgeR#
		tryCatch({
			edgeR_cds <- DGEList(rawdata, group = condition)
			edgeR_cds <- calcNormFactors(edgeR_cds)
			edgeR_cds <- estimateCommonDisp(edgeR_cds)
			edgeR_cds <- estimateTagwiseDisp(edgeR_cds)
			res <- exactTest(edgeR_cds, pair = c("A", "B"))$table
			pval <- res$PValue
			padj <- p.adjust(pval, method = "BH")
			res <- cbind(pval, padj)
			pval_list[["er"]] <- as.matrix(res)
			rm(edgeR_cds, res, pval, padj)
			gc()
		}, error = function(e) { printerror(e, "edgeR") })
	} else if (program == "sSeq") {
		#sSeq#
		tryCatch({
			as.character(condition) -> sSeq_condition
			res <- nbTestSH(rawdata, sSeq_condition, condA = unique(sSeq_condition)[1], condB = unique(sSeq_condition)[2])
			pval <- res$pval
			padj <- p.adjust(pval, method = "BH")
			res <- cbind(pval, padj)
			pval_list[["ss"]] <- as.matrix(res)
			rm(res, sSeq_condition, pval, padj)
			gc()
		}, error = function(e) { printerror(e, "sSeq") })
	} else if (program == "EBSeq") {
		#EBSeq
		tryCatch({
			Sizes <- MedianNorm(rawdata)
			EBOut <- EBTest(Data = rawdata, Conditions = condition, sizeFactors = Sizes, maxround = 5)
			data.frame(pval = 1 - GetPP(EBOut)) -> temp0
			temp1 <- rawdata
			merge(temp1, temp0, all.x = TRUE, by.x = 0, by.y = 0) -> temp2
			pval <- temp2[, "pval"]
			names(pval) <- temp2[, "Row.names"]
			pval <- pval[rownames(rawdata)]
			padj <- pval
			res <- cbind(pval, padj)
			pval_list[["eb"]] <- as.matrix(res)
			rm(temp0, temp1, temp2, EBOut, Sizes, res, pval, padj)
			gc()
		}, error = function(e) { printerror(e, "EBSeq") })
	} else { stop("please select a program: DESeq2, edgeR, EBSeq or sSeq") }
	pval_list
}


paired_eval <- function(rawdata, Patient, Tissue) {
	pval_list <- list()
	if (program == "DESeq2") {
		#DESeq2#
		tryCatch({
			dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = data.frame(Patient, Tissue), design = ~Patient + Tissue)
			dds <- estimateSizeFactors(dds)
			dds <- estimateDispersions(dds)
			dds <- nbinomWaldTest(dds)
			res <- results(dds)
			pval <- res$pval
			padj <- res$padj
			res <- cbind(pval, padj)
			pval_list[["ds2"]] <- as.matrix(res)
			rm(res, pval, padj)
			gc()
		}, error = function(e) { printerror(e, "DESeq2") })
	}
	#     else if (program == "DESeq") {
	# 	#DESeq#
	# 	tryCatch({
	# 		DESeq_cds <- newCountDataSet(rawdata, data.frame(Patient, Tissue))
	# 		DESeq_cds <- estimateSizeFactors(DESeq_cds)
	# 		DESeq_cds <- estimateDispersions(DESeq_cds, method = "pooled-CR", modelFrame = data.frame(Patient, Tissue), modelFormula = count ~ Tissue + Patient)
	# 		fit1 <- fitNbinomGLMs(DESeq_cds, count ~ Tissue + Patient)
	# 		fit0 <- fitNbinomGLMs(DESeq_cds, count ~ Patient)
	# 		pval <- nbinomGLMTest(fit1, fit0)
	# 		padj <- p.adjust(pval, method = "BH")
	# 		res <- cbind(pval, padj)
	# 		pval_list[["ds"]] <- as.matrix(res)
	# 		rm(res, pval, padj)
	# 		gc()
	# 	}, error = function(e) { printerror(e, "DESeq") })
	# }
	else if (program == "edgeR") {
		#edgeR#
		tryCatch({
			design <- model.matrix(~Patient + Tissue)
			edgeR_cds <- DGEList(rawdata)
			edgeR_cds <- calcNormFactors(edgeR_cds)
			edgeR_cds <- estimateGLMCommonDisp(edgeR_cds, design)
			edgeR_cds <- estimateGLMTrendedDisp(edgeR_cds, design)
			edgeR_cds <- estimateGLMTagwiseDisp(edgeR_cds, design)
			fit <- glmFit(edgeR_cds, design)
			res <- glmLRT(fit)$table
			pval <- res$PValue
			padj <- p.adjust(pval, method = "BH")
			res <- cbind(pval, padj)
			pval_list[["er"]] <- as.matrix(res)
			rm(res, pval, padj)
			gc()
		}, error = function(e) { printerror(e, "edgeR") })
	} else if (program == "sSeq") {
		#sSeq#
		tryCatch({
			as.character(Tissue) -> sSeq_condition
			temp <- order(sSeq_condition)
			sSeq_condition <- sSeq_condition[temp]
			sSeq_Patient <- Patient[temp]
			sSeq_rawdata <- rawdata[, temp]
			res <- nbTestSH(sSeq_rawdata, sSeq_condition, condA = unique(sSeq_condition)[1], condB = unique(sSeq_condition)[2], coLevels = data.frame(sSeq_Patient), pairedDesign = TRUE, pairedDesign.dispMethod = "pooled")
			pval <- res$pval
			padj <- p.adjust(pval, method = "BH")
			res <- cbind(pval, padj)
			pval_list[["ss"]] <- as.matrix(res)
			rm(res, pval, padj)
			gc()
		}, error = function(e) { printerror(e, "edgeR") })
	} else { stop("please select a program: DESeq2, DESeq, edgeR or sSeq") }
	pval_list
}


#estimate params for preliminary data
#designtype = "paired" or "one factor"
#rawdata = matrix of count data
#conditions = condition of samples, in order of columns of rawdata
##for one factor, vector of two conditions
#pairing = pairing information for paired experiments only
#returns: list of values required for simulation
##
estimate_params <- function(rawdata, condition, pairing = NULL, designtype) {
	condition[condition == unique(condition)[1]] <- "A"
	condition[condition == unique(condition)[2]] <- "B"
	sort_order <- order(condition)
	condition <- condition[sort_order]
	rawdata <- rawdata[, sort_order]
	if (designtype == "one factor") {
		params <- of_estimate_params(rawdata, condition)
		de <- of_DE_call(rawdata[-params$nofit,], condition)
		params[["de"]] <- de
	} else if (designtype == "paired") {
		pairing <- pairing[sort_order]
		params <- paired_estimate_params(rawdata, pairing, condition)
		de <- paired_DE_call(rawdata[-params$nofit], pairing, condition)
		params[["de"]] <- de
	} else { stop("choose designtype, paired or one factor") }
	#list(y=y, fc=fc[-nofit,], dispsCR = dispsCR[-nofit], sample_data=sample_data, nofit=nofit)
	params
}

#budget = max budget | 3000
#per_sample_price = preperation price per sample | 241
#lane_size = size of lane | 150000000;
#lane_price = lane price | 1331;
#mapping_proportion = proportion of genes that map to transcriptome after filtering | 0.2;
#sims = number of simulations to run | 5
#params = estimation of parameters from estimate_params
#program = program to estimate power
RS_simulation <- function(budget = 3000, per_sample_price = 241, lane_size = 150e6, lane_price = 1331, mapping_proportion = 0.2,
						  sims = 5, params, designtype = "one factor", nmax = 20, nmin = 2, program) {
	attach(params)
	result_matrix <- matrix(ncol = sims, nrow = 0)
	n <- max(2, nmin)
	meanlibsize <- (lane_size * mapping_proportion) * (budget - (2 * n * per_sample_price)) / lane_price;
	while (meanlibsize > 100000 & n <= nmax) {
		results <- vector()
		for (randomseed in 1:sims) {
			#if meanlibsize specified (not zero) then change the mean libsize in the sample_data
			if (meanlibsize != "0") {
				temp <- log(meanlibsize)
				sample_data$libsize <- sample_data$libsize - (mean(max(sample_data$libsize), min(sample_data$libsize)) - temp)
			}
			if (designtype == "one factor") {
				m_sim <- of_simdata(disps = dispsCR, libsizes = sample_data$libsize, fc = fc, n = n, de, randomseed = randomseed)
				condition <- c(rep("A", n), rep("B", n))
				suppressMessages(pval_list <- of_eval(m_sim, condition, program))
			} else if (designtype == "paired") {
				m_sim <- paired_simdata(disps = dispsCR, libsizes = sample_data$libsize, fc_patient = fc[, 1:(dim(fc)[2] - 1)], fc = fc[, dim(fc)[2]], n = n, de, randomseed = randomseed)
				Patient <- factor(sort(c(1:n, 1:n)))
				Tissue <- factor(rep(c("N", "T"), n))
				suppressMessages(pval_list <- paired_eval(m_sim, Patient, Tissue))
			}
			n_pos <- length(which(de & !is.na(pval_list[[1]])))
			n_tp <- length(which(pval_list[[1]] < 0.05 & de))
			results[randomseed] <- n_tp / n_pos
		}
		result_matrix <- rbind(result_matrix, results)
		rownames(result_matrix)[dim(result_matrix)[1]] <- n
		n <- n + 1
		meanlibsize <- (lane_size * mapping_proportion) * (budget - (2 * n * per_sample_price)) / lane_price;
	}
	detach(params)
	colnames(result_matrix) <- 1:sims
	result_matrix
}
