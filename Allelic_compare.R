######################################################
# Code 2: Allelic analysis of GSE147507 RNA-seq data
######################################################
# Code 2.1: comparison within each series
dat <- read.table("GSE147507_stage1/allelicSNPgloabl_GSE147507.txt.gz", header=F)
snpQuery <- dat[dat[,4] + dat[,5] >= 7, ]   # Depath >= 7

# series <- c("Series2.A549", "Series3.A549", "Series4.A549", "Series8.A549", "Series9.NHBE")
query <- c("Series1.NHBE", "Series2", "Series3", "Series4", "Series5", "Series6", "Series7", "Series8", "Series9.NHBE$", "Series16")
for (series in query){
	indexed <- grep(series, snpQuery[,10])
	snpSeries <- snpQuery[indexed, ]

	snp_stat <- tapply(droplevels(snpSeries[,9]), droplevels(snpSeries[,1]), summary)
	num <- length(snp_stat)

	snp_kept <- vector()
	for (x in 1:num){
		snp_count <- snp_stat[[eval(names(snp_stat)[x])]]
		count = 0
		cells <- names(snp_count)
		for (m in 1:length(cells)){
			if (snp_count[cells[m]] >=2){
				count = count + 1
			}
		}
	
		# at least 2 conditions
		if (count >= 2){
			snp_kept <- append(snp_kept, names(snp_stat)[x])
		}
	}
	
	series <- sub("\\$", "", series)
	print(paste("Number of SNP retained:", length(snp_kept), "in", series, sep=" "))
	
	snpKeep <- subset(snpSeries, snpSeries[,1] %in% snp_kept)
	outRdata <- paste("GSE147507_stage1/AllelicRatio_SNPkeep_", series, "_GSE147507.RData", sep="")
	save(snpKeep, file = outRdata)

	load("allelicFunction.RData")

	snpID <- as.character(unique(snpKeep[,1]))
	size <- length(snpID)

	# Kruskal-Wallis rank sum test
	compResults <- matrix(nrow = size, ncol = 6)
	colnames(compResults) <- c("chi2_ratio", "p_ratio", "chi2_fc", "p_fc", "chi2_directed", "p_directed")
	rownames(compResults) <- snpID

	for (m in 1:size){
		eachTest <- subset(snpKeep, snpKeep[,1] %in% snpID[m])
		cell_group = eachTest[,9]
		
		allelicRatio = apply(eachTest[,4:5], 1, Adjust_Allelic_Ratio)
		res_anova <- kruskal.test(allelicRatio ~ cell_group)
		compResults[m, 1] = res_anova$statistic
		compResults[m, 2] = res_anova$p.value
	
		fc = apply(eachTest[,4:5], 1, Allelic_Fold_Change)
		res_anova_fc <- kruskal.test(fc ~ cell_group)
		compResults[m, 3] = res_anova_fc$statistic
		compResults[m, 4] = res_anova_fc$p.value
	
		direct = apply(eachTest[,4:5], 1, RefvsAlt_Ratio)
		res_anova_direct <- kruskal.test(direct ~ cell_group)
		compResults[m, 5] = res_anova_direct$statistic
		compResults[m, 6] = res_anova_direct$p.value	
	}
	outTest <- paste("GSE147507_stage1/resultKruskal_Allelic_", series, "_GSE147507.txt", sep="")
	write.table(compResults, file = outTest, row.names=T, col.names=T, quote=FALSE, sep="\t")

	sigP05 <- subset(compResults, compResults[,2] < 0.051 & compResults[,4] < 0.051 & compResults[,6] < 0.051)
	outSig <- paste("GSE147507_stage1/resultKruskal_Allelic_", series, "_GSE147507_P05.txt", sep="")
	write.table(sigP05, file = outSig, row.names=T, col.names=T, quote=FALSE, sep="\t")
}
