###############################
# Code 1: Process meta-data
###############################
dat <- read.table("GSE147507_stage1/dataInfo_GSE147507.txt", header=F, sep="\t")

demo <- dat[!duplicated(dat[,3]), c(3:5,9)]

library("GEOquery")

gse = getGEO(filename="GSE147507-GPL18573_series_matrix.txt.gz")

info <- pData(phenoData(gse))

matched <- match(demo[,1], rownames(info))
demo_add <- cbind(demo, info[matched, ])

sel_col <- c("V3", "V4", "V5", "V9", "title", "characteristics_ch1.3")
demo_final <- demo_add[,sel_col]
save(demo_final, file = "GSE147507_stage1/demoInfo_GSE147507.RData")

### End ###
