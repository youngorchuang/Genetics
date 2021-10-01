load("GSE147507_stage1/Allelic5SNPinOtherVirus_GSE147507.RData")

rsid <- c("rs1050240", "rs41541216", "rs41545916", "rs41541515", "rs41551014")

index_id <- grep("Series9.NHBE", dataMerge[,3])
index_cntl <- grep("Mock", dataMerge[,4])
index_ifn <- grep("IFNB", dataMerge[,4])

index_1 <- intersect(index_id, index_ifn)
index_2 <- intersect(index_id, index_cntl)
index_over <- union(index_1, index_2)

dataSer9 <- dataMerge[index_over, ]
dataSer9$time <- "0h"

index_ifn <- grep("IFNB", dataSer9[,4])
dataSer9$time[index_ifn] <- do.call(rbind, strsplit(as.character(dataSer9[index_ifn,3]), split="IFNB"))[,2]
dataSer9$time <- as.factor(gsub("h", "", dataSer9$time))

dataSer9_mean <- as.data.frame(tapply(dataSer9[,2], list(droplevels(dataSer9[,5]), droplevels(dataSer9[,1])), mean))
calSD <- as.data.frame(tapply(dataSer9[,2], list(droplevels(dataSer9[,5]), droplevels(dataSer9[,1])), sd))

dataStacked <- stack(dataSer9_mean)
dataStackSD <- stack(calSD)
dataStackSD[is.na(dataStackSD[,1]), 1] <- 0   # one value only having one replicate

num <- dim(dataSer9_mean)[2]    # replication based on the column size
time_rep <- rep(as.numeric(rownames(dataSer9_mean)), num)

dataPlot = data.frame(snp = as.character(dataStacked[,2]), ratio = as.numeric(dataStacked[,1]), lower = as.numeric(dataStacked[,1] - as.numeric(dataStackSD[,1])), upper = as.numeric(dataStacked[,1] + as.numeric(dataStackSD[,1])), series = as.numeric(time_rep), fc = 1,stringsAsFactors=FALSE)


# Mean fold change relative to the 0h of mock treatment
snpUniq <- unique(dataPlot[,1])
for(id in snpUniq){
	matched <- grep(id, dataPlot[,1])
	
	index_mock <- which(dataPlot[matched, "series"] %in% "0")
	for (x in matched){
		dataPlot[x,"fc"] <- dataPlot[x, "ratio"]/dataPlot[matched[index_mock],"ratio"]
	}
}

library("ggplot2")

ggplot(dataPlot, aes(x=as.factor(series))) + geom_bar(aes(y=ratio, fill=series), stat = "identity", position = "dodge", width = 0.8) + geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge", width = 0.25) + geom_point(data=dataPlot, aes(y=(fc + 15), group = 1), size=2, shape=21) + geom_line(data=dataPlot, aes(y=(fc + 15), group = 1, color=I(4)), alpha=0.7, size=1) + geom_hline(yintercept=16, color='grey', linetype="dashed") + facet_wrap(~ snp , scales = "free_x", ncol=1) + xlab("Duration times (h) of IFNb treatment") + ylab("Allelic FC\n") + scale_y_continuous(name = "Allelic FC\nRef vs Alt alleles", sec.axis = sec_axis(~ ., name = "Fold change\nrelative to 0h", breaks = c(15, 17, 19), labels = seq(0, 4, 2)))  + theme(legend.position = "none", panel.background = element_blank(), text = element_text(size=11), axis.line = element_line(colour="black"), axis.text.x = element_text(size = 10), 
   axis.title.y.right = element_text(color = 4, size=12)
)  # axis.line.y.right = element_line(colour=4)
## Block ##
