#GSE (Series) data set (Contains GSM and GPL) 
getGSEset <- function(GSEset) {
	gse <- getGEO(GSEset,GSEMatrix = F)
	#Meta data 
	print(head(Meta(gse)))

	#List GSM (Sample record)
	print(names(GSMList(gse)))
	Table(GSMList(gse)[[1]])

	#List GPL (annotation/platform)
	print(names(GPLList(gse)))
	Table(GPLList(gse)[[1]])

	#Create dataframe with just the IDs
	IDs <- Table(GSMList(gse)[[1]])
	IDs$VALUE <- NULL

	for (gsm in GSMList(gse)) {
	  temp <- Table(gsm)
	  colnames(temp)[2] <- gsm@header$geo_accession
	  IDs <- merge(IDs, temp, by = "ID_REF",all=T)
	  rm(temp)
	}

	write.csv(IDs, file="GSE_sample_RMA_signal_log2.csv")
	
	phenoData <- pData(phenoData(gse[[1]]))
	write.csv(phenoData, file="GSE_phenoData.csv")
	
	featData <- fData(gse[[1]])
	write.csv(featData, file = "GSE_featData.csv")
}


# #Changing to expression Data set
# gse <- getGEO("GSE39645")
# show(gse)
# #Get phenotype Data
# phenoData <- pData(phenoData(gse[[1]]))
# write.csv(phenoData, file="GSE_phenoData.csv")
# #Get feature Data
# #Inconsistent (fData(featureData(gse[[1]]))) doesn't work
# featData <- fData(gse[[1]])
# write.csv(featData, file = "GSE_featData.csv")

