#Jenny Smith
#January 5, 2017
#Original Script Author: Hamid Bolouri

#purpose: convert RPKM rna-seq to TPM (transcript per million) 


getwd()

#set working directory
# setwd("H:/RNAseq_BCCA28Apr2016/gene_coverage_GSC-1367/")
setwd("H:/RNA_seq_Analysis/2017.02.23_LAML_TCGA_MMRN1/mRNAseq/BCCA_Illumina_data")

LAML <- read.table("geneLevelExpn/RPKM_TCGA_LAML.gene.expression.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


head(LAML[,1:20])
dim(LAML) #20443  


#create an object for all the files to be converted.
allFiles <- dir(pattern = ".transcript.normalized")					

#number of files
length(allFiles)
#494 files

#ensure the list of files is loaded correctly
head(allFiles)

#############################################################
#Convert 1 file at first.
oneFile <- read.table(file = "A37371_1_lanes.stranded_collapsed_coverage.transcript.normalized", header = FALSE, sep = "\t")

#check that the file loaded correctly.
head(oneFile)

dim(oneFile)
#35190 rows, 18 columns
#column names can be found in gene_coverage_analysis_filetypes_readme.pdf
#column 15 is the RPKM gene transcripts, column 16 is the gene symbol
#column 9 is the raw counts, but they are fractional values - not intergers. Cannot be directly used in differential expression analysis. 


#variable to hold the RPKM counts
oneFileRPKM = oneFileID = NULL

#ID is the numeric identifier from BCCA.
oneFileID <- gsub("(A[0-9]+)_1_.+", "\\1", allFiles[1])
oneFileRPKM <- as.data.frame(oneFile[,15])
oneFilegeneSymbol <- oneFile[,16]
colnames(oneFileRPKM) <- oneFileID
oneFileRPKM <- cbind(oneFilegeneSymbol, oneFileRPKM)

# Based on https://groups.google.com/forum/#!topic/rsem-users/W9RQrZIOzA4
#conversion factor is the sum of RPKM for all genes for one patient divided by
#one million.
convertFactor <- sum(oneFileRPKM[,2])/ 1E6

#TPM is the RPKM value for each gene divided by the conversion factor.
oneFileTPM <- oneFileRPKM[,2] / convertFactor

head(oneFileTPM)
####################################################





# i=0
# allRPKM = IDs = allFracCounts = NULL
# for (F in allFiles) {
#   aFile <- read.delim(F, sep="\t", header=FALSE, as.is=TRUE)
# 	ID <- gsub("(A[0-9]+).+", "\\1", F)
# 	IDs <- c(IDs, ID)
# 	RPKM <- aFile[,15]
# 	fractionalCounts <- aFile[,9]#NB: These are fractional counts (not interger raw counts)
# 	allRPKM <- cbind(allRPKM, RPKM)
# 	allFracCounts <- cbind(allFracCounts, fractionalCounts)
# 
# 	i <- i + 1
# 	print(i)
# }



allRPKM <- as.data.frame(sapply(LAML[which(LAML[1,] == "RPKM")][-1,], as.double))
allMLN <- as.data.frame(sapply(LAML[which(LAML[1,] == "median_length_normalized")][-1,], as.double))
allRawcounts <- as.data.frame(sapply(LAML[which(LAML[1,] == "raw_counts")][-1,], as.double))

dim(allRPKM) #180

head(LAML[,1:5])
head(allRPKM[,1:5])

head(allMLN[,1:5])
head(allRawcounts[,1:5])




#make the column names the BCCA identifiers. 
# colnames(allRPKM) <- IDs
# colnames(allRawCounts) <- IDs
# colnames(allMLN) <- IDs
# 
# head(allRPKM)
# head(allFracCounts)

# Based on https://groups.google.com/forum/#!topic/rsem-users/W9RQrZIOzA4
conversionFactor <- apply(allRPKM, 2, sum) / 1E6


head(conversionFactor)
tail(conversionFactor)


#Convert the allRPKM to the TPM normalized values. 
allTPM <- allRPKM
for (i in 1:dim(allRPKM)[2]) {
  allTPM[ ,i] <- allRPKM[ ,i] / conversionFactor[i]
  print(i)
}

head(allRPKM[,1:5])
head(allTPM[,1:5])



#Convert the BCCA sample IDs to the patient identifiers (USI)
#Since the BCCA sample IDs and USIs are in the same order, 
#just replace the BCCA sample IDs with the patient USIs for
#the column names. 
# idMap <- read.table(file = "ID.map.txt", sep = "\t", header = TRUE)
# 
# head(idMap)
# head(allRPKM)
# head(allTPM)
# 
# colnames(allRPKM) <- idMap$Patient
# colnames(allTPM) <- idMap$Patient
# colnames(allFracCounts) <- idMap$Patient
# 
# head(allRPKM)
# head(allTPM)
# head(allFracCounts)
# 
# dim(allRPKM)
# dim(allTPM)
# dim(allFracCounts)


#gene symbols for one of the the RNA-seq datasets will be used for the 
#gene symbol for the RPKM, TPM, and rawCounts files since they are in the same order in
#all of the RNA-seq datasets. 

# geneSymbol <- aFile[,16]
# 
# head(aFile)
# head(geneSymbol)
# length(geneSymbol)

geneSymbol <- LAML$Hybridization.REF[-1]
head(geneSymbol) #20442


sapply(allTPM, class)

colSums(allTPM[,2:ncol(allTPM)])







#add the gene symbols to the data frames. 
allRPKM <- data.frame(geneSymbol, allRPKM)
allTPM <- data.frame(geneSymbol, allTPM)
allRawcounts <- data.frame(geneSymbol, allRawcounts)
allMLN <- data.frame(geneSymbol, allMLN)


head(allRPKM[,1:10])
head(allTPM[,1:10])
head(allRawcounts[,1:10])
head(allMLN[,1:10])


dim(allRPKM)
dim(allTPM)
dim(allRawcounts)
dim(allMLN)


#write the files to the csv file. 
write.csv(allRPKM, file="LAML_TCGA_RPKM_geneExpression.csv", row.names=TRUE, 
          col.names=TRUE, quote=FALSE)

write.csv(allTPM, file="LAML_TCGA_TPM_geneExpression.csv", row.names=TRUE, 
          col.names=TRUE, quote=FALSE)

write.csv(allRawcounts, file="LAML_TCGA_RawCounts_geneExpression.csv", row.names=TRUE, 
          col.names=TRUE, quote=FALSE)

write.csv(allMLN, file = "LAML_TCGA_MedianLengthNormalized_geneExpression.csv", row.names = FALSE, quote = FALSE)


save(allRPKM, file="LAML_TCGA_RPKM.RData")
save(allTPM, file="LAML_TCGA_TPM.RData")
save(allRawCounts, file="LAML_TCGA_RawCounts.RData")




####################################################
#Jenny Smith
#Purpose: Practice creating a function pullout RPKMs from multiple files and convert to TPM.

#function to pullout the RPKMs from multiple files.  
RPKM_List <- function(x){
  RPKM = ID = NULL
  aFile <- read.delim(F, sep="\t", header=FALSE, as.is=TRUE)
  #ID <- gsub("(A[0-9]+)_1_.+", "\\1", x)
  #geneSymbol <- aFile[,16]
  RPKM <- aFile[,15]
  #allRPKM <- cbind(allRPKM, RPKM)
  #IDs <- c(IDs, ID)
  #colnames(RPKM) <- ID #use rbind instead here?
  #RPKM <- rbind(ID, RPKM)
  return(RPKM)
}

#Function for the TPM conversion. 
RPKM_to_TPM <- function(x){
  conversionFact <- sum(x) / 1E6
  TPM <- x / conversionFact
  return(TPM)
}

#make a character vector with the BCCA sample IDs
IDs <- gsub("(A[0-9]+).+", "\\1", allFiles)

head(IDs)

#Call the function to select the RPKMs for each file 
test1 <- RPKM_List(allFiles[1])

head(test1)
class(test1)

#Use sapply to call the function on many files. 
test2 <- sapply(allFiles[1:10], RPKM_List)

head(test2)
dim(test2)

#Add the column names to the the dataframe
colnames(test2) <- IDs[1:10]
head(test2)

#call the TPM function to convert 
test3 <- RPKM_to_TPM(test2[,1])

#use the sapply to make a 
test4 <- sapply(test2, RPKM_to_TPM)

#######################################################

