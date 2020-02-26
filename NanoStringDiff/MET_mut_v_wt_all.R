###install NanoStringDiff
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("NanoStringDiff")

###specify directory and path to counts file
directory <- "C:/Users/jasmi/Documents/Bild Lab/NanoString/Salgia Project/Salgia MET Project/R Code/NanoString diff raw data all MET MUT and WT/"
path <- paste(directory,"MET_mut_vs_wt_expression.csv", sep="/")

###load the phenotype data
designs=read.csv(file="C:/Users/jasmi/Documents/Bild Lab/NanoString/Salgia Project/Salgia MET Project/R Code/NanoString diff raw data all MET MUT and WT/MET_mut_vs_wt_phenotype_labels.csv", header=TRUE, row.names=1, sep=",")
designs

###Call NanoStringDiff and create the NS set
library("NanoStringDiff")
NanoStringData=createNanoStringSetFromCsv(path,header=TRUE,designs)

###Make a matrix to prepare for normalization, that only considers one phenotype label, MET status
pheno=pData(NanoStringData)
group=pheno$group
design.full=model.matrix(~0+group)
design.full

###Normalization of the data
NanoStringData=estNormalizationFactors(NanoStringData)
positiveFactor(NanoStringData)
negativeFactor(NanoStringData)
housekeepingFactor(NanoStringData)

###perform glm ratio test of MET MUT to control WT
###control group gets contrast -1 designation, group of interest gets +1
result=glm.LRT(NanoStringData,design.full,contrast=c(1,-1))
head(result$table)
str(result)

###write result out to csv
write.csv(result[["table"]], file = "MET_MUT_vs_WT_all_nanostringdiff.csv")

