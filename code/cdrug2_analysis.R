options(stringsAsFactors=FALSE)

saveres <- file.path("output")

source(file.path("code", "cdrug2_foo.R"))

### install all the libraries at once
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("VennDiagram", "Hmisc", "xtable", "RColorBrewer", "pROC", "Biobase", "genefu", "PharmacoGx", "xlsx"))
library(VennDiagram)
library(Hmisc)
library(xtable)
library(RColorBrewer)
library(pROC)
library(Biobase)
library(genefu)
library(xlsx)
library(psych)
# install the latest devel version of the PharmacoGx package
# library(devtools)
# devtools::install_github("bhklab/PharmacoGx", ref="master")
library(PharmacoGx)


#################################################
## get pharmacogenomic datasets
#################################################

### download curated pharmacogenomic data for GDSC and CCLE
GDSC <- PharmacoGx::downloadPSet("GDSC", saveDir=file.path(saveres, "PSets"))
CCLE <- PharmacoGx::downloadPSet("CCLE", saveDir=file.path(saveres, "PSets")) 

common <- intersectPSet(pSets = list("CCLE"=CCLE, "GDSC"=GDSC), intersectOn = c("cell.lines", "drugs"), strictIntersect = TRUE)
drugs <- c("paclitaxel", "17-AAG", "PD-0325901", "AZD6244", "TAE684", "AZD0530", "PD-0332991", "Crizotinib", "PLX4720", "Nutlin-3", "lapatinib", "Nilotinib", "PHA-665752", "Erlotinib", "Sorafenib")

##check if snp fingure printing outliers are in the common cells screened for common drugs between studies
snp.outliers <- c(
  "LC-1F",
  "HCC1937",
  "MDA-MB-468",
  "HuH-7",
  "SW403",
  "COR-L51",
  "MOG-G-CCM",
  "NB4")

require(VennDiagram) || stop("Library gdata is not available!")
mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")

### venn diagram of common cell lines between studies
pdf(file.path(saveres, "celllines.pdf"), height=4, width=4)
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@cell), area2=nrow(GDSC@cell), cross.area=nrow(common$CCLE@cell), fill=c(mycol[1], mycol[2]), lty="blank",cex=1.5, cat.cex=1, cat.col = c("black", "black"))
dev.off()

### venn diagram of common drugs between studies
pdf(file.path(saveres, "drugs.pdf"), height=4, width=4)
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@drug), area2=nrow(GDSC@drug), cross.area=nrow(common$CCLE@drug), fill=c(mycol[1], mycol[2]), lty="blank",cex=1.5, cat.cex=1, cat.col = c("black", "black"))
dev.off()                 


##venn diagram of common tissues between studies
#mycol <- RColorBrewer::brewer.pal(n=7, name="Set3")
mycol <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f",
           "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99")#,"#b15928")
pdf(file = file.path(saveres, "tissues.pdf"), height=8, width=10)
temp <- c(table(common$GDSC@cell[,"tissueid"])[15], table(common$GDSC@cell[,"tissueid"])[1:14],   table(common$GDSC@cell[,"tissueid"])[16] ,table(common$GDSC@cell[,"tissueid"])[18:23], table(common$GDSC@cell[,"tissueid"])[17])               
pie(temp, labels = gsub("_", " ", capitalize(names(temp))), col=mycol, main="", radius=1, cex=0.8)#,theta=pi/4, labelcex = .8, labelrad = 1.1)
#pie(Isoforms_No, lbls, col=rainbow(length(lbls)), main = "Distribution of Genes based on the number of isoforms")
dev.off()


### plotting noisy curves, curves might be noisy in one or both studies(Supplementary File 2)
ccle.filter <- PharmacoGx:::filterNoisyCurves(common$CCLE, nthread=detectCores())
lapply(ccle.filter, length)

gdsc.filter <- PharmacoGx:::filterNoisyCurves(common$GDSC, nthread=detectCores())
lapply(gdsc.filter, length)


pdf(file.path(saveres, "Supplementary_File_2.pdf"), height=15.5, width=12.5)
par(mfrow = c(4,3))
exps <- NULL
for (exp in ccle.filter$noisy)
{
  cell <- sensitivityInfo(common$CCLE)[exp, "cellid"]
  drug <- sensitivityInfo(common$CCLE)[exp, "drugid"]
  exps <- c(exps, paste(cell, drug, sep="_"))
  drugDoseResponseCurve(pSets=common, drug=drug, cellline=cell)
  legend("bottomleft", legend="CCLE", bty="n")
}

for (exp in gdsc.filter$noisy) {
  cell <- sensitivityInfo(common$GDSC)[exp, "cellid"]
  drug <- sensitivityInfo(common$GDSC)[exp, "drugid"]
  if (!(paste(cell, drug, sep="_") %in% exps)) {
    drugDoseResponseCurve(pSets=common, drug=drug, cellline=cell)
  }
  legend("bottomleft", legend="GDSC", bty="n")
  
}
dev.off()
### plotting GDSC noisy curves, curves might be noisy in one or both studies


pdf(file.path(saveres, "GDSC_noisy.pdf"), height=15.5, width=12.5)
par(mfrow = c(4,3))
for (exp in gdsc.filter$noisy) {
  cell <- sensitivityInfo(GDSC)[exp, "cellid"]
  drug <- sensitivityInfo(GDSC)[exp, "drugid"]
  drugDoseResponseCurve(pSets=GDSC, drug=drug, cellline=cell)
}
dev.off()


####resistance/sensitive cases
cases <- rbind(
  c("CAL-85-1", "17-AAG"),
  c("HT-29", "PLX4720"),
  c("COLO-320-HSR", "AZD6244"),
  c("HT-1080", "PD-0332991"))

for (i in 1:nrow(cases)) {
  pdf(file.path(saveres, sprintf("%s_%s.pdf", cases[i,2], cases[i,1])), height=5, width=5)
  drugDoseResponseCurve(pSets=common, drug=cases[i,2], cellline=cases[i,1], legends.label="auc_published", plot.type="Fitted", ylim=c(0,130))
  dev.off()
}

###filtered out cases
cases <- rbind(
  c("KNS-62", "17-AAG"),
  c("HCC70", "PD-0332991"),
  c("EFM-19", "PD-0325901"),
  c("LS-513", "Nutlin-3"))

for (i in 1:nrow(cases)) {
  pdf(file.path(saveres, sprintf("%s_%s.pdf", cases[i,2], cases[i,1])), height=5, width=5)
  drugDoseResponseCurve(pSets=common, drug=cases[i,2], cellline=cases[i,1])
  dev.off()
}


###Internal concordance in GDSC
cases <- rbind(
  c("A498", "AZD6482"),
  c("BPH-1", "AZD6482"),
  c("KURAMOCHI","AZD6482"),
  c("NCI-H1092","AZD6482"))
for (i in 1:nrow(cases)) {
  pdf(file.path(saveres, sprintf("%s_%s.pdf", cases[i,2], cases[i,1])), height=5, width=5)
  drugDoseResponseCurve(pSets=GDSC, drug=cases[i,2], cellline=cases[i,1], summarize.replicates = FALSE, plot.type = "Fitted", mycol=RColorBrewer::brewer.pal(n=4, name="Set1")[c(3,4)])
  exp <- which(sensitivityInfo(GDSC)$cellid==cases[i,1] & sensitivityInfo(GDSC)$drugid==cases[i,2])
  ABC <- computeABC(conc1 = GDSC@sensitivity$raw[exp[1],,"Dose"], conc2 = GDSC@sensitivity$raw[exp[2],,"Dose"], viability1 = GDSC@sensitivity$raw[exp[1],,"Viability"], viability2 = GDSC@sensitivity$raw[exp[2],,"Viability"])
  legend("bottomright", legend = sprintf("ABC= %s", round(ABC, digits=2)), bty="n")
  dev.off()
}

############remove noisy experiments and mislabeled cell lines (snp.outliers)
myf <- file.path(saveres, "PSets", "common_clarified.RData")
if(!file.exists(myf)){
  common$CCLE@sensitivity$info <- common$CCLE@sensitivity$info[ccle.filter$ok, ]
  common$CCLE@sensitivity$raw <- common$CCLE@sensitivity$raw[ccle.filter$ok, , ]
  common$CCLE@sensitivity$profiles <- common$CCLE@sensitivity$profiles[ccle.filter$ok, ]
  
  common$GDSC@sensitivity$info <- common$GDSC@sensitivity$info[gdsc.filter$ok, ]
  common$GDSC@sensitivity$raw <- common$GDSC@sensitivity$raw[gdsc.filter$ok, , ]
  common$GDSC@sensitivity$profiles <- common$GDSC@sensitivity$profiles[gdsc.filter$ok, ]
  
  cells <- intersectList(phenoInfo(common$CCLE, "rna")$cellid, phenoInfo(common$GDSC, "rna2")$cellid, unique(sensitivityInfo(common$CCLE)$cellid), unique(sensitivityInfo(common$GDSC)$cellid))
  cells <- setdiff(cells, snp.outliers)
  common.clarified <- intersectPSet(pSets = list("CCLE"=common$CCLE, "GDSC"=common$GDSC), intersectOn = c("cell.lines", "drugs"), cells=cells)
  save(common.clarified, file=myf)
} else{
  load(myf)
}
 
myf <- file.path(saveres, "signatures_data.RData")
if(!file.exists(myf)){
  drugs <- c("paclitaxel", "17-AAG", "PD-0325901", "AZD6244", "TAE684", "AZD0530", "PD-0332991", "Crizotinib", "PLX4720", "Nutlin-3", "lapatinib", "Nilotinib", "PHA-665752", "Erlotinib", "Sorafenib")
  features <- intersect(rownames(featureInfo(CCLE, "rna")), rownames(featureInfo(GDSC, "rna2")))
  save(drugs, features, common.clarified, file=myf)  
} else{
  load(myf)
}

### plotting all the common experiments between studies(Supplementary File 3)
pdf(file.path(saveres, "Supplementary_File_3.pdf"), height=15.5, width=12.5)
par(mfrow = c(4,3))
for (exp in rownames(sensitivityInfo(common.clarified$CCLE))) {
  cell <- sensitivityInfo(common.clarified$CCLE)[exp, "cellid"]
  drug <- sensitivityInfo(common.clarified$CCLE)[exp, "drugid"]
  drugDoseResponseCurve(pSets=common.clarified, drug=drug, cellline=cell, plot.type="Fitted", ylim=c(0,130))
}
dev.off() 


############remove noisy experiments and mislabeled cell lines (snp.outliers) from actual pSets
myf <- file.path(saveres, "PSets", "GDSC_clarified.RData")
if(!file.exists(myf)) {
  gdsc.filter <- PharmacoGx:::filterNoisyCurves(GDSC, nthread=detectCores())
  lapply(gdsc.filter, length)
  #$noisy
  #[1] 2315
  
  #$ok
  #[1] 77588
  
  GDSC.clarified <- GDSC
  GDSC.clarified@sensitivity$info <- GDSC.clarified@sensitivity$info[gdsc.filter$ok, ]
  GDSC.clarified@sensitivity$raw <- GDSC.clarified@sensitivity$raw[gdsc.filter$ok, , ]
  GDSC.clarified@sensitivity$profiles <- GDSC.clarified@sensitivity$profiles[gdsc.filter$ok, ]
  save(GDSC.clarified, file=myf)
} else {
  load(myf)
}

myf <- file.path(saveres, "PSets", "CCLE_clarified.RData")
if(!file.exists(myf)) {
  
  ccle.filter <- PharmacoGx:::filterNoisyCurves(CCLE, nthread=detectCores())
  lapply(ccle.filter, length)
  # $noisy
  # [1] 123
  # 
  # $ok
  # [1] 11547
  
  CCLE.clarified <- CCLE
  CCLE.clarified@sensitivity$info <- CCLE.clarified@sensitivity$info[ccle.filter$ok, ]
  CCLE.clarified@sensitivity$raw <- CCLE.clarified@sensitivity$raw[ccle.filter$ok, , ]
  CCLE.clarified@sensitivity$profiles <- CCLE.clarified@sensitivity$profiles[ccle.filter$ok, ]
  save(CCLE.clarified, file=myf)
}else {
  load(myf)
}

### plotting all internal replicates in GDSC 

replicates <- sensitivityInfo(GDSC.clarified)[which(sensitivityInfo(GDSC.clarified)$drugid == "AZD6482"),]
replicates <- replicates[which(duplicated(replicates$cellid)),] #577 it was 617 before removing the noisy ones

ABC <- NULL
for (exps in 1:nrow(replicates)) {
  cell <- replicates[exps, "cellid"]
  drug <- replicates[exps, "drugid"]
  exp <-   which(sensitivityInfo(GDSC.clarified)$cellid==cell & sensitivityInfo(GDSC.clarified)$drugid==drug)
  ABC <- c(ABC, computeABC(conc1 = GDSC.clarified@sensitivity$raw[exp[1],,"Dose"], conc2 = GDSC.clarified@sensitivity$raw[exp[2],,"Dose"], viability1 = GDSC.clarified@sensitivity$raw[exp[1],,"Viability"], viability2 = GDSC.clarified@sensitivity$raw[exp[2],,"Viability"]))
}
biologucal_replicates_ABC <- ABC[which(!is.na(ABC))]
pdf(file.path(saveres, "Supplementary_File_4.pdf"), height=15.5, width=12.5)
par(mfrow = c(4,3))

for (i in order(ABC, decreasing = T)) {
  exp <- rownames(replicates)[i]
  cell <- replicates[exp, "cellid"]
  drug <- replicates[exp, "drugid"]
  drugDoseResponseCurve(pSets=GDSC, drug=drug, cellline=cell, summarize.replicates = FALSE, plot.type = "Fitted", mycol=RColorBrewer::brewer.pal(n=4, name="Set1")[c(3,4)])
  legend("bottomleft", legend = sprintf("ABC= %s", round(ABC[i], digits=2)), bty="n")
}

dev.off()


##boxplot of ABC across celllines and across drugs

ABC <- NULL
for (exp.ccle in rownames(sensitivityInfo(common.clarified$CCLE))) {
  cell <- sensitivityInfo(common.clarified$CCLE)[exp.ccle, "cellid"]
  drug <- sensitivityInfo(common.clarified$CCLE)[exp.ccle, "drugid"]
  exp.gdsc <-   rownames(sensitivityInfo(common.clarified$GDSC))[which(sensitivityInfo(common.clarified$GDSC)$cellid==cell & sensitivityInfo(common.clarified$GDSC)$drugid==drug)]
  ABC <- c(ABC, computeABC(conc1 = common.clarified$CCLE@sensitivity$raw[exp.ccle,,"Dose"], 
                           conc2 = common.clarified$GDSC@sensitivity$raw[exp.gdsc,,"Dose"], 
                           viability1 = common.clarified$CCLE@sensitivity$raw[exp.ccle,,"Viability"], 
                           viability2 = common.clarified$GDSC@sensitivity$raw[exp.gdsc,,"Viability"]))
}
names(ABC) <- rownames(sensitivityInfo(common.clarified$CCLE))
ABC <- cbind(ABC, "cellid"=sensitivityInfo(common.clarified$CCLE)[names(ABC),"cellid"], "drugid"=sensitivityInfo(common.clarified$CCLE)[names(ABC),"drugid"])

Inter_studies_ABC <- ABC[which(!is.na(ABC[,"ABC"])),"ABC"]
###wilcoxon test
biologucal_replicates_ABC <- as.numeric(biologucal_replicates_ABC)
Inter_studies_ABC <- as.numeric(Inter_studies_ABC)
tt <- matrix(NA, ncol=2, nrow=max(length(Inter_studies_ABC), length(biologucal_replicates_ABC)))
colnames(tt) <- c("Intra", "Inter")
tt[1:length(biologucal_replicates_ABC),"Intra"] <- biologucal_replicates_ABC
tt[1:length(Inter_studies_ABC),"Inter"] <- Inter_studies_ABC

pdf(file.path(saveres, "inter_intra_abc.pdf"), height=7, width=7)
par(mar=c(9,5,5,2))
boxplot(tt, las = 2, col = "gray", cex.lab=1, cex.axis=1, pch=19, ylab="ABC", outpch=20, outcex=0.5)
dev.off()
test <- wilcox.test(Inter_studies_ABC, biologucal_replicates_ABC, paired=FALSE)
test$p.value
###box plot for ABC
tt <- as.data.frame(matrix(NA, ncol=length(drugNames(common.clarified$CCLE)), nrow=length(cellNames(common.clarified$CCLE))))
colnames(tt)=drugNames(common.clarified$CCLE)
rownames(tt)=cellNames(common.clarified$CCLE)
for (ii in 1:nrow(ABC)) {
  tt[ABC[ii, "cellid"], ABC[ii, "drugid"]] <- as.numeric(ABC[ii, "ABC"])
}
emp.cell <- NULL
for(i in 1:nrow(tt)) {
  if(all(is.na(tt[i,]))) {
    emp.cell <- c(emp.cell,i)
  }
}
if(!is.null(emp.cell)){
  tt <- tt[-emp.cell,]
}
##drugs
pdf(file.path(saveres, "drugs_abc.pdf"), height=5, width=10)
par(mar=c(9,5,5,2))
oo <- order(apply(tt, 2, median, na.rm=TRUE), decreasing=TRUE)
boxplot(tt[,oo], las = 2, col = "gray", cex.lab=1, cex.axis=1, pch=19, ylab="ABC", outpch=20, outcex=0.5)
dev.off()

##cells
pdf(file.path(saveres, "cells_ABC.pdf"), height=5, width=20)
par(mar=c(9,5,5,2))
oo <- order(apply(tt, 1, median, na.rm=TRUE), decreasing=TRUE)
boxplot(t(tt[oo,]), las = 2, col = "gray", cex.lab=.5, cex.axis=.5, pch='.', ylab="ABC", outpch=20, outcex=0.5)
dev.off()

##histogram
ttt <- do.call(c,tt)
pdf(file.path(saveres, "ABC_hist.pdf"),width=5,height=5)
hist(ttt,breaks=50, main="",xlab="ABC")
dev.off()


## since all the kept experiments are common experiments it's just needed to check for one of the studies
unique(sensitivityInfo(common.clarified$GDSC)[which(sensitivityInfo(common.clarified$GDSC)$cellid %in% snp.outliers),"cellid"])
sort(unique(sensitivityInfo(common.clarified$GDSC)[which(sensitivityInfo(common.clarified$GDSC)$cellid %in% snp.outliers),"drugid"]))


###biomarkers consistancy

load(file.path(saveres, "PSets","signatures_data.RData"))
load(file.path(saveres, "PSets", "GDSC_clarified.RData"))
load(file.path(saveres, "PSets", "CCLE_clarified.RData"))

myf <- file.path(saveres, "PSets", "Sigs", "common_ccle_sig_rna.RData")
if(!file.exists(myf))
{
  common.ccle.sig.rna <- drugSensitivitySig(pSet=common.clarified$CCLE, mDataType="rna", drugs=drugs, features=features,sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(common.ccle.sig.rna, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "common_gdsc_sig_rna2.RData")
if(!file.exists(myf))
{
  common.gdsc.sig.rna2 <- drugSensitivitySig(pSet=common.clarified$GDSC, mDataType="rna2", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(common.gdsc.sig.rna2, file=myf)
  
}else{
  load(myf)
}

continuous.common <- biomarkers.validation(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=continuous.common, method="continuous", cell="common")
continuous.common.100 <- biomarkers.validation(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous", drugs=drugs, fdr.cut.off=0.05, nperm=100, top.ranked=100)
plotValidation(validation.result=continuous.common.100, method="continuous_top100", cell="common")


myf <- file.path(saveres, "PSets", "Sigs", "common_gdsc_sig_rna2_binary.RData")
if(!file.exists(myf))
{
  common.gdsc.sig.rna2.bin <- drugSensitivitySig(pSet=common.clarified$GDSC, mDataType="rna2", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", sensitivity.cutoff=0.2, nthread=4)
  save(common.gdsc.sig.rna2.bin, file=myf)
  
}else{
  load(myf)
}
myf <- file.path(saveres, "PSets", "Sigs", "common_ccle_sig_rna_binary.RData")
if(!file.exists(myf))
{
  common.ccle.sig.rna.bin <- drugSensitivitySig(pSet=common.clarified$CCLE, mDataType="rna", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", sensitivity.cutoff=0.2, nthread=4)
  save(common.ccle.sig.rna.bin, file=myf)
}else{
  load(myf)
}



binary.commom <- biomarkers.validation(ccle.sig.rna=common.ccle.sig.rna.bin, gdsc.sig.rna=common.gdsc.sig.rna2.bin, cell="common", method="binary", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=binary.commom, method="binary", cell="common")
binary.commom.100 <- biomarkers.validation(ccle.sig.rna=common.ccle.sig.rna.bin, gdsc.sig.rna=common.gdsc.sig.rna2.bin, cell="common", method="binary", drugs=drugs, fdr.cut.off=0.05, nperm=100, top.ranked=100)
plotValidation(validation.result=binary.commom.100, method="binary_top100", cell="common")


test <- wilcox.test(continuous.common.100$validation, binary.commom.100$validation, alternative="greater", paired=TRUE)
test$p.value


load(file.path(saveres, "PSets", "GDSC_clarified.RData"))
load(file.path(saveres, "PSets", "CCLE_clarified.RData"))
load(file.path(saveres, "signatures_data.RData"))

myf <- file.path(saveres, "PSets", "Sigs", "ccle_sig_rna.RData")
if(!file.exists(myf))
{
  ccle.sig.rna <- drugSensitivitySig(pSet=CCLE.clarified, mDataType="rna", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(ccle.sig.rna, file=myf)
}else{
  load(myf)
}


myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_rna2.RData")
if(!file.exists(myf))
{
  gdsc.sig.rna2 <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="rna2", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  
  save(gdsc.sig.rna2, file=myf)
}else{
  load(myf)
}

cnv.fetures <- intersect(rownames(fData(CCLE@molecularProfiles$cnv)), rownames(fData(GDSC@molecularProfiles$cnv)))

myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_cnv.RData")
if(!file.exists(myf))
{
  
  gdsc.sig.cnv <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="cnv", drugs=drugs, features=cnv.fetures, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(gdsc.sig.cnv, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "ccle_sig_cnv.RData")
if(!file.exists(myf))
{
  ccle.sig.cnv <- drugSensitivitySig(pSet=CCLE.clarified, mDataType="cnv", drugs=drugs, features=cnv.fetures, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(ccle.sig.cnv, file=myf)
} else{
  load(myf)
}
continuous.all <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=continuous.all, method="continuous", cell="all")
continuous.all.100 <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous", drugs=drugs, fdr.cut.off=0.05, nperm=100, top.ranked=100)
plotValidation(validation.result=continuous.all.100, method="continuous_top100", cell="all")



myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_rna2_binary.RData")
if(!file.exists(myf))
{
  gdsc.sig.rna2.bin <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="rna2", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  
  save(gdsc.sig.rna2.bin, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "ccle_sig_rna_binary.RData")
if(!file.exists(myf))
{
  ccle.sig.rna.bin <- drugSensitivitySig(pSet=CCLE.clarified, mDataType="rna", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", sensitivity.cutoff=0.2, nthread=4)
  save(ccle.sig.rna.bin, file=myf)
}else{
  load(myf)
}

binary.all <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna.bin, gdsc.sig.rna=gdsc.sig.rna2.bin, cell="all", method="binary", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=binary.all, method="binary", cell="all")
binary.all.100 <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna.bin, gdsc.sig.rna=gdsc.sig.rna2.bin, cell="all", method="binary", drugs=drugs, fdr.cut.off=0.05, nperm=100, top.ranked=100)
plotValidation(validation.result=binary.all.100, method="binary_top100", cell="all")


test <- wilcox.test(continuous.all.100$validation, binary.all.100$validation, alternative="greater", paired=TRUE)
test$p.value


test <- wilcox.test(continuous.all.100$validation, continuous.common.100$validation, alternative="greater", paired=TRUE)
test$p.value

myf <- file.path(saveres, "PSets", "Sigs", "gdsc_mutation.RData")
if(!file.exists(myf))
{
  gdsc.sig.mutation <- drugSensitivitySig(pSet=GDSC, mDataType="mutation", sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4)
  
  save(gdsc.sig.rna2.bin, file=myf)
}else{
  load(myf)
}

features.mutation.gdsc <- sapply(rownames(molecularProfiles(GDSC.clarified, "mutation")), function(x){if(!all(is.na(molecularProfiles(GDSC.clarified, "mutation")[x,]))){x}},simplify=T)
features.mutation.gdsc <- do.call(c, features.mutation.gdsc)
features.mutation.ccle <- sapply(rownames(molecularProfiles(CCLE.clarified, "mutation")), function(x){if(!all(is.na(molecularProfiles(CCLE.clarified, "mutation")[x,]))){x}},simplify=T)
features.mutation <- intersect(features.mutation.ccle, features.mutation.gdsc)

myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_mutation.RData")
if(!file.exists(myf))
{
  gdsc.sig.mutation <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="mutation", drugs=drugs, features=features.mutation, sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4)
  
  save(gdsc.sig.mutation, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "ccle_sig_mutation.RData")
if(!file.exists(myf))
{
  ccle.sig.mutation <- drugSensitivitySig(pSet=CCLE.clarified, mDataType="mutation", drugs=drugs, features=features.mutation, sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4)
  
  save(ccle.sig.mutation, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_fusion.RData")
if(!file.exists(myf))
{
  gdsc.sig.fusion <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="fusion", drugs=drugs, sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4)
  
  save(gdsc.sig.fusion, file=myf)
}else{
  load(myf)
}

###Try several false discovery rate on all continuous
fdrBarplot(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous", drugs=drugs)
fdrBarplot(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous", drugs=drugs)

fdrBarplot(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous_top100", drugs=drugs, top.ranked=100)
fdrBarplot(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous_top100", drugs=drugs, top.ranked=100)

topBarplot(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous_topx", drugs=drugs, fdr.cut.off=0.05)
topBarplot(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous_topx", drugs=drugs, fdr.cut.off=0.05)


#### Sactter plots for estimates for all continuous
estimatesScatterplot(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous", drugs=drugs, features=features, fdr.cut.off=0.05)
estimatesScatterplot(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous", drugs=drugs, features=features, fdr.cut.off=0.05)
estimatesScatterplot(ccle.sig.rna=ccle.sig.rna.bin, gdsc.sig.rna=gdsc.sig.rna2.bin, cell="all", method="binary", drugs=drugs, features=features, fdr.cut.off=0.05)
estimatesScatterplot(ccle.sig.rna=common.ccle.sig.rna.bin, gdsc.sig.rna=common.gdsc.sig.rna2.bin, cell="common", method="binary", drugs=drugs, features=features, fdr.cut.off=0.05)

###Proportion of dataset-specific biomarkers
continuous_common <- biomarkersSpecificity(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous", drugs=drugs, features=features, fdr.cut.off=0.05)
continuous_all <- biomarkersSpecificity(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous", drugs=drugs, features=features, fdr.cut.off=0.05)
binary_common <- biomarkersSpecificity(ccle.sig.rna=common.ccle.sig.rna.bin, gdsc.sig.rna=common.gdsc.sig.rna2.bin, cell="common", method="binary", drugs=drugs, features=features, fdr.cut.off=0.05)
binary_all <- biomarkersSpecificity(ccle.sig.rna=ccle.sig.rna.bin, gdsc.sig.rna=gdsc.sig.rna2.bin, cell="all", method="binary", drugs=drugs, features=features, fdr.cut.off=0.05)


tt <- rbind("Continuous Common"=continuous_common, "Continuous All"=continuous_all, "Binary Common"=binary_common, "Binary All"=binary_all)

mycol <- RColorBrewer::brewer.pal(n=7, name="Set3")

test <- wilcox.test(continuous_common, binary_common, alternative="greater", paired=TRUE)
test$p.value

test <- wilcox.test(continuous_all, binary_all, alternative="greater", paired=TRUE)
test$p.value


pdf(file=file.path(saveres, "jacard_index.pdf"), height=10 , width=15)
par(mar=c(8,5,2,2))
barplot(tt, beside=TRUE, ylim=c(0,50), las=2, col=mycol[1:4], ylab= "Jaccard Index")
legend("topright", legend=rownames(tt), col= mycol[1:4], bty="n", pch=15)
dev.off()

knownBiomarkersCheck(ccle.sig.rna=ccle.sig.rna, 
                     gdsc.sig.rna=gdsc.sig.rna2, 
                     ccle.sig.mutation=ccle.sig.mutation, 
                     gdsc.sig.mutation=gdsc.sig.mutation, 
                     gdsc.sig.fusion=gdsc.sig.fusion,
                     ccle.sig.cnv=ccle.sig.cnv,
                     gdsc.sig.cnv=gdsc.sig.cnv,
                     cell="all", method="continuous")
knownBiomarkersCheck(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2.bin, cell="common", method="continuous")


###Put all biomarkers in excel file
load(file.path(saveres, "signatures_data.RData"))

drugBasedBiomarkers(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous", drugs=drugs, features=features, cut.off=0.05)
drugBasedBiomarkers(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous", drugs=drugs, features=features, cut.off=0.05)
drugBasedBiomarkers(ccle.sig.rna=ccle.sig.mutation, gdsc.sig.rna=gdsc.sig.mutation, cell="mutation", method="continuous", drugs=drugs, features=features.mutation, cut.off=0.05)
drugBasedBiomarkers(ccle.sig.rna=ccle.sig.cnv, gdsc.sig.rna=gdsc.sig.cnv, cell="cnv", method="continuous", drugs=drugs, features=cnv.fetures, cut.off=0.05)


all.types.biomarkers <- integrateDrugBasedBiomarkers(method="continuous", drugs=drugs, cut.off=0.05)

integrateEstimatesScatterplot(biomarkers=all.types.biomarkers, method="continuous", drugs=drugs, fdr.cut.off=0.05)

fdrBarplot(biomarkers=all.types.biomarkers, method="continuous", drugs=drugs)
integrate.all.validation <- integrateBiomarkersValidation(biomarkers=all.types.biomarkers, method="continuous", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=integrate.all.validation, method="integrate", cell="all")


