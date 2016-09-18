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
library(psych)
library(utils)
library(mclust)
library(epibasix)
library(lsa)
library(xtable)
library(WriteXLS)
library(gplots)
library(vcd)
library(caret)
library(abind)
library(parallel)
library(piano)
library(mRMRe)

# install the latest devel version of the PharmacoGx package
# library(devtools)
# devtools::install_github("bhklab/PharmacoGx", ref="master")
library(PharmacoGx)

#################################################
## global parameters

## set method for downloading
# options(download.file.method="auto")
options(download.file.method="wget")
## change to curl, wget or internal depending on your system

## prevent strings to be converted into factors
options(stringsAsFactors=FALSE)

## set random seed to ensuer reproducibility of the resuls
set.seed(54321)

## number of cpu cores available for the analysis pipeline
## set to 'NULL' if all the available cores should be used
nbcore <- 8
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)

## list of characters to be removed from row and column names
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## directory where all the analysis results will be stored
saveres <- "Output"
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }
  
## number of permutation
nperm <- 1e4

consistency.stats <- list("full"=c("Pearson corr for full data"="pcc.full", "Spearman corr for full data"="scc.full", "Somer's Dxy for full data"="dxy.full"), "sens"=c("Pearson corr for sensitive data"="pcc.sens", "Spearman corr for sensitive data"="scc.sens", "Somer's Dxy for sensitive data"="dxy.sens"), "bin"=c("Matthew corr for binary data"="mcc", "Cramer's V for binary data"="cramerv", "Informedness for binary data"="inform"))

## reasonable cutoff for affymetrix microarrays = 5
affy.cutoff <- 5
## reasonable cutoff for rnaseq = 1
cutoff.rnaseq.ccle <- 1
## reasonable cutoff for AUC = 0.2
auc.cutoff <- 0.2
auc.cytotoxic.cutoff <- 0.5
## reasonable cutoff for -log10(IC50 in nanooM) = 2
ic50.cutoff <- 3
ic50.cytotoxic.cutoff <- 4
## reasonable cutoffs for cnv
cutoff.cnv.amplification <- 1.25
cutoff.cnv.deletion <- 0.75


#################################################
## get pharmacogenomic datasets
#################################################

### download curated pharmacogenomic data for GDSC and CCLE
GDSC <- PharmacoGx::downloadPSet("GDSC", saveDir=file.path(saveres, "PSets"))
CCLE <- PharmacoGx::downloadPSet("CCLE", saveDir=file.path(saveres, "PSets")) 

common <- PharmacoGx::intersectPSet(pSets = list("CCLE"=CCLE, "GDSC"=GDSC), intersectOn = c("cell.lines", "drugs"), strictIntersect = TRUE)
drugs <- c("paclitaxel", "17-AAG", "PD-0325901", "AZD6244", "TAE684", "AZD0530", "PD-0332991", "Crizotinib", "PLX4720", "Nutlin-3", "lapatinib", "Nilotinib", "PHA-665752", "Erlotinib", "Sorafenib")

## drugs in common
## which drugs are cytotoxic? 1 = cytotoxic, 0.5 = cytotoxic-like, 0 = targeted agents
drug.cytotoxic <- c("paclitaxel"=1, "17-AAG"=0.5, "PD-0325901"=0.5, "AZD6244"=0.5, "TAE684"=0, "AZD0530"=0, "PD-0332991"=0, "Crizotinib"=0, "PLX4720"=0, "Nutlin-3"=0, "lapatinib"=0, "Nilotinib"=0, "PHA-665752"=0, "Erlotinib"=0, "Sorafenib"=0)
drugix <- names(drug.cytotoxic)

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
  cell <- PharmacoGx::sensitivityInfo(common$CCLE)[exp, "cellid"]
  drug <- PharmacoGx::sensitivityInfo(common$CCLE)[exp, "drugid"]
  exps <- c(exps, paste(cell, drug, sep="_"))
  PharmacoGx::drugDoseResponseCurve(pSets=common, drug=drug, cellline=cell)
  legend("bottomleft", legend="CCLE", bty="n")
}

for (exp in gdsc.filter$noisy) {
  cell <- sensitivityInfo(common$GDSC)[exp, "cellid"]
  drug <- sensitivityInfo(common$GDSC)[exp, "drugid"]
  if (!(paste(cell, drug, sep="_") %in% exps)) {
    PharmacoGx::drugDoseResponseCurve(pSets=common, drug=drug, cellline=cell)
  }
  legend("bottomleft", legend="GDSC", bty="n")
  
}
dev.off()

ex.gdsc.filter <- PharmacoGx:::filterNoisyCurves(GDSC, nthread=detectCores())
lapply(ex.gdsc.filter, length)
### plotting GDSC noisy curves, curves might be noisy in one or both studies
pdf(file.path(saveres, "GDSC_noisy.pdf"), height=15.5, width=12.5)
par(mfrow = c(4,3))
for (exp in ex.gdsc.filter$noisy) {
  cell <- sensitivityInfo(GDSC)[exp, "cellid"]
  drug <- sensitivityInfo(GDSC)[exp, "drugid"]
  PharmacoGx::drugDoseResponseCurve(pSets=GDSC, drug=drug, cellline=cell)
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
  PharmacoGx::drugDoseResponseCurve(pSets=common, drug=cases[i,2], cellline=cases[i,1], legends.label="ic50_published", plot.type="Fitted", ylim=c(0,130))
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
  PharmacoGx::drugDoseResponseCurve(pSets=common, drug=cases[i,2], cellline=cases[i,1])
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
  PharmacoGx::drugDoseResponseCurve(pSets=GDSC, drug=cases[i,2], cellline=cases[i,1], summarize.replicates = FALSE, plot.type = "Fitted", mycol=RColorBrewer::brewer.pal(n=4, name="Set1")[c(3,4)])
  exp <- which(sensitivityInfo(GDSC)$cellid==cases[i,1] & sensitivityInfo(GDSC)$drugid==cases[i,2])
  ABC <- PharmacoGx::computeABC(conc1 = GDSC@sensitivity$raw[exp[1],,"Dose"], conc2 = GDSC@sensitivity$raw[exp[2],,"Dose"], viability1 = GDSC@sensitivity$raw[exp[1],,"Viability"], viability2 = GDSC@sensitivity$raw[exp[2],,"Viability"])
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
  
  cells <- PharmacoGx::intersectList(phenoInfo(common$CCLE, "rna")$cellid, phenoInfo(common$GDSC, "rna2")$cellid, unique(sensitivityInfo(common$CCLE)$cellid), unique(sensitivityInfo(common$GDSC)$cellid))
  cells <- setdiff(cells, snp.outliers)
  common.clarified <- PharmacoGx::intersectPSet(pSets = list("CCLE"=common$CCLE, "GDSC"=common$GDSC), intersectOn = c("cell.lines", "drugs"), cells=cells)
  save(common.clarified, file=myf)
} else{
  load(myf)
}
 

### plotting all the common experiments between studies(Supplementary File 3)
pdf(file.path(saveres, "Supplementary_File_3.pdf"), height=15.5, width=12.5)
par(mfrow = c(4,3))
for (exp in rownames(sensitivityInfo(common.clarified$CCLE))) {
  cell <- PharmacoGx::sensitivityInfo(common.clarified$CCLE)[exp, "cellid"]
  drug <- PharmacoGx::sensitivityInfo(common.clarified$CCLE)[exp, "drugid"]
  PharmacoGx::drugDoseResponseCurve(pSets=common.clarified, drug=drug, cellline=cell, plot.type="Fitted", ylim=c(0,130))
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
  cells <- cellNames(GDSC.clarified)
  cells <- setdiff(cells, snp.outliers)
  GDSC.clarified <- subsetTo(GDSC.clarified, cells=cells)
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
  cells <- cellNames(CCLE.clarified)
  cells <- setdiff(cells, snp.outliers)
  CCLE.clarified <- subsetTo(CCLE.clarified, cells=cells)
  save(CCLE.clarified, file=myf)
}else {
  load(myf)
}

myf <- file.path(saveres, "signatures_data.RData")
if(!file.exists(myf)){
  drugs <- c("paclitaxel", "17-AAG", "PD-0325901", "AZD6244", "TAE684", "AZD0530", "PD-0332991", "Crizotinib", "PLX4720", "Nutlin-3", "lapatinib", "Nilotinib", "PHA-665752", "Erlotinib", "Sorafenib")
  features <- intersect(rownames(featureInfo(CCLE.clarified, "rna")), rownames(featureInfo(GDSC.clarified, "rna2")))
  save(drugs, features, common.clarified, file=myf)  
} else{
  load(myf)
}

### plotting all internal replicates in GDSC 

replicates <- PharmacoGx::sensitivityInfo(GDSC.clarified)[which(sensitivityInfo(GDSC.clarified)$drugid == "AZD6482"),]
replicates <- replicates[which(duplicated(replicates$cellid)),] #577 it was 617 before removing the noisy ones

ABC <- NULL
for (exps in 1:nrow(replicates)) {
  cell <- replicates[exps, "cellid"]
  drug <- replicates[exps, "drugid"]
  exp <-   which(PharmacoGx::sensitivityInfo(GDSC.clarified)$cellid==cell & sensitivityInfo(GDSC.clarified)$drugid==drug)
  ABC <- c(ABC, PharmacoGx::computeABC(conc1 = GDSC.clarified@sensitivity$raw[exp[1],,"Dose"], conc2 = GDSC.clarified@sensitivity$raw[exp[2],,"Dose"], viability1 = GDSC.clarified@sensitivity$raw[exp[1],,"Viability"], viability2 = GDSC.clarified@sensitivity$raw[exp[2],,"Viability"]))
}
biological_replicates_ABC <- ABC[which(!is.na(ABC))]
pdf(file.path(saveres, "Supplementary_File_4.pdf"), height=15.5, width=12.5)
par(mfrow = c(4,3))

for (i in order(ABC, decreasing = T)) {
  exp <- rownames(replicates)[i]
  cell <- replicates[exp, "cellid"]
  drug <- replicates[exp, "drugid"]
  PharmacoGx::drugDoseResponseCurve(pSets=GDSC.clarified, drug=drug, cellline=cell, summarize.replicates = FALSE, plot.type = "Fitted", mycol=RColorBrewer::brewer.pal(n=4, name="Set1")[c(3,4)])
  legend("bottomleft", legend = sprintf("ABC= %s", round(ABC[i], digits=2)), bty="n")
}

dev.off()


##boxplot of ABC across celllines and across drugs

ABC <- NULL
for (exp.ccle in rownames(sensitivityInfo(common.clarified$CCLE))) {
  cell <- PharmacoGx::sensitivityInfo(common.clarified$CCLE)[exp.ccle, "cellid"]
  drug <- PharmacoGx::sensitivityInfo(common.clarified$CCLE)[exp.ccle, "drugid"]
  exp.gdsc <-   rownames(PharmacoGx::sensitivityInfo(common.clarified$GDSC))[which(PharmacoGx::sensitivityInfo(common.clarified$GDSC)$cellid==cell & PharmacoGx::sensitivityInfo(common.clarified$GDSC)$drugid==drug)]
  ABC <- c(ABC, PharmacoGx::computeABC(conc1 = common.clarified$CCLE@sensitivity$raw[exp.ccle,,"Dose"], 
                           conc2 = common.clarified$GDSC@sensitivity$raw[exp.gdsc,,"Dose"], 
                           viability1 = common.clarified$CCLE@sensitivity$raw[exp.ccle,,"Viability"], 
                           viability2 = common.clarified$GDSC@sensitivity$raw[exp.gdsc,,"Viability"]))
}
names(ABC) <- rownames(PharmacoGx::sensitivityInfo(common.clarified$CCLE))
ABC <- cbind(ABC, "cellid"=PharmacoGx::sensitivityInfo(common.clarified$CCLE)[names(ABC),"cellid"], "drugid"=PharmacoGx::sensitivityInfo(common.clarified$CCLE)[names(ABC),"drugid"])

Inter_studies_ABC <- ABC[which(!is.na(ABC[,"ABC"])),"ABC"]
###wilcoxon test
biological_replicates_ABC <- as.numeric(biological_replicates_ABC)
Inter_studies_ABC <- as.numeric(Inter_studies_ABC)
tt <- matrix(NA, ncol=2, nrow=max(length(Inter_studies_ABC), length(biological_replicates_ABC)))
colnames(tt) <- c("Intra", "Inter")
tt[1:length(biological_replicates_ABC),"Intra"] <- biological_replicates_ABC
tt[1:length(Inter_studies_ABC),"Inter"] <- Inter_studies_ABC

pdf(file.path(saveres, "inter_intra_abc.pdf"), height=7, width=7)
par(mar=c(9,5,5,2))
boxplot(tt, las = 2, col = "gray", cex.lab=1, cex.axis=1, pch=19, ylab="ABC", outpch=20, outcex=0.5)
dev.off()
test <- wilcox.test(Inter_studies_ABC, biological_replicates_ABC, paired=FALSE)
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
unique(PharmacoGx::sensitivityInfo(common.clarified$GDSC)[which(PharmacoGx::sensitivityInfo(common.clarified$GDSC)$cellid %in% snp.outliers),"cellid"])
sort(unique(PharmacoGx::sensitivityInfo(common.clarified$GDSC)[which(PharmacoGx::sensitivityInfo(common.clarified$GDSC)$cellid %in% snp.outliers),"drugid"]))



## min and max concentrations tested in both datasets
range.concentration.gdsc <- range.concentration.ccle <- range.concentration <- NULL
## GDSC
for (i in 1:nrow(drugInfo(GDSC.clarified))) {
  ## GDSC
  dix <- !is.na(sensitivityInfo(GDSC.clarified)[ , "drugid"]) & sensitivityInfo(GDSC.clarified)[ , "drugid"] == rownames(drugInfo(GDSC.clarified))[i]
  tt <- GDSC.clarified@sensitivity$raw[dix, ,"Dose"]
  range.concentration.gdsc <- rbind(range.concentration.gdsc, c(min(as.numeric(tt[ , 1]), na.rm=TRUE), max(as.numeric(tt[ , ncol(tt)]), na.rm=TRUE)))
}
dimnames(range.concentration.gdsc) <- list(rownames(drugInfo(GDSC.clarified)), c("Dose.min", "Dose.max"))
## CCLE
for (i in 1:nrow(drugInfo(CCLE.clarified))) {
  ## CCLE
  dix <- !is.na(sensitivityInfo(CCLE.clarified)[ , "drugid"]) & sensitivityInfo(CCLE.clarified)[ , "drugid"] == rownames(drugInfo(CCLE.clarified))[i]
  tt <- CCLE.clarified@sensitivity$raw[dix, ,"Dose"]
  range.concentration.ccle <- rbind(range.concentration.ccle, c(min(as.numeric(tt[ , 1]), na.rm=TRUE), max(as.numeric(tt[ , ncol(tt)]), na.rm=TRUE)))
}
dimnames(range.concentration.ccle) <- list(rownames(drugInfo(CCLE.clarified)), c("Dose.min", "Dose.max"))

## distribution of expression data
grDevices::pdf(file=file.path(saveres, sprintf("distribution_rna.pdf")), height=4, width=12)
par(mfrow=c(2, 3), cex=0.8, las=1)

yylim <- c(0, 0.55)

## common genes
gix <- paste("geneid", intersectList(list(featureInfo(pSet=common.clarified$GDSC, mDataType="rna2")[ , "EnsemblGeneId"], featureInfo(pSet=common.clarified$CCLE, mDataType="rna")[ , "EnsemblGeneId"], featureInfo(pSet=common.clarified$CCLE, mDataType="rnaseq")[ , "EnsemblGeneId"])), sep=".")

mdd <- molecularProfiles(pSet=common.clarified$GDSC, mDataType="rna2")
fdd <- paste("geneid", featureInfo(pSet=common.clarified$GDSC, mDataType="rna2")[ , "EnsemblGeneId"], sep=".")
mdd <- mdd[!is.na(fdd) & is.element(fdd, gix), , drop=FALSE]
## identify cutoff
myf <- file.path(saveres, "mclust_rna_gdsc.RData")
if (!file.exists(myf)) {
  mclust.rna.gdsc <- Mclust(data=as.numeric(mdd), G=2, modelNames="V")
  cutoff.rna.gdsc <- max(1, ceiling(qnorm(p=0.1, mean=mclust.rna.gdsc$parameters$mean[2], sd=sqrt(mclust.rna.gdsc$parameters$variance$sigmasq[2]))))
  save(list=c("mclust.rna.gdsc", "cutoff.rna.gdsc"), compress=TRUE, file=myf)
} else {
  load(myf)
}
hist(as.numeric(mdd), breaks=100, main="Distribution of GDSC Affymetrix data", freq=FALSE, xlab="Affymetrix HG-U219 expression values", ylim=yylim)
abline(v=cutoff.rna.gdsc, col="red", lty=2, lwd=1)

mdd <- molecularProfiles(pSet=common.clarified$CCLE, mDataType="rna")
fdd <- paste("geneid", featureInfo(pSet=common.clarified$CCLE, mDataType="rna")[ , "EnsemblGeneId"], sep=".")
mdd <- mdd[!is.na(fdd) & is.element(fdd, gix), , drop=FALSE]
## identify cutoff
myf <- file.path(saveres, "mclust_rna_ccle.RData")
if (!file.exists(myf)) {
  mclust.rna.ccle <- Mclust(data=as.numeric(mdd), G=2, modelNames="V")
  cutoff.rna.ccle <- max(1, ceiling(qnorm(p=0.1, mean=mclust.rna.ccle$parameters$mean[2], sd=sqrt(mclust.rna.ccle$parameters$variance$sigmasq[2]))))
  save(list=c("mclust.rna.ccle", "cutoff.rna.ccle"), compress=TRUE, file=myf)
} else {
  load(myf)
}
hist(as.numeric(mdd), breaks=100, main="Distribution of CCLE Affymetrix data", freq=FALSE, xlab="Affymetrix HG-U133PLUS2 expression values", ylim=yylim)
abline(v=cutoff.rna.ccle, col="red", lty=2, lwd=1)

mdd <- molecularProfiles(pSet=common.clarified$CCLE, mDataType="rnaseq")
fdd <- paste("geneid", featureInfo(pSet=common.clarified$CCLE, mDataType="rnaseq")[ , "EnsemblGeneId"], sep=".")
mdd <- mdd[!is.na(fdd) & is.element(fdd, gix), , drop=FALSE]
## identify cutoff
myf <- file.path(saveres, "mclust_rnaseq_ccle.RData")
if (!file.exists(myf)) {
  mclust.rnaseq.ccle <- Mclust(data=as.numeric(mdd), G=2, modelNames="V")
  cutoff.rnaseq.ccle <- max(1, ceiling(qnorm(p=0.1, mean=mclust.rnaseq.ccle$parameters$mean[2], sd=sqrt(mclust.rnaseq.ccle$parameters$variance$sigmasq[2]))))
  save(list=c("mclust.rnaseq.ccle", "cutoff.rnaseq.ccle"), compress=TRUE, file=myf)
} else {
  load(myf)
}
hist(as.numeric(mdd), breaks=100, main="Distribution of CCLE RNA-seq data", freq=FALSE, xlab="Illumina RNA-seq expression values", ylim=yylim)
abline(v=cutoff.rnaseq.ccle, col="red", lty=2, lwd=1)

## distribution of cnv data
hist(as.numeric(molecularProfiles(pSet=common.clarified$CCLE, mDataType="cnv")), breaks=150, main="Distribution of CCLE CNV data", freq=FALSE, xlab="CNV ratios")
abline(v=cutoff.cnv.amplification, col="red", lty=2)
abline(v=cutoff.cnv.deletion, col="red", lty=2)

hist(as.numeric(molecularProfiles(pSet=common.clarified$GDSC, mDataType="cnv")), breaks=150, main="Distribution of GDSC CNV data", freq=FALSE, xlab="CNV ratios")
abline(v=cutoff.cnv.amplification, col="red", lty=2)
abline(v=cutoff.cnv.deletion, col="red", lty=2)

dev.off()

## distribution of sensitivity data
grDevices::pdf(file=file.path(saveres, sprintf("distribution_sensitivity.pdf")), height=8, width=8)
par(mfrow=c(2, 2), cex=0.8, las=1)

hist(as.numeric(summarizeSensitivityProfiles(pSet=common.clarified$GDSC, sensitivity.measure="auc_published", summary.stat="median")), breaks=50, main="Distribution of GDSC AUC data", freq=FALSE, xlab="AUC values as published", xlim=c(0, 1))
abline(v=auc.cutoff, col="red", lty=2)
abline(v=auc.cytotoxic.cutoff, col="orange", lty=2)

hist(as.numeric(summarizeSensitivityProfiles(pSet=common.clarified$CCLE, sensitivity.measure="auc_published", summary.stat="median")), breaks=50, main="Distribution of CCLE AUC data", freq=FALSE, xlab="AUC values as published", xlim=c(0, 1))
abline(v=auc.cutoff, col="red", lty=2)
abline(v=auc.cytotoxic.cutoff, col="orange", lty=2)

hist(- log10(as.numeric(summarizeSensitivityProfiles(pSet=common.clarified$GDSC, sensitivity.measure="ic50_published", summary.stat="median") / 1000)), breaks=50, main="Distribution of GDSC IC50 data", freq=FALSE, xlab="- log10(IC50 nanoMolar) values as published", xlim=c(0, 8))
abline(v=ic50.cutoff, col="red", lty=2)
abline(v=ic50.cytotoxic.cutoff, col="orange", lty=2)

hist(- log10(as.numeric(summarizeSensitivityProfiles(pSet=common.clarified$CCLE, sensitivity.measure="ic50_published", summary.stat="median")) / 1000), breaks=50, main="Distribution of CCLE IC50 data", freq=FALSE, xlab="- log10(IC50 nanoMolar) values as published", xlim=c(0, 8))
abline(v=ic50.cutoff, col="red", lty=2)
abline(v=ic50.cytotoxic.cutoff, col="orange", lty=2)

dev.off()


#################################################
## MAD of drug sensitivity data
#################################################

## MAD of cytotoxic drugs vs the rest in full studies
drug.mad.gdsc <- apply(PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC.clarified, sensitivity.measure="auc_published", summary.stat="median"), 1, mad, na.rm=TRUE)
iix.gdsc <- factor(!is.na(PharmacoGx::drugInfo(GDSC.clarified)[ , "Drug.class.II"]) & PharmacoGx::drugInfo(GDSC.clarified)[ , "Drug.class.II"] == "Cytotoxic", levels=c(TRUE, FALSE))
levels(iix.gdsc) <- c("Cytotoxic", "Targeted")
drug.mad.ccle <- apply(PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE.clarified, sensitivity.measure="auc_published", summary.stat="median"), 1, mad, na.rm=TRUE)
iix.ccle <- factor(!is.na(PharmacoGx::drugInfo(CCLE.clarified)[ , "Class"]) & PharmacoGx::drugInfo(CCLE.clarified)[ , "Class"] == "Cytotoxic", levels=c(TRUE, FALSE))
levels(iix.ccle) <- c("Cytotoxic", "Targeted")

## find optimal cutoff to descriminate cytotoxic vs targeted drugs based on AUC MAD
rroc.gdsc <- pROC::roc(formula=drug.type ~ drug.mad, data=data.frame("drug.mad"=drug.mad.gdsc, "drug.type"=iix.gdsc))
threshold.gdsc <- pROC::coords(roc=rroc.gdsc, x="best", best.method="youden")[1]
rroc.ccle <- pROC::roc(formula=drug.type ~ drug.mad, data=data.frame("drug.mad"=drug.mad.ccle, "drug.type"=iix.ccle))
threshold.ccle <- pROC::coords(roc=rroc.ccle, x="best", best.method="youden")[1]
mad.cytotoxic.cutoff <- round(mean(threshold.gdsc, threshold.ccle) * 100) / 100


pdf(file.path(saveres, "auc_mad_vs_cytotoxic_full.pdf"), width=10, height=5)
par(mfrow=c(1, 2), mar=c(5, 4, 3, 2) + 0.1)
yylim <- c(0, ceiling(max(c(drug.mad.gdsc, drug.mad.ccle), na.rm=TRUE) * 10) / 10)
boxplot(drug.mad.gdsc ~ iix.gdsc, ylim=yylim, col="lightgrey", border="black", pars=list(outcol="black", outpch=20, outcex=0.5), ylab="AUC published (GDSC)", main="")
abline(h=mad.cytotoxic.cutoff, lty=2, col="red")
mtext(side=3, at=par("usr")[1] - 0.33, text=LETTERS[1], line=2, font=2, cex=0.8)
boxplot(drug.mad.ccle ~ iix.ccle, ylim=yylim, col="lightgrey", border="black", pars=list(outcol="black", outpch=20, outcex=0.5), ylab="AUC published (CCLE)", main="")
abline(h=mad.cytotoxic.cutoff, lty=2, col="red")
mtext(side=3, at=par("usr")[1] - 0.33, text=LETTERS[2], line=2, font=2, cex=0.8)
dev.off()

## MAD of AUC in common cell lines in both studies
drug.mad.gdsc <- apply(summarizeSensitivityProfiles(pSet=common.clarified$GDSC, sensitivity.measure="auc_published", summary.stat="median"), 1, mad, na.rm=TRUE)[drugix]
drug.mad.ccle <- apply(summarizeSensitivityProfiles(pSet=common.clarified$CCLE, sensitivity.measure="auc_published", summary.stat="median"), 1, mad, na.rm=TRUE)[drugix]

pdf(file.path(saveres, "auc_mad_vs_mad_common.pdf"), width=5, height=5)
par(mar=c(5, 4, 1, 2) + 0.1, cex=0.8)
xxlim <- c(floor(min(drug.mad.gdsc, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.gdsc, na.rm=TRUE) * 1200) / 1000)
yylim <- c(floor(min(drug.mad.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.ccle, na.rm=TRUE) * 1200) / 1000)
llim <- c(min(xxlim[1], yylim[1]), max(xxlim[2], yylim[2]))
plot(x=drug.mad.gdsc, y=drug.mad.ccle, xlim=llim, ylim=llim, pch=20, col=blues9[7], xlab="MAD of AUC published (GDSC)", ylab="MAD of AUC published (CCLE)")
abline(a=0, b=1, col="darkgrey", lwd=0.5)
abline(h=mad.cytotoxic.cutoff, col="red", lty=2, lwd=0.5)
abline(v=mad.cytotoxic.cutoff, col="red", lty=2, lwd=0.5)
text(x=drug.mad.gdsc, y=drug.mad.ccle, labels=drugix, cex=0.7, font=1, srt=45, pos=4)
# wordcloud::textplot(x=drug.mad.gdsc, y=drug.mad.ccle, words=drugix, new=FALSE, cex=0.8, srt=45)
dev.off()


################################################
## consistency of published vs recomputed drug sensitivity data
################################################

myf <- file.path(saveres, "auc_ic50_published_recomputed.RData")
if (!file.exists(myf)) {
  ## GDSC AUC
  x0.gdsc <- PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC.clarified, sensitivity.measure="auc_published", summary.stat="median")
  x1.gdsc <- PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC.clarified, sensitivity.measure="auc_recomputed", summary.stat="median")
  ## GDSC IC50
  xx0.gdsc <- PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC.clarified, sensitivity.measure="ic50_published", summary.stat="median")
  xx1.gdsc <- summarizeSensitivityProfiles(pSet=GDSC.clarified, sensitivity.measure="ic50_recomputed", summary.stat="median")
  for (i in 1:nrow(xx0.gdsc)) {
    ## truncate ic50= and ic50=Inf by the min and max concentrations tested in each study
    ## ic50 = 0 when the drug dose-response curve starts below 50% viability
    ## ic50 = Inf when 50% viability is never reached
    xx0.gdsc.i <- xx0.gdsc[i, ]
    xx0.gdsc.i[xx0.gdsc.i < range.concentration.gdsc[i, 1]] <- range.concentration.gdsc[i, 1]
    xx0.gdsc.i[xx0.gdsc.i > range.concentration.gdsc[i, 2]] <- range.concentration.gdsc[i, 2]
    xx0.gdsc[i, ] <- - log10(xx0.gdsc.i / 1000)
    xx1.gdsc.i <- xx1.gdsc[i, ]
    xx1.gdsc.i[xx1.gdsc.i < range.concentration.gdsc[i, 1]] <- range.concentration.gdsc[i, 1]
    xx1.gdsc.i[xx1.gdsc.i > range.concentration.gdsc[i, 2]] <- range.concentration.gdsc[i, 2]
    xx1.gdsc[i, ] <- - log10(xx1.gdsc.i / 1000)
  }
  ## CCLE AUC
  x0.ccle <- PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE.clarified, sensitivity.measure="auc_published", summary.stat="median")
  x1.ccle <- PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE.clarified, sensitivity.measure="auc_recomputed", summary.stat="median")
  ## CCLE IC50
  xx0.ccle <- PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE.clarified, sensitivity.measure="ic50_published", summary.stat="median")
  xx1.ccle <- PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE.clarified, sensitivity.measure="ic50_recomputed", summary.stat="median")
  for (i in 1:nrow(xx0.ccle)) {
    ## truncate ic50= and ic50=Inf by the min and max concentrations tested in each study
    ## ic50 = 0 when the drug dose-response curve starts below 50% viability
    ## ic50 = Inf when 50% viability is never reached
    xx0.ccle.i <- xx0.ccle[i, ]
    xx0.ccle.i[xx0.ccle.i < range.concentration.ccle[i, 1]] <- range.concentration.ccle[i, 1]
    xx0.ccle.i[xx0.ccle.i > range.concentration.ccle[i, 2]] <- range.concentration.ccle[i, 2]
    xx0.ccle[i, ] <- - log10(xx0.ccle.i / 1000)
    xx1.ccle.i <- xx1.ccle[i, ]
    xx1.ccle.i[xx1.ccle.i < range.concentration.ccle[i, 1]] <- range.concentration.ccle[i, 1]
    xx1.ccle.i[xx1.ccle.i > range.concentration.ccle[i, 2]] <- range.concentration.ccle[i, 2]
    xx1.ccle[i, ] <- - log10(xx1.ccle.i / 1000)
  }

  dd <- list("GDSC.AUC.PUBLISHED"=data.frame(x0.gdsc, check.names=FALSE), "GDSC.AUC.RECOMPUTED"=data.frame(x1.gdsc, check.names=FALSE), "GDSC.IC50.PUBLISHED"=data.frame(xx0.gdsc, check.names=FALSE), "GDSC.IC50.RECOMPUTED"=data.frame(xx1.gdsc, check.names=FALSE), "CCLE.AUC.PUBLISHED"=data.frame(x0.ccle, check.names=FALSE), "CCLE.AUC.RECOMPUTED"=data.frame(x1.ccle, check.names=FALSE), "CCLE.IC50.PUBLISHED"=data.frame(xx0.ccle, check.names=FALSE), "CCLE.IC50.RECOMPUTED"=data.frame(xx1.ccle, check.names=FALSE))
  WriteXLS::WriteXLS(x="dd", row.names=TRUE, ExcelFileName=file.path(saveres, "auc_ic50_published_recomputed.xlsx"))
  save(list=c("x0.gdsc", "x1.gdsc", "xx0.gdsc", "xx1.gdsc", "x0.ccle", "x1.ccle", "xx0.ccle", "xx1.ccle"), compress=TRUE, file=myf)
} else {
  load(myf)
}

grDevices::pdf(file=file.path(saveres, sprintf("sensitivity_published_recomputed.pdf")), height=10, width=10)
par(mfrow=c(2, 2), mar=c(5, 4, 3, 2) + 0.1, cex=0.8, las=1)
  
mylim <- c(0, 1)
mysub <- NULL
myScatterPlot(x=as.numeric(x0.gdsc), y=as.numeric(x1.gdsc), method=c("transparent"), transparency=0.05, pch=16, minp=50, xlim=mylim, ylim=mylim, main="Published vs. recomputed AUC values in GDSC", cex.sub=0.7, sub=mysub, xlab="Published AUC in GDSC", ylab="Recomputed AUC in GDSC")
tt <- cor.test(x=as.numeric(x0.gdsc), y=as.numeric(x1.gdsc), method="spearman", use="complete.obs", exact=FALSE)
legend("bottomright", legend=sprintf("SCC = %.2g", tt$estimate), col=c("white"), pch=NA, bty="n", cex=1.2)
mtext(side=3, at=par("usr")[1]-0.13, text=LETTERS[1], line=2, font=2, cex=0.8)

mylim <- c(0, 1)
mysub <- NULL
myScatterPlot(x=as.numeric(x0.ccle), y=as.numeric(x1.ccle), method=c("transparent"), transparency=0.05, pch=16, minp=50, xlim=mylim, ylim=mylim, main="Published vs. recomputed AUC values in CCLE", cex.sub=0.7, sub=mysub, xlab="Piblished AUC in CCLE", ylab="Recomputed AUC in CCLE")
tt <- cor.test(x=as.numeric(x0.ccle), y=as.numeric(x1.ccle), method="spearman", use="complete.obs", exact=FALSE)
legend("bottomright", legend=sprintf("SCC = %.2g", tt$estimate), col=c("white"), pch=NA, bty="n", cex=1.2)
mtext(side=3, at=par("usr")[1]-0.11, text=LETTERS[2], line=2, font=2, cex=0.8)

mylim <- range(c(xx0.gdsc, xx1.gdsc), na.rm=TRUE)
mysub <- NULL
myScatterPlot(x=as.numeric(xx0.gdsc), y=as.numeric(xx1.gdsc), method=c("transparent"), transparency=0.05, pch=16, minp=50, xlim=mylim, ylim=mylim, main="Published vs. recomputed IC50 values in GDSC", cex.sub=0.7, sub=mysub, xlab="- log10(published IC50) in GDSC", ylab="- log10(recomputed IC50) in GDSC")
tt <- cor.test(x=as.numeric(xx0.gdsc), y=as.numeric(xx1.gdsc), method="spearman", use="complete.obs", exact=FALSE)
legend("bottomright", legend=sprintf("SCC = %.2g", tt$estimate), col=c("white"), pch=NA, bty="n", cex=1.2)
mtext(side=3, at=par("usr")[1]-1, text=LETTERS[3], line=2, font=2, cex=0.8)

mylim <- range(c(xx0.ccle, xx1.ccle), na.rm=TRUE)
mysub <- NULL
myScatterPlot(x=as.numeric(xx0.ccle), y=as.numeric(xx1.ccle), method=c("transparent"), transparency=0.05, pch=16, minp=50, xlim=mylim, ylim=mylim, main="Published vs. recomputed IC50 values in CCLE", cex.sub=0.7, sub=mysub, xlab="- log10(published IC50) in CCLE", ylab="- log10(recomputed IC50) in CCLE")
tt <- cor.test(x=as.numeric(xx0.ccle), y=as.numeric(xx1.ccle), method="spearman", use="complete.obs", exact=FALSE)
legend("bottomright", legend=sprintf("SCC = %.2g", tt$estimate), col=c("white"), pch=NA, bty="n", cex=1.2)
mtext(side=3, at=par("usr")[1]-0.5, text=LETTERS[4], line=2, font=2, cex=0.8)

dev.off()

################################################
## consistency statistics per drug
################################################

consisn <- NULL

########################
## AUC as published
########################

auc.gdsc.common <- PharmacoGx::summarizeSensitivityProfiles(pSet=common.clarified$GDSC, sensitivity.measure="auc_published", summary.stat="median")
auc.ccle.common <- PharmacoGx::summarizeSensitivityProfiles(pSet=common.clarified$CCLE, sensitivity.measure="auc_published", summary.stat="median")

## compute consistency for auc published
myf <- file.path(saveres, sprintf("auc_%s_consistency.RData", "published"))
if (!file.exists(myf)) {
  consis <- computeConsistencySensitivity(x=auc.gdsc.common, y=auc.ccle.common, type="auc", drugs=drug.cytotoxic, cutoff=auc.cutoff, cutoff.cytotoxic=auc.cytotoxic.cutoff)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn <- c(consisn, list("AUC.PUBLISHED"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
## scatterplots
grDevices::pdf(file=file.path(saveres, sprintf("auc_%s_bin_lines.pdf", "published")), height=13, width=12)
par(mfrow=c(4, 4), mar=c(5, 4, 3, 2) + 0.1, cex=0.8, las=1)
for (i in drugix) {
  ## binarization
  cc <- ifelse(drug.cytotoxic[i] == 1, auc.cytotoxic.cutoff, auc.cutoff)
  auc.gdsc.common.bin <- factor(ifelse(auc.gdsc.common[i, ] > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))
  auc.ccle.common.bin <- factor(ifelse(auc.ccle.common[i, ] > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))  
  ccix <- complete.cases(auc.gdsc.common[i, ], auc.ccle.common[i, ])
  # mylim <- range(c(auc.gdsc.common[i, ccix], auc.ccle.common[i, ccix]), na.rm=TRUE)
  # if (mylim[2] < (auc.cutoff + 0.1)) {
  #  mylim[2] <- ceiling(auc.cutoff + 0.1)
  # }
  # if (mylim[1] > (auc.cutoff - 0.1)) {
  #  mylim[1] <- floor(auc.cutoff - 0.1)
  # }
  mylim <- c(0, 1)
  mycol <- rep("#2B83BA", length(auc.gdsc.common.bin))
  mycol[complete.cases(auc.gdsc.common.bin, auc.ccle.common.bin) & auc.gdsc.common.bin == "resistant" & auc.ccle.common.bin == "resistant"] <- gplots::col2hex("#D7191C")
  mycol[complete.cases(auc.gdsc.common.bin, auc.ccle.common.bin) & auc.gdsc.common.bin == "sensitive" & auc.ccle.common.bin == "resistant"] <- gplots::col2hex("#FDAE61")
  mycol[complete.cases(auc.gdsc.common.bin, auc.ccle.common.bin) & auc.gdsc.common.bin == "resistant" & auc.ccle.common.bin == "sensitive"] <- gplots::col2hex("#ABDDA4")
  # mysub <- sprintf("PCC=%.2g(%.2g), SCC=%.2g(%.2g), DXY=%.2g(%.2g)\nMCC=%.2g, CRAMERV=%.2g, INFORMEDNESS=%.2g", consis[i, "pcc.full", "estimate"], consis[i, "pcc.sens", "estimate"], consis[i, "scc.full", "estimate"], consis[i, "scc.sens", "estimate"], consis[i, "dxy.full", "estimate"], consis[i, "dxy.sens", "estimate"], consis[i, "mcc", "estimate"], consis[i, "cramerv", "estimate"], consis[i, "inform", "estimate"])
  mysub <- NULL
  myScatterPlot(x=auc.gdsc.common[i, ], y=auc.ccle.common[i, ], method=c("transparent"), transparency=0.50, pch=16, minp=50, col=mycol, xlim=mylim, ylim=mylim, main=i, cex.sub=0.7, sub=mysub, xlab="AUC GDSC", ylab="AUC CCLE")
  abline(v=cc, col="grey", lty=2)
  abline(h=cc, col="grey", lty=2)
  ## add diagonal line
  abline(a=0, b=1, col="grey", lty=1, lwd=0.5)
}
## legend
plot.new()
par(mar=c(0.1, 0.1, 0.1, 0.1))
legend("center", legend=c("Both sensitive", "", "Both resistant", "", "GDSC sensitive / CCLE resistant", "", "GDSC resistant / CCLE sensitive"), col=c("#2B83BA", "white", "#D7191C", "white", "#FDAE61", "white", "#ABDDA4"), pch=16, bty="n", cex=1, pt.cex=1.5)
dev.off()

auc.published.consistency <- consis[drugix, , , drop=FALSE]
## save consistency statistics
xt <- auc.published.consistency[ , unlist(consistency.stats), "estimate"]
colnames(xt) <- toupper(colnames(xt))
xtable::print.xtable(xtable::xtable(xt, digits=2), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=file.path(saveres, "auc_published_consistency.tex"), append=FALSE)

########################
## AUC recomputed
########################

auc.gdsc.common <- summarizeSensitivityProfiles(pSet=common.clarified$GDSC, sensitivity.measure="auc_recomputed", summary.stat="median")
auc.ccle.common <- summarizeSensitivityProfiles(pSet=common.clarified$CCLE, sensitivity.measure="auc_recomputed", summary.stat="median")

myf <- file.path(saveres, sprintf("auc_%s_consistency.RData", "recomputed"))
if (!file.exists(myf)) {
  consis <- computeConsistencySensitivity(x=auc.gdsc.common, y=auc.ccle.common, type="auc", drugs=drug.cytotoxic, cutoff=auc.cutoff, cutoff.cytotoxic=auc.cytotoxic.cutoff)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn <- c(consisn, list("AUC.RECOMPUTED"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
## scatterplots
grDevices::pdf(file=file.path(saveres, sprintf("auc_%s_bin_lines.pdf", "recomputed")), height=13, width=12)
par(mfrow=c(4, 4), mar=c(5, 4, 3, 2) + 0.1, cex=0.8, las=1)
for (i in drugix) {
  ## binarization
  cc <- ifelse(drug.cytotoxic[i] == 1, auc.cytotoxic.cutoff, auc.cutoff)
  auc.gdsc.common.bin <- factor(ifelse(auc.gdsc.common[i, ] > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))
  auc.ccle.common.bin <- factor(ifelse(auc.ccle.common[i, ] > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))  
  ccix <- complete.cases(auc.gdsc.common[i, ], auc.ccle.common[i, ])
  mylim <- c(0, 1)
  mycol <- rep("#2B83BA", length(auc.gdsc.common.bin))
  mycol[complete.cases(auc.gdsc.common.bin, auc.ccle.common.bin) & auc.gdsc.common.bin == "resistant" & auc.ccle.common.bin == "resistant"] <- gplots::col2hex("#D7191C")
  mycol[complete.cases(auc.gdsc.common.bin, auc.ccle.common.bin) & auc.gdsc.common.bin == "sensitive" & auc.ccle.common.bin == "resistant"] <- gplots::col2hex("#FDAE61")
  mycol[complete.cases(auc.gdsc.common.bin, auc.ccle.common.bin) & auc.gdsc.common.bin == "resistant" & auc.ccle.common.bin == "sensitive"] <- gplots::col2hex("#ABDDA4")
  mysub <- NULL
  myScatterPlot(x=auc.gdsc.common[i, ], y=auc.ccle.common[i, ], method=c("transparent"), transparency=0.50, pch=16, minp=50, col=mycol, xlim=mylim, ylim=mylim, main=i, cex.sub=0.7, sub=mysub, xlab="AUC GDSC", ylab="AUC CCLE")
  abline(v=cc, col="grey", lty=2)
  abline(h=cc, col="grey", lty=2)
  ## add diagonal line
  abline(a=0, b=1, col="grey", lty=1, lwd=0.5)
}
## legend
plot.new()
par(mar=c(0.1, 0.1, 0.1, 0.1))
legend("center", legend=c("Both sensitive", "", "Both resistant", "", "GDSC sensitive / CCLE resistant", "", "GDSC resistant / CCLE sensitive"), col=c("#2B83BA", "white", "#D7191C", "white", "#FDAE61", "white", "#ABDDA4"), pch=16, bty="n", cex=1, pt.cex=1.5)
dev.off()

auc.recomputed.consistency <- consis[drugix, , , drop=FALSE]
## save consistency statistics
xt <- auc.recomputed.consistency[ , unlist(consistency.stats), "estimate"]
colnames(xt) <- toupper(colnames(xt))
xtable::print.xtable(xtable::xtable(xt, digits=2), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=file.path(saveres, "auc_recomputed_consistency.tex"), append=FALSE)

########################
## AUC recomputed star
########################
common.star <- PharmacoGx::intersectPSet(pSets = list("CCLE"=CCLE, "GDSC"=GDSC), intersectOn = c("cell.lines", "drugs", "concentrations"), strictIntersect=TRUE)

auc.gdsc.common <- summarizeSensitivityProfiles(pSet=common.star$GDSC, sensitivity.measure="auc_recomputed", summary.stat="median")
auc.ccle.common <- summarizeSensitivityProfiles(pSet=common.star$CCLE, sensitivity.measure="auc_recomputed", summary.stat="median")

## compute consistency for auc star
myf <- file.path(saveres, sprintf("auc_%s_consistency.RData", "star"))
if (!file.exists(myf)) {
  consis <- computeConsistencySensitivity(x=auc.gdsc.common, y=auc.ccle.common, type="auc", drugs=drug.cytotoxic, cutoff=auc.cutoff, cutoff.cytotoxic=auc.cytotoxic.cutoff)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn <- c(consisn, list("AUC.STAR.RECOMPUTED"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
## scatterplots
grDevices::pdf(file=file.path(saveres, sprintf("auc_%s_bin_lines.pdf", "star")), height=13, width=12)
par(mfrow=c(4, 4), mar=c(5, 4, 3, 2) + 0.1, cex=0.8, las=1)
for (i in drugix) {
  ## binarization
  cc <- ifelse(drug.cytotoxic[i] == 1, auc.cytotoxic.cutoff, auc.cutoff)
  auc.gdsc.common.bin <- factor(ifelse(auc.gdsc.common[i, ] > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))
  auc.ccle.common.bin <- factor(ifelse(auc.ccle.common[i, ] > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))  
  ccix <- complete.cases(auc.gdsc.common[i, ], auc.ccle.common[i, ])
  mylim <- c(0, 1)
  mycol <- rep("#2B83BA", length(auc.gdsc.common.bin))
  mycol[complete.cases(auc.gdsc.common.bin, auc.ccle.common.bin) & auc.gdsc.common.bin == "resistant" & auc.ccle.common.bin == "resistant"] <- gplots::col2hex("#D7191C")
  mycol[complete.cases(auc.gdsc.common.bin, auc.ccle.common.bin) & auc.gdsc.common.bin == "sensitive" & auc.ccle.common.bin == "resistant"] <- gplots::col2hex("#FDAE61")
  mycol[complete.cases(auc.gdsc.common.bin, auc.ccle.common.bin) & auc.gdsc.common.bin == "resistant" & auc.ccle.common.bin == "sensitive"] <- gplots::col2hex("#ABDDA4")
  mysub <- NULL
  myScatterPlot(x=auc.gdsc.common[i, ], y=auc.ccle.common[i, ], method=c("transparent"), transparency=0.50, pch=16, minp=50, col=mycol, xlim=mylim, ylim=mylim, main=i, cex.sub=0.7, sub=mysub, xlab="AUC GDSC", ylab="AUC CCLE")
  abline(v=cc, col="grey", lty=2)
  abline(h=cc, col="grey", lty=2)
  ## add diagonal line
  abline(a=0, b=1, col="grey", lty=1, lwd=0.5)
}
## legend
plot.new()
par(mar=c(0.1, 0.1, 0.1, 0.1))
legend("center", legend=c("Both sensitive", "", "Both resistant", "", "GDSC sensitive / CCLE resistant", "", "GDSC resistant / CCLE sensitive"), col=c("#2B83BA", "white", "#D7191C", "white", "#FDAE61", "white", "#ABDDA4"), pch=16, bty="n", cex=1, pt.cex=1.5)
dev.off()


auc.star.consistency <- consis[drugix, , , drop=FALSE]
## save consistency statistics
xt <- auc.star.consistency[ , unlist(consistency.stats), "estimate"]
colnames(xt) <- toupper(colnames(xt))
xtable::print.xtable(xtable::xtable(xt, digits=2), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=file.path(saveres, "auc_star_consistency.tex"), append=FALSE)

########################
## Combined barplot for AUC values
########################

## barplots for consistency
pdf(file.path(saveres, "barplot_auc_consistency.pdf"), height=12, width=12)
par(mfrow=c(length(consistency.stats), max(sapply(consistency.stats, length))), mar=c(6, 4, 3, 0), xaxt="n")

for (i in 1:length(consistency.stats)) {
  for (j in 1:length(consistency.stats[[i]])) {
    ## auc publihed
    xx <- auc.published.consistency[ , consistency.stats[[i]][[j]], "estimate"]
    pp <- auc.published.consistency[ , consistency.stats[[i]][[j]], "p"]
    names(xx) <- names(pp) <- rownames(auc.published.consistency)
    pp[xx < 0] <- 1
    xx[xx < 0] <- 0
    res1 <- list("estimate"=xx, "p"=pp)
    ## auc recomputed
    xx <- auc.recomputed.consistency[ , consistency.stats[[i]][[j]], "estimate"]
    pp <- auc.recomputed.consistency[ , consistency.stats[[i]][[j]], "p"]
    names(xx) <- names(pp) <- rownames(auc.recomputed.consistency)
    pp[xx < 0] <- 1
    xx[xx < 0] <- 0
    res2 <- list("estimate"=xx, "p"=pp)
    ## auc recomputed star
    xx <- auc.star.consistency[ , consistency.stats[[i]][[j]], "estimate"]
    pp <- auc.star.consistency[ , consistency.stats[[i]][[j]], "p"]
    names(xx) <- names(pp) <- rownames(auc.recomputed.consistency)
    pp[xx < 0] <- 1
    xx[xx < 0] <- 0
    res3 <- list("estimate"=xx, "p"=pp)
    ## plot
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(res1$estimate, res2$estimate, res3$estimate), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=3), ylab=toupper(consistency.stats[[i]][[j]]), ylim=yylim, angle=c(180, -45, 0), density=c(100, 20, 30), main=names(consistency.stats[[i]])[j])
    legend("topleft", legend=c("AUC published", "AUC recomputed", "AUC*"), fill=c("black", "black", "black"), density=c(100, 20, 30), bty="n", cex=1)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ] + 1.95, y=res1$estimate, pos=2, labels=ifelse(res1$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
    text(x=mp[2, ] + 2.05, y=res2$estimate, pos=2, labels=ifelse(res2$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
    text(x=mp[3, ] + 2.20, y=res3$estimate, pos=2, labels=ifelse(res3$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
    if (j == 1) {
      mtext(side=3, at=par("usr")[1]-6, text=LETTERS[i], line=2, font=2, cex=0.8)
    }
  }
}
dev.off()


########################
## IC50 as published
########################

ic50.gdsc.common <- summarizeSensitivityProfiles(pSet=common.clarified$GDSC, sensitivity.measure="ic50_published", summary.stat="median")
ic50.ccle.common <- summarizeSensitivityProfiles(pSet=common.clarified$CCLE, sensitivity.measure="ic50_published", summary.stat="median")

## compute consistency for ic50 published
myf <- file.path(saveres, sprintf("ic50_%s_consistency.RData", "published"))
if (!file.exists(myf)) {
  consis <- computeConsistencySensitivity(x=ic50.gdsc.common, y=ic50.ccle.common, type="ic50", drugs=drug.cytotoxic, concentrations=list(range.concentration.gdsc, range.concentration.ccle), cutoff=ic50.cutoff, cutoff.cytotoxic=ic50.cytotoxic.cutoff)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn <- c(consisn, list("IC50.PUBLISHED"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
## scatterplots
grDevices::pdf(file=file.path(saveres, sprintf("ic50_%s_bin_lines.pdf", "published")), height=13, width=12)
par(mfrow=c(4, 4), mar=c(5, 4, 3, 2) + 0.1, cex=0.8, las=1)
for (i in drugix) {
  ## truncate ic50= and ic50=Inf by the min and max concentrations tested in each study
  ## ic50 = 0 when the drug dose-response curve starts below 50% viability
  ## ic50 = Inf when 50% viability is never reached
  ic50.gdsc.common.i <- ic50.gdsc.common[i, ]
  ic50.gdsc.common.i[ic50.gdsc.common.i < range.concentration.gdsc[i, 1]] <- range.concentration.gdsc[i, 1]
  ic50.gdsc.common.i[ic50.gdsc.common.i > range.concentration.gdsc[i, 2]] <- range.concentration.gdsc[i, 2]
  ic50.gdsc.common.i <- - log10(ic50.gdsc.common.i / 1000)
  ic50.ccle.common.i <- ic50.ccle.common[i, ]
  ic50.ccle.common.i[ic50.ccle.common.i < range.concentration.ccle[i, 1]] <- range.concentration.ccle[i, 1]
  ic50.ccle.common.i[ic50.ccle.common.i > range.concentration.ccle[i, 2]] <- range.concentration.ccle[i, 2]
  ic50.ccle.common.i <- - log10(ic50.ccle.common.i / 1000)
  ## binarization
  cc <- ifelse(drug.cytotoxic[i] == 1, ic50.cytotoxic.cutoff, ic50.cutoff)
  ic50.gdsc.common.bin <- factor(ifelse(ic50.gdsc.common.i > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))
  ic50.ccle.common.bin <- factor(ifelse(ic50.ccle.common.i > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))  
  ccix <- complete.cases(ic50.gdsc.common.i, ic50.ccle.common.i)
  mylim <- range(c(ic50.gdsc.common.i[ccix], ic50.ccle.common.i[ccix]), na.rm=TRUE)
  mylim <- c(2, 6.25)
  mycol <- rep("#2B83BA", length(ic50.gdsc.common.bin))
  mycol[complete.cases(ic50.gdsc.common.bin, ic50.ccle.common.bin) & ic50.gdsc.common.bin == "resistant" & ic50.ccle.common.bin == "resistant"] <- gplots::col2hex("#D7191C")
  mycol[complete.cases(ic50.gdsc.common.bin, ic50.ccle.common.bin) & ic50.gdsc.common.bin == "sensitive" & ic50.ccle.common.bin == "resistant"] <- gplots::col2hex("#FDAE61")
  mycol[complete.cases(ic50.gdsc.common.bin, ic50.ccle.common.bin) & ic50.gdsc.common.bin == "resistant" & ic50.ccle.common.bin == "sensitive"] <- gplots::col2hex("#ABDDA4")
  mysub <- NULL
  myScatterPlot(x=ic50.gdsc.common.i, y=ic50.ccle.common.i, method=c("transparent"), transparency=0.50, pch=16, minp=50, col=mycol, xlim=mylim, ylim=mylim, main=i, cex.sub=0.7, sub=mysub, xlab="- log10(IC50) GDSC", ylab="- log10(IC50 CCLE)")
  abline(v=cc, col="grey", lty=2)
  abline(h=cc, col="grey", lty=2)
  ## add diagonal line
  abline(a=0, b=1, col="grey", lty=1, lwd=0.5)
}
## legend
plot.new()
par(mar=c(0.1, 0.1, 0.1, 0.1))
legend("center", legend=c("Both sensitive", "", "Both resistant", "", "GDSC sensitive / CCLE resistant", "", "GDSC resistant / CCLE sensitive"), col=c("#2B83BA", "white", "#D7191C", "white", "#FDAE61", "white", "#ABDDA4"), pch=16, bty="n", cex=1, pt.cex=1.5)
dev.off()

ic50.published.consistency <- consis[drugix, , , drop=FALSE]

## save consistency statistics
xt <- ic50.published.consistency[ , unlist(consistency.stats), "estimate"]
colnames(xt) <- toupper(colnames(xt))
xtable::print.xtable(xtable::xtable(xt, digits=2), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=file.path(saveres, "ic50_published_consistency.tex"), append=FALSE)

########################
## IC50 recomputed
########################

ic50.gdsc.common <- summarizeSensitivityProfiles(pSet=common.clarified$GDSC, sensitivity.measure="ic50_recomputed", summary.stat="median")
ic50.ccle.common <- summarizeSensitivityProfiles(pSet=common.clarified$CCLE, sensitivity.measure="ic50_recomputed", summary.stat="median")

## compute consistency for ic50 recomputed
myf <- file.path(saveres, sprintf("ic50_%s_consistency.RData", "recomputed"))
if (!file.exists(myf)) {
  consis <- computeConsistencySensitivity(x=ic50.gdsc.common, y=ic50.ccle.common, type="ic50", drugs=drug.cytotoxic, concentrations=list(range.concentration.gdsc, range.concentration.ccle), cutoff=ic50.cutoff, cutoff.cytotoxic=ic50.cytotoxic.cutoff)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn <- c(consisn, list("IC50.RECOMPUTED"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
## scatterplots
grDevices::pdf(file=file.path(saveres, sprintf("ic50_%s_bin_lines.pdf", "recomputed")), height=13, width=12)
par(mfrow=c(4, 4), mar=c(5, 4, 3, 2) + 0.1, cex=0.8, las=1)
for (i in drugix) {
  ## truncate ic50= and ic50=Inf by the min and max concentrations tested in each study
  ## ic50 = 0 when the drug dose-response curve starts below 50% viability
  ## ic50 = Inf when 50% viability is never reached
  ic50.gdsc.common.i <- ic50.gdsc.common[i, ]
  ic50.gdsc.common.i[ic50.gdsc.common.i < range.concentration.gdsc[i, 1]] <- range.concentration.gdsc[i, 1]
  ic50.gdsc.common.i[ic50.gdsc.common.i > range.concentration.gdsc[i, 2]] <- range.concentration.gdsc[i, 2]
  ic50.gdsc.common.i <- - log10(ic50.gdsc.common.i / 1000)
  ic50.ccle.common.i <- ic50.ccle.common[i, ]
  ic50.ccle.common.i[ic50.ccle.common.i < range.concentration.ccle[i, 1]] <- range.concentration.ccle[i, 1]
  ic50.ccle.common.i[ic50.ccle.common.i > range.concentration.ccle[i, 2]] <- range.concentration.ccle[i, 2]
  ic50.ccle.common.i <- - log10(ic50.ccle.common.i / 1000)
  ## binarization
  cc <- ifelse(drug.cytotoxic[i] == 1, ic50.cytotoxic.cutoff, ic50.cutoff)
  ic50.gdsc.common.bin <- factor(ifelse(ic50.gdsc.common.i > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))
  ic50.ccle.common.bin <- factor(ifelse(ic50.ccle.common.i > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))  
  ccix <- complete.cases(ic50.gdsc.common.i, ic50.ccle.common.i)
  mylim <- range(c(ic50.gdsc.common.i[ccix], ic50.ccle.common.i[ccix]), na.rm=TRUE)
  mylim <- c(2, 6.25)
  mycol <- rep("#2B83BA", length(ic50.gdsc.common.bin))
  mycol[complete.cases(ic50.gdsc.common.bin, ic50.ccle.common.bin) & ic50.gdsc.common.bin == "resistant" & ic50.ccle.common.bin == "resistant"] <- gplots::col2hex("#D7191C")
  mycol[complete.cases(ic50.gdsc.common.bin, ic50.ccle.common.bin) & ic50.gdsc.common.bin == "sensitive" & ic50.ccle.common.bin == "resistant"] <- gplots::col2hex("#FDAE61")
  mycol[complete.cases(ic50.gdsc.common.bin, ic50.ccle.common.bin) & ic50.gdsc.common.bin == "resistant" & ic50.ccle.common.bin == "sensitive"] <- gplots::col2hex("#ABDDA4")
  mysub <- NULL
  myScatterPlot(x=ic50.gdsc.common.i, y=ic50.ccle.common.i, method=c("transparent"), transparency=0.50, pch=16, minp=50, col=mycol, xlim=mylim, ylim=mylim, main=i, cex.sub=0.7, sub=mysub, xlab="- log10(IC50) GDSC", ylab="- log10(IC50 CCLE)")
  # plot(x=ic50.gdsc.common.i, y=ic50.ccle.common.i, method=c("transparent"), transparency=0.50, pch=16, minp=50, col=mycol, xlim=mylim, ylim=mylim, main=i, cex.sub=0.7, sub=mysub, xlab="IC50 GDSC", ylab="IC50 CCLE")
  abline(v=cc, col="grey", lty=2)
  abline(h=cc, col="grey", lty=2)
  ## add diagonal line
  abline(a=0, b=1, col="grey", lty=1, lwd=0.5)
}
## legend
plot.new()
par(mar=c(0.1, 0.1, 0.1, 0.1))
legend("center", legend=c("Both sensitive", "", "Both resistant", "", "GDSC sensitive / CCLE resistant", "", "GDSC resistant / CCLE sensitive"), col=c("#2B83BA", "white", "#D7191C", "white", "#FDAE61", "white", "#ABDDA4"), pch=16, bty="n", cex=1, pt.cex=1.5)
dev.off()

ic50.recomputed.consistency <- consis[drugix, , , drop=FALSE]

## save consistency statistics
xt <- ic50.recomputed.consistency[ , unlist(consistency.stats), "estimate"]
colnames(xt) <- toupper(colnames(xt))
xtable::print.xtable(xtable::xtable(xt, digits=2), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=file.path(saveres, "ic50_recomputed_consistency.tex"), append=FALSE)

########################
## Combined barplot for IC50 recomputed
########################

## barplots for consistency
pdf(file.path(saveres, "barplot_ic50_consistency.pdf"), height=12, width=12)
par(mfrow=c(length(consistency.stats), max(sapply(consistency.stats, length))), mar=c(6, 4, 3, 0), xaxt="n")
for (i in 1:length(consistency.stats)) {
  for (j in 1:length(consistency.stats[[i]])) {
    ## ic50 published
    xx <- ic50.published.consistency[ , consistency.stats[[i]][[j]], "estimate"]
    pp <- ic50.published.consistency[ , consistency.stats[[i]][[j]], "p"]
    names(xx) <- names(pp) <- rownames(ic50.published.consistency)
    pp[xx < 0] <- 1
    xx[xx < 0] <- 0
    res1 <- list("estimate"=xx, "p"=pp)
    ## ic50 recomputed
    xx <- ic50.recomputed.consistency[ , consistency.stats[[i]][[j]], "estimate"]
    pp <- ic50.recomputed.consistency[ , consistency.stats[[i]][[j]], "p"]
    names(xx) <- names(pp) <- rownames(ic50.recomputed.consistency)
    pp[xx < 0] <- 1
    xx[xx < 0] <- 0
    res2 <- list("estimate"=xx, "p"=pp)
    yylim <- c(0, 1)
    mp <- barplot(height=rbind(res1$estimate, res2$estimate), beside=TRUE, space=c(0.1, 2), col=rep(rainbow(length(xx), v=0.9), each=2), ylab=toupper(consistency.stats[[i]][[j]]), ylim=yylim, angle=c(45, -45), density=c(100, 40), main=names(consistency.stats[[i]])[j])
    legend("topleft", legend=c("IC50 published", "IC50 recomputed"), fill=c("black", "black"), density=c(100, 40), bty="n", cex=1)
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
    text(x=mp[1, ] + 1.55, y=rbind(res1$estimate, res2$estimate)[1, ], pos=2, labels=ifelse(res1$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
    text(x=mp[2, ] + 1.55, y=rbind(res1$estimate, res2$estimate)[2, ], pos=2, labels=ifelse(res2$p < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
    if (j == 1) {
      mtext(side=3, at=par("usr")[1]-6, text=LETTERS[i], line=2, font=2, cex=0.8)
    }
  }
}
dev.off()


WriteXLS::WriteXLS(x="consisn", row.names=TRUE, ExcelFileName=file.path(saveres, "consis_auc_ic50.xlsx"))

################################################
## consistency of gene expression
################################################

consisn2 <- NULL

########################
## ccle affy vs gdsc affy
########################

myf <- file.path(saveres, sprintf("rna_common.RData"))
if (!file.exists(myf)) {
  rna.gdsc.common <- exprs(summarizeMolecularProfiles(pSet=common.clarified$GDSC, mDataType="rna2", summary.stat="median"))
  rna.ccle.common <- exprs(summarizeMolecularProfiles(pSet=common.clarified$CCLE, mDataType="rna", summary.stat="median"))
  ## common genes
  gix <- intersect(rownames(rna.gdsc.common), rownames(rna.ccle.common))
  rna.gdsc.common <- rna.gdsc.common[gix, ]
  rna.ccle.common <- rna.ccle.common[gix, ]
  save(list=c("rna.gdsc.common", "rna.ccle.common", "gix"), compress=TRUE, file=myf)
} else {
  load(myf)
}

## compute consistency
myf <- file.path(saveres, sprintf("rna_consistency.RData"))
if (!file.exists(myf)) {
  mystats <- c("pcc.full", "pcc.high", "scc.full", "scc.high", "dxy.full", "dxy.high", "cosine", "mcc", "kappa", "cramerv", "inform")
  consis <- array(NA, dim=c(length(gix), length(mystats), 2), dimnames=list(gix, mystats, c("estimate", "p")))
  pb <- utils::txtProgressBar(min=0, max=length(gix), style=3)
  ii <- 1
  for (i in gix) {
    ## binarization
    ## consistency measures
    iix <- !(is.na(rna.gdsc.common.bin) | is.na(rna.ccle.common.bin)) & !(rna.gdsc.common.bin == "low" & rna.ccle.common.bin == "low")
    ccix <- complete.cases(rna.gdsc.common[i, ], rna.ccle.common[i, ])
    ## pcc full
    res <- try(cor.test(rna.gdsc.common[i, ], rna.ccle.common[i, ], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.full", ] <- c(res$estimate, res$p.value)
    }
    ## pcc sensitive
    res <- try(cor.test(rna.gdsc.common[i, iix], rna.ccle.common[i, iix], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.high", ] <- c(res$estimate, res$p.value)
    }
    ## scc full
    res <- try(cor.test(rna.gdsc.common[i, ], rna.ccle.common[i, ], method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.full", ] <- c(res$estimate, res$p.value)
    }
    ## scc sensitive
    res <- try(cor.test(rna.gdsc.common[i, iix], rna.ccle.common[i, iix], method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.high", ] <- c(res$estimate, res$p.value)
    }
    ## dxy
    res <- try(mRMRe::correlate(X=rna.gdsc.common[i, ], Y=rna.ccle.common[i, ], method="cindex", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "dxy.full", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## dxy sensitive
    res <- try(mRMRe::correlate(X=rna.gdsc.common[i, iix], Y=rna.ccle.common[i, iix], method="cindex", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "dxy.high", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## cosine
    res <- try(PharmacoGx::cosinePerm(x=rna.gdsc.common[i, ccix] - cutoff.rna.gdsc, y=rna.ccle.common[i, ccix] - cutoff.rna.ccle, nperm=0, nthread=nbcore, alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "cosine", ] <- c(res$estimate, res$p.value)
    }
    ## mcc
    res <- try(PharmacoGx::mcc(x=rna.gdsc.common.bin, y=rna.ccle.common.bin, nperm=0, alternative="greater", nthread=nbcore), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "mcc", ] <- c(res$estimate, res$p.value)
    }
    ## kappa
    tt <- table(x=rna.gdsc.common.bin, y=rna.ccle.common.bin)
    res <- try(epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "kappa", "estimate"] <- res$kappa
    }
    ## cramerv
    res <- try(vcd::assocstats(x=tt))
    if (class(res) != "try-error") {
      consis[i, "cramerv", ] <- c(res$cramer, res$chisq_tests[1, 3])
    }
    ## Informedness: 2 * balanced accuracy - 1
    res <- try(caret::sensitivity(data=rna.gdsc.common.bin, reference=rna.ccle.common.bin) + caret::specificity(data=rna.gdsc.common.bin, reference=rna.ccle.common.bin) - 1)
    if (class(res) != "try-error") {
      consis[i, "inform", "estimate"] <- res
    }
    utils::setTxtProgressBar(pb, ii)
    ii <- ii + 1
  }
  close(pb)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn2 <- c(consisn2, list("GE.ARRAYS"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
rna.consistency <- consis

########################
## ccle affy vs ccle rnaseq
########################

myf <- file.path(saveres, sprintf("rnaseq_common.RData"))
if (!file.exists(myf)) {
  ## gene expression data from CCLE
  rna.ccle.common <- exprs(summarizeMolecularProfiles(pSet=CCLE.clarified, mDataType="rna", summary.stat="median"))
  rnaseq.ccle.common <- exprs(summarizeMolecularProfiles(pSet=CCLE.clarified, mDataType="rnaseq", summary.stat="median"))
  ## common genes
  gix <- intersect(rownames(rna.ccle.common), rownames(rnaseq.ccle.common))
  rna.ccle.common <- rna.ccle.common[gix, ]
  rnaseq.ccle.common <- rnaseq.ccle.common[gix, ]
  save(list=c("rna.ccle.common", "rnaseq.ccle.common", "gix"), compress=TRUE, file=myf)
} else {
  load(myf)
}

## compute consistency
myf <- file.path(saveres, sprintf("rnaseq_consistency.RData"))
if (!file.exists(myf)) {
  mystats <- c("pcc.full", "pcc.high", "scc.full", "scc.high", "dxy.full", "dxy.high", "cosine", "mcc", "kappa", "cramerv", "inform")
  consis <- array(NA, dim=c(length(gix), length(mystats), 2), dimnames=list(gix, mystats, c("estimate", "p")))
  pb <- utils::txtProgressBar(min=0, max=length(gix), style=3)
  ii <- 1
  for (i in gix) {
    ## binarization
    cc <- cutoff.rnaseq.ccle
    rnaseq.ccle.common.bin <- factor(ifelse(rnaseq.ccle.common[i, ] > cc, "high", "low"), levels=c("low","high"))
    cc <- cutoff.rna.ccle
    rna.ccle.common.bin <- factor(ifelse(rna.ccle.common[i, ] > cc, "high", "low"), levels=c("low","high"))
    cc <- cutoff.rna.gdsc
    rna.gdsc.common.bin <- factor(ifelse(rna.gdsc.common[i, ] > cc, "high", "low"), levels=c("low","high"))
    ## consistency measures
    iix <- !(is.na(rnaseq.ccle.common.bin) | is.na(rna.ccle.common.bin)) & !(rnaseq.ccle.common.bin == "low" & rna.ccle.common.bin == "low")
    ccix <- complete.cases(rnaseq.ccle.common[i, ], rna.ccle.common[i, ])
    ## pcc full
    res <- try(cor.test(rnaseq.ccle.common[i, ], rna.ccle.common[i, ], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.full", ] <- c(res$estimate, res$p.value)
    }
    ## pcc high
    res <- try(cor.test(rnaseq.ccle.common[i, iix], rna.ccle.common[i, iix], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.high", ] <- c(res$estimate, res$p.value)
    }
    ## scc full
    res <- try(cor.test(rnaseq.ccle.common[i, ], rna.ccle.common[i, ], method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.full", ] <- c(res$estimate, res$p.value)
    }
    ## scc high
    res <- try(cor.test(rnaseq.ccle.common[i, iix], rna.ccle.common[i, iix], method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.high", ] <- c(res$estimate, res$p.value)
    }
    ## dxy
    res <- try(mRMRe::correlate(X=rnaseq.ccle.common[i, ], Y=rna.ccle.common[i, ], method="cindex", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "dxy.full", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## dxy high
    res <- try(mRMRe::correlate(X=rnaseq.ccle.common[i, iix], Y=rna.ccle.common[i, iix], method="cindex", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "dxy.high", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## cosine
    res <- try(PharmacoGx::cosinePerm(x=rnaseq.ccle.common[i, ccix] - cutoff.rnaseq.ccle, y=rna.ccle.common[i, ccix] - cutoff.rna.ccle, nperm=0, nthread=nbcore, alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "cosine", ] <- c(res$estimate, res$p.value)
    }
    ## mcc
    res <- try(PharmacoGx::mcc(x=rnaseq.ccle.common.bin, y=rna.ccle.common.bin, nperm=0, nthread=nbcore), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "mcc", ] <- c(res$estimate, res$p.value)
    }
    ## kappa
    tt <- table(x=rnaseq.ccle.common.bin, y=rna.ccle.common.bin)
    res <- try(epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "kappa", "estimate"] <- res$kappa
    }
    ## cramerv
    res <- try(vcd::assocstats(x=tt))
    if (class(res) != "try-error") {
      consis[i, "cramerv", ] <- c(res$cramer, res$chisq_tests[1, 3])
    }
    ## Informedness: 2 * balanced accuracy - 1
    res <- try(caret::sensitivity(data=rnaseq.ccle.common.bin, reference=rna.ccle.common.bin) + caret::specificity(data=rnaseq.ccle.common.bin, reference=rna.ccle.common.bin) - 1)
    if (class(res) != "try-error") {
      consis[i, "inform", "estimate"] <- res
    }
    utils::setTxtProgressBar(pb, ii)
    ii <- ii + 1
  }
  close(pb)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn2 <- c(consisn2, list("GE.CCLE.ARRAY.RNASEQ"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
rnaseq.consistency <- consis

########################
## ccle rnaseq vs gdsc affy
########################

myf <- file.path(saveres, sprintf("rna_rnaseq_common.RData"))
if (!file.exists(myf)) {
  rna.gdsc.common <- exprs(summarizeMolecularProfiles(pSet=common.clarified$GDSC, mDataType="rna2", summary.stat="median"))
  rnaseq.ccle.common <- exprs(summarizeMolecularProfiles(pSet=common.clarified$CCLE, mDataType="rnaseq", summary.stat="median"))
  ## common genes
  gix <- intersect(rownames(rna.gdsc.common), rownames(rnaseq.ccle.common))
  rna.gdsc.common <- rna.gdsc.common[gix, ]
  rnaseq.ccle.common <- rnaseq.ccle.common[gix, ]
  save(list=c("rna.gdsc.common", "rnaseq.ccle.common", "gix"), compress=TRUE, file=myf)
} else {
  load(myf)
}

## compute consistency
myf <- file.path(saveres, sprintf("rna_rnaseq_consistency.RData"))
if (!file.exists(myf)) {
  mystats <- c("pcc.full", "pcc.high", "scc.full", "scc.high", "dxy.full", "dxy.high", "cosine", "mcc", "kappa", "cramerv", "inform")
  consis <- array(NA, dim=c(length(gix), length(mystats), 2), dimnames=list(gix, mystats, c("estimate", "p")))
  pb <- utils::txtProgressBar(min=0, max=length(gix), style=3)
  ii <- 1
  for (i in gix) {
    ## binarization
    cc <- cutoff.rna.gdsc
    rna.gdsc.common.bin <- factor(ifelse(rna.gdsc.common[i, ] > cc, "high", "low"), levels=c("low","high"))
    cc <- cutoff.rnaseq.ccle
    rnaseq.ccle.common.bin <- factor(ifelse(rnaseq.ccle.common[i, ] > cc, "high", "low"), levels=c("low","high"))
    ## consistency measures
    iix <- !(is.na(rna.gdsc.common.bin) | is.na(rnaseq.ccle.common.bin)) & !(rna.gdsc.common.bin == "low" & rnaseq.ccle.common.bin == "low")
    ccix <- complete.cases(rna.gdsc.common[i, ], rnaseq.ccle.common[i, ])
    ## pcc full
    res <- try(cor.test(rna.gdsc.common[i, ], rnaseq.ccle.common[i, ], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.full", ] <- c(res$estimate, res$p.value)
    }
    ## pcc sensitive
    res <- try(cor.test(rna.gdsc.common[i, iix], rnaseq.ccle.common[i, iix], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.high", ] <- c(res$estimate, res$p.value)
    }
    ## scc full
    res <- try(cor.test(rna.gdsc.common[i, ], rnaseq.ccle.common[i, ], method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.full", ] <- c(res$estimate, res$p.value)
    }
    ## scc sensitive
    res <- try(cor.test(rna.gdsc.common[i, iix], rnaseq.ccle.common[i, iix], method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.high", ] <- c(res$estimate, res$p.value)
    }
    ## dxy
    res <- try(mRMRe::correlate(X=rna.gdsc.common[i, ], Y=rnaseq.ccle.common[i, ], method="cindex", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "dxy.full", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## dxy sensitive
    res <- try(mRMRe::correlate(X=rna.gdsc.common[i, iix], Y=rnaseq.ccle.common[i, iix], method="cindex", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "dxy.high", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## cosine
    res <- try(PharmacoGx::cosinePerm(x=rna.gdsc.common[i, ccix] - cutoff.rna.gdsc, y=rnaseq.ccle.common[i, ccix] - cutoff.rnaseq.ccle, nperm=0, nthread=nbcore, alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "cosine", ] <- c(res$estimate, res$p.value)
    }
    ## mcc
    res <- try(PharmacoGx::mcc(x=rna.gdsc.common.bin, y=rnaseq.ccle.common.bin, nperm=0, alternative="greater", nthread=nbcore), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "mcc", ] <- c(res$estimate, res$p.value)
    }
    ## kappa
    tt <- table(x=rna.gdsc.common.bin, y=rnaseq.ccle.common.bin)
    res <- try(epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "kappa", "estimate"] <- res$kappa
    }
    ## cramerv
    res <- try(vcd::assocstats(x=tt))
    if (class(res) != "try-error") {
      consis[i, "cramerv", ] <- c(res$cramer, res$chisq_tests[1, 3])
    }
    ## Informedness: 2 * balanced accuracy - 1
    res <- try(caret::sensitivity(data=rna.gdsc.common.bin, reference=rnaseq.ccle.common.bin) + caret::specificity(data=rna.gdsc.common.bin, reference=rnaseq.ccle.common.bin) - 1)
    if (class(res) != "try-error") {
      consis[i, "inform", "estimate"] <- res
    }
    utils::setTxtProgressBar(pb, ii)
    ii <- ii + 1
  }
  close(pb)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn2 <- c(consisn2, list("GE.ARRAY.RNASEQ"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
rna.rnaseq.consistency <- consis


################################################
## consistency of cnv profiles
################################################
  
myf <- file.path(saveres, sprintf("cnv_common.RData"))
if (!file.exists(myf)) {
  cnv.gdsc.common <- exprs(summarizeMolecularProfiles(pSet=common.clarified$GDSC, mDataType="cnv", summary.stat="median"))
  cnv.ccle.common <- exprs(summarizeMolecularProfiles(pSet=common.clarified$CCLE, mDataType="cnv", summary.stat="median"))
  ## common genes
  gix <- intersect(rownames(cnv.gdsc.common), rownames(cnv.ccle.common))
  cnv.gdsc.common <- cnv.gdsc.common[gix, ]
  cnv.ccle.common <- cnv.ccle.common[gix, ]
  save(list=c("cnv.gdsc.common", "cnv.ccle.common", "gix"), compress=TRUE, file=myf)
} else {
  load(myf)
}

## compute consistency
myf <- file.path(saveres, sprintf("cnv_consistency.RData"))
if (!file.exists(myf)) {
  mystats <- c("pcc.full", "pcc.high", "scc.full", "scc.high", "dxy.full", "dxy.high", "cosine", "mcc", "kappa", "cramerv", "inform")
  consis <- array(NA, dim=c(length(gix), length(mystats), 2), dimnames=list(gix, mystats, c("estimate", "p")))
  pb <- utils::txtProgressBar(min=0, max=length(gix), style=3)
  ii <- 1
  for (i in gix) {
    ## binarization
    cc1 <- cutoff.cnv.amplification
    cc2 <- cutoff.cnv.deletion
    
    cnv.gdsc.common.bin <- factor(ifelse(cnv.gdsc.common[i, ] > cc1 | cnv.gdsc.common[i, ] < cc2, "variation", "normal"), levels=c("normal","variation"))
    cnv.ccle.common.bin <- factor(ifelse(cnv.ccle.common[i, ] > cc1 | cnv.ccle.common[i, ] < cc2 , "variation", "normal"), levels=c("normal","variation"))
    ## consistency measures
    iix <- !(is.na(cnv.gdsc.common.bin) | is.na(cnv.ccle.common.bin)) & !(cnv.gdsc.common.bin == "normal" & cnv.ccle.common.bin == "normal")
    ccix <- complete.cases(cnv.gdsc.common[i, ], cnv.ccle.common[i, ])
    ## pcc full
    res <- try(cor.test(cnv.gdsc.common[i, ], cnv.ccle.common[i, ], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.full", ] <- c(res$estimate, res$p.value)
    }
    ## pcc sensitive
    res <- try(cor.test(cnv.gdsc.common[i, iix], cnv.ccle.common[i, iix], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.high", ] <- c(res$estimate, res$p.value)
    }
    ## scc full
    res <- try(cor.test(cnv.gdsc.common[i, ], cnv.ccle.common[i, ], method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.full", ] <- c(res$estimate, res$p.value)
    }
    ## scc sensitive
    res <- try(cor.test(cnv.gdsc.common[i, iix], cnv.ccle.common[i, iix], method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.high", ] <- c(res$estimate, res$p.value)
    }
    ## dxy
    res <- try(mRMRe::correlate(X=cnv.gdsc.common[i, ], Y=cnv.ccle.common[i, ], method="cindex", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "dxy.full", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## dxy sensitive
    res <- try(mRMRe::correlate(X=cnv.gdsc.common[i, iix], Y=cnv.ccle.common[i, iix], method="cindex", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "dxy.high", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## cosine
    x <- cnv.gdsc.common[i, ccix]
    x[x > cc2] <- x[x > cc2] - cc1
    
    y <- cnv.ccle.common[i, ccix]
    y[y > cc2] <- y[y > cc2] - cc1
    
    
    res <- try(PharmacoGx::cosinePerm(x=x, y=y, nperm=0, nthread=nbcore, alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "cosine", ] <- c(res$estimate, res$p.value)
    }
    ## mcc
    res <- try(PharmacoGx::mcc(x=cnv.gdsc.common.bin, y=cnv.ccle.common.bin, nperm=0, alternative="greater", nthread=nbcore), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "mcc", ] <- c(res$estimate, res$p.value)
    }
    ## kappa
    tt <- table(x=cnv.gdsc.common.bin, y=cnv.ccle.common.bin)
    res <- try(epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "kappa", "estimate"] <- res$kappa
    }
    ## cramerv
    res <- try(vcd::assocstats(x=tt))
    if (class(res) != "try-error") {
      consis[i, "cramerv", ] <- c(res$cramer, res$chisq_tests[1, 3])
    }
    ## Informedness: 2 * balanced accuracy - 1
    res <- try(caret::sensitivity(data=cnv.gdsc.common.bin, reference=cnv.ccle.common.bin) + caret::specificity(data=cnv.gdsc.common.bin, reference=cnv.ccle.common.bin) - 1)
    if (class(res) != "try-error") {
      consis[i, "inform", "estimate"] <- res
    }
    utils::setTxtProgressBar(pb, ii)
    ii <- ii + 1
  }
  close(pb)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn2 <- c(consisn2, list("CNV"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
cnv.consistency <- consis

################################################
## consistency of mutation profiles
################################################
  
myf <- file.path(saveres, sprintf("mutation_common.RData"))
if (!file.exists(myf)) {
  mut.gdsc.common <- exprs(summarizeMolecularProfiles(pSet=common.clarified$GDSC, mDataType="mutation", summary.stat="or"))
  nn <- dimnames(mut.gdsc.common)
  mut.gdsc.common <- apply(mut.gdsc.common, 2, as.numeric)
  dimnames(mut.gdsc.common) <- nn
  mut.ccle.common <- exprs(summarizeMolecularProfiles(pSet=common.clarified$CCLE, mDataType="mutation", summary.stat="or"))
  nn <- dimnames(mut.ccle.common)
  mut.ccle.common <- apply(mut.ccle.common, 2, as.numeric)
  dimnames(mut.ccle.common) <- nn
  ## common genes
  gix <- intersect(rownames(mut.gdsc.common), rownames(mut.ccle.common))
  mut.gdsc.common <- mut.gdsc.common[gix, ]
  mut.ccle.common <- mut.ccle.common[gix, ]
  save(list=c("mut.gdsc.common", "mut.ccle.common", "gix"), compress=TRUE, file=myf)
} else {
  load(myf)
}

## compute consistency
myf <- file.path(saveres, sprintf("mutation_consistency.RData"))
if (!file.exists(myf)) {
  mystats <- c("mcc", "kappa", "cramerv", "inform")
  consis <- array(NA, dim=c(length(gix), length(mystats), 2), dimnames=list(gix, mystats, c("estimate", "p")))
  pb <- utils::txtProgressBar(min=0, max=length(gix), style=3)
  ii <- 1
  for (i in gix) {

    ## mcc
    res <- try(PharmacoGx::mcc(x=factor(mut.gdsc.common[i, ], levels=c(0, 1)), y=factor(mut.ccle.common[i, ], levels=c(0, 1)), nperm=0, nthread=nbcore), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "mcc", ] <- c(res$estimate, res$p.value)
    }
    ## kappa
    tt <- table(x=mut.gdsc.common[i, ], y=mut.ccle.common[i, ])
    res <- try(epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "kappa", "estimate"] <- res$kappa
    }
    ## cramerv
    res <- try(vcd::assocstats(x=tt))
    if (class(res) != "try-error") {
      consis[i, "cramerv", ] <- c(res$cramer, res$chisq_tests[1, 3])
    }
    ## Informedness: 2 * balanced accuracy - 1
    res <- try(caret::sensitivity(data=factor(mut.gdsc.common[i, ], levels=c(0, 1)), reference=factor(mut.ccle.common[i, ], levels=c(0, 1))) + caret::specificity(data=factor(mut.gdsc.common[i, ], levels=c(0, 1)), reference=factor(mut.ccle.common[i, ], levels=c(0, 1))) - 1)
    if (class(res) != "try-error") {
      consis[i, "inform", "estimate"] <- res
    }
    utils::setTxtProgressBar(pb, ii)
    ii <- ii + 1
  }
  close(pb)
  save(list=c("consis"), compress=TRUE, file=myf)
} else {
  load(myf)
}
consisn2 <- c(consisn2, list("MUTATION"=data.frame(consis[ , , "estimate"], check.names=FALSE)))
mut.consistency <- consis

WriteXLS::WriteXLS(x="consisn2", row.names=TRUE, ExcelFileName=file.path(saveres, "consis_molecular.xlsx"))

########################
## boxplot for consistency of gene expression
########################

pdf(file.path(saveres, "boxplot_rna_sensitivity_consistency.pdf"), height=12, width=12)
par(mfrow=c(length(consistency.stats), max(sapply(consistency.stats, length))), mar=c(7, 5, 3, 1), xaxt="n")
wp.all <- nn.all <- nn2.all <- NULL
for (i in 1:length(consistency.stats)) {
  for (j in 1:length(consistency.stats[[i]])) {
    nn <- consistency.stats[[i]][[j]]
    nn2 <- names(consistency.stats[[i]])[j]
    xx <- list(
      "ge.ccle.array.rnaseq"=rnaseq.consistency[ , gsub("sens", "high", nn), "estimate"],
      "ge.arrays"=rna.consistency[ , gsub("sens", "high", nn), "estimate"],
      "ge.array.rnaseq"=rna.rnaseq.consistency[ , gsub("sens", "high", nn), "estimate"],
      "cnv"=cnv.consistency[ , gsub("sens", "high", nn), "estimate"],
      " "=NA,
      "auc.published"=auc.published.consistency[ , nn, "estimate"],
      "auc.recomputed"=auc.recomputed.consistency[ , nn, "estimate"],
      "auc.star"=auc.star.consistency[ , nn, "estimate"],
      " "=NA,
      "ic50.published"=ic50.published.consistency[ , nn, "estimate"],
      "ic50.recomputed"=ic50.recomputed.consistency[ , nn, "estimate"]
      )
    yylim <- c(-1, 1)
    bp <- boxplot(xx, ylim=yylim, col="lightgrey", border="black", outline=FALSE, pars=list(outcol="black", outpch=20, outcex=0.5), main=gsub("sensitive", "high expression/sensitive", names(consistency.stats[[i]])[j]), ylab=toupper(unlist(strsplit(consistency.stats[[i]][[j]], "[.]"))[1]))
    axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
    text(x=1:length(xx), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2, cex=0.8)
    if (j == 1) {
      mtext(side=3, at=par("usr")[1]-2, text=LETTERS[i], line=2, font=2, cex=0.8)
    }
    
    ## pairwise wilcoxon test
    xx2 <- xx[!is.na(xx)]
    cc <- combn(x=length(xx2), m=2)
    rr <- apply(cc, 2, function (x, y) {
      rr <- wilcox.test(x=y[[x[1]]], y=y[[x[2]]], exact=FALSE)
      return (rr$p.value)
    }, y=xx2)
    wp <- wp3 <- matrix(NA, nrow=length(xx2), ncol=length(xx2), dimnames=list(names(xx2), names(xx2)))
    wp[t(cc)] <- rr
    wp <- wp[!apply(is.na(wp), 1, all), !apply(is.na(wp), 2, all)]
    wp.all <- c(wp.all, list(wp))
    nn.all <- c(nn.all, nn)
    nn2.all <- c(nn2.all, nn2)
  }
}
names(wp.all) <- nn.all
dev.off()

pdf(file.path(saveres, "heatmap_pvalue_rna_sensitivity_consistency.pdf"), height=12, width=12)
par(mfrow=c(length(consistency.stats), max(sapply(consistency.stats, length))), mar=c(9, 8, 3, 1), xaxt="n", yaxt="n")
pp <- NULL
j <- 1
for (i in 1:length(wp.all)) {
  wp <- wp.all[[i]]
  wp2 <- apply(wp, c(1,2), function (x) {
    return (sprintf("%.1E", x))
  })
  pp <- c(pp, list(data.frame(wp2)))
  ## plot image
  wp3 <- wp
  wp3[] <- NA
  wp3[!is.na(wp) & wp < 0.05] <- 1
  wp3[!is.na(wp) & wp >= 0.05 & wp < 0.10] <- 2
  wp3[!is.na(wp) & wp >= 0.10 & wp < 0.50] <- 3
  wp3[!is.na(wp) & wp >= 0.50] <- 4

  mycol <- c("#d73027", "#fee090", "#e0f3f8")
  image(x=1:nrow(wp3), y=1:ncol(wp3), z=wp3[ , ncol(wp3):1], col=mycol, axes=FALSE, xlab="", ylab="", main=gsub("sensitive", "high expression/sensitive", nn2.all[i]), cex.main=1.1)
  # axis(3, at=1:ncol(wp3), labels=colnames(wp3), srt=45, xpd=NA, font=2, cex=0.8, tick=FALSE)
  # axis(2, at=1:nrow(wp3), labels=rownames(wp3), srt=45, xpd=NA, font=2, cex=0.8, tick=FALSE)
  text(x=par("usr")[1] + 1:nrow(wp3) - 0.5, y=par("usr")[3] - (par("usr")[4] * 0.04), pos=2, labels=toupper(rownames(wp3)), srt=90, xpd=NA, font=2, cex=0.8)
  text(x=par("usr")[1], y=par("usr")[3] + 1:ncol(wp3) - 0.3, pos=2, labels=rev(toupper(colnames(wp3))), srt=0, xpd=NA, font=2, cex=0.8)
  grid(nx=nrow(wp3), ny=ncol(wp3), col="gray", lty=1, lwd=0.8)
  legend("topright", legend=c("p < 0.05", "0.05 <= p < 0.10", "p >= 0.10"), col=mycol, pch=15, pt.cex=1.75, box.lwd=0, box.col="white", bg="white")
  if ((i %% max(sapply(consistency.stats, length))) == 1) {
    message(LETTERS[j])
    mtext(side=3, at=par("usr")[1]-2, text=LETTERS[j], line=2, font=2, cex=0.8)
    j <- j + 1
  }
}
names(pp) <- names(wp.all)
dev.off()

## save consistency statistics
WriteXLS::WriteXLS(x="pp", ExcelFileName=file.path(saveres, "significance_consistency.xls"))



#################################################
###biomarkers consistancy

load(file.path(saveres, "PSets","signatures_data.RData"))
load(file.path(saveres, "PSets", "GDSC_clarified.RData"))
load(file.path(saveres, "PSets", "CCLE_clarified.RData"))
load(file.path(saveres, "PSets", "common_clarified.RData"))

myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_rna2.RData")
if(!file.exists(myf))
{
  gdsc.sig.rna2 <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="rna2", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  
  save(gdsc.sig.rna2, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_rna2_binary.RData")
if(!file.exists(myf))
{
  gdsc.sig.rna2.bin <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="rna2", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  
  save(gdsc.sig.rna2.bin, file=myf)
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

myf <- file.path(saveres, "PSets", "Sigs", "common_gdsc_sig_rna2_binary.RData")
if(!file.exists(myf))
{
  common.gdsc.sig.rna2.bin <- drugSensitivitySig(pSet=common.clarified$GDSC, mDataType="rna2", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", sensitivity.cutoff=0.2, nthread=4)
  save(common.gdsc.sig.rna2.bin, file=myf)
  
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "ccle_sig_rna.RData")
if(!file.exists(myf))
{
  ccle.sig.rna <- drugSensitivitySig(pSet=CCLE.clarified, mDataType="rna", drugs=drugs, features=features, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(ccle.sig.rna, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "ccle_sig_rna_binary.RData")
if(!file.exists(myf))
{
  ccle.sig.rna.bin <- drugSensitivitySig(pSet=CCLE.clarified, mDataType="rna", drugs=drugs, features=features[1:20], sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", sensitivity.cutoff=0.2, nthread=4)
  save(ccle.sig.rna.bin, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "common_ccle_sig_rna.RData")
if(!file.exists(myf))
{
  common.ccle.sig.rna <- drugSensitivitySig(pSet=common.clarified$CCLE, mDataType="rna", drugs=drugs, features=features,sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(common.ccle.sig.rna, file=myf)
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

continuous.all <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=continuous.all, method="continuous", cell="all")
continuous.all.100 <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous", drugs=drugs, fdr.cut.off=0.05, nperm=100, top.ranked=100)
plotValidation(validation.result=continuous.all.100, method="continuous_top100", cell="all")

binary.all <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna.bin, gdsc.sig.rna=gdsc.sig.rna2.bin, cell="all", method="binary", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=binary.all, method="binary", cell="all")
binary.all.100 <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna.bin, gdsc.sig.rna=gdsc.sig.rna2.bin, cell="all", method="binary", drugs=drugs, fdr.cut.off=0.05, nperm=100, top.ranked=100)
plotValidation(validation.result=binary.all.100, method="binary_top100", cell="all")


test <- wilcox.test(continuous.all.100$validation, binary.all.100$validation, alternative="greater", paired=TRUE)
test$p.value


test <- wilcox.test(continuous.all.100$validation, continuous.common.100$validation, alternative="greater", paired=TRUE)
test$p.value

continuous.common <- biomarkers.validation(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=continuous.common, method="continuous", cell="common")
continuous.common.100 <- biomarkers.validation(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous", drugs=drugs, fdr.cut.off=0.05, nperm=100, top.ranked=100)
plotValidation(validation.result=continuous.common.100, method="continuous_top100", cell="common")


binary.commom <- biomarkers.validation(ccle.sig.rna=common.ccle.sig.rna.bin, gdsc.sig.rna=common.gdsc.sig.rna2.bin, cell="common", method="binary", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=binary.commom, method="binary", cell="common")
binary.commom.100 <- biomarkers.validation(ccle.sig.rna=common.ccle.sig.rna.bin, gdsc.sig.rna=common.gdsc.sig.rna2.bin, cell="common", method="binary", drugs=drugs, fdr.cut.off=0.05, nperm=100, top.ranked=100)
plotValidation(validation.result=binary.commom.100, method="binary_top100", cell="common")


test <- wilcox.test(continuous.common.100$validation, binary.commom.100$validation, alternative="greater", paired=TRUE)
test$p.value



cnv.fetures <- intersect(rownames(fData(CCLE.clarified@molecularProfiles$cnv)), rownames(fData(GDSC.clarified@molecularProfiles$cnv)))

myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_cnv.RData")
if(!file.exists(myf))
{
  
  gdsc.sig.cnv <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="cnv", drugs=drugs, features=cnv.fetures, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(gdsc.sig.cnv, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_cnv_binary.RData")
if(!file.exists(myf))
{
  
  gdsc.sig.cnv.bin <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="cnv", drugs=drugs, features=cnv.fetures, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  save(gdsc.sig.cnv.bin, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "common_gdsc_sig_cnv.RData")
if(!file.exists(myf))
{
  
  common.gdsc.sig.cnv <- drugSensitivitySig(pSet=common.clarified$GDSC, mDataType="cnv", drugs=drugs, features=cnv.fetures, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(common.gdsc.sig.cnv, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "common_gdsc_sig_cnv_binary.RData")
if(!file.exists(myf))
{
  
  common.gdsc.sig.cnv.bin <- drugSensitivitySig(pSet=common.clarified$GDSC, mDataType="cnv", drugs=drugs, features=cnv.fetures, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  save(common.gdsc.sig.cnv.bin, file=myf)
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

myf <- file.path(saveres, "PSets", "Sigs", "ccle_sig_cnv_binary.RData")
if(!file.exists(myf))
{
  ccle.sig.cnv.bin <- drugSensitivitySig(pSet=CCLE.clarified, mDataType="cnv", drugs=drugs, features=cnv.fetures, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  save(ccle.sig.cnv.bin, file=myf)
} else{
  load(myf)
}

myf <- file.path("PSets", "Sigs", "common_ccle_sig_cnv.RData")
if(!file.exists(myf))
{
  
  common.ccle.sig.cnv <- drugSensitivitySig(pSet=common.clarified$CCLE, mDataType="cnv", drugs=drugs, features=cnv.fetures, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4)
  save(common.ccle.sig.cnv, file=myf)
}else{
  load(myf)
}

myf <- file.path("PSets", "Sigs", "common_ccle_sig_cnv_binary.RData")
if(!file.exists(myf))
{
  
  common.ccle.sig.cnv.bin <- drugSensitivitySig(pSet=common.clarified$CCLE, mDataType="cnv", drugs=drugs, features=cnv.fetures, sensitivity.measure="auc_published", molecular.summary.stat="median", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  save(common.ccle.sig.cnv.bin, file=myf)
}else{
  load(myf)
}


myf <- file.path(saveres, "PSets", "Sigs", "gdsc_mutation.RData")
if(!file.exists(myf))
{
  gdsc.sig.mutation <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="mutation", sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4)
  
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

myf <- file.path(saveres, "PSets", "Sigs", "gdsc_sig_mutation_binary.RData")
if(!file.exists(myf))
{
  gdsc.sig.mutation.bin <- drugSensitivitySig(pSet=GDSC.clarified, mDataType="mutation", drugs=drugs, features=features.mutation, sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  
  save(gdsc.sig.mutation.bin, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "common_gdsc_sig_mutation.RData")
if(!file.exists(myf))
{
  common.gdsc.sig.mutation <- drugSensitivitySig(pSet=common.clarified$GDSC, mDataType="mutation", drugs=drugs, features=features.mutation, sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4)
  
  save(common.gdsc.sig.mutation, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "common_gdsc_sig_mutation_binary.RData")
if(!file.exists(myf))
{
  common.gdsc.sig.mutation.bin <- drugSensitivitySig(pSet=common.clarified$GDSC, mDataType="mutation", drugs=drugs, features=features.mutation, sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  
  save(common.gdsc.sig.mutation.bin, file=myf)
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

myf <- file.path(saveres, "PSets", "Sigs", "ccle_sig_mutation_binary.RData")
if(!file.exists(myf))
{
  ccle.sig.mutation.bin <- drugSensitivitySig(pSet=CCLE.clarified, mDataType="mutation", drugs=drugs, features=features.mutation, sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  
  save(ccle.sig.mutation.bin, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "common_ccle_sig_mutation.RData")
if(!file.exists(myf))
{
  common.ccle.sig.mutation <- drugSensitivitySig(pSet=common.clarified$CCLE, mDataType="mutation", drugs=drugs, features=features.mutation, sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4)
  
  save(common.ccle.sig.mutation, file=myf)
}else{
  load(myf)
}

myf <- file.path(saveres, "PSets", "Sigs", "common_ccle_sig_mutation_binary.RData")
if(!file.exists(myf))
{
  common.ccle.sig.mutation.bin <- drugSensitivitySig(pSet=common.clarified$CCLE, mDataType="mutation", drugs=drugs, features=features.mutation, sensitivity.measure="auc_published", molecular.summary.stat="or", sensitivity.summary.stat="median", nthread=4, sensitivity.cutoff=0.2)
  
  save(common.ccle.sig.mutation.bin, file=myf)
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

##all continuous
all.types.biomarkers <- integrateDrugBasedBiomarkers(method="continuous_all", drugs=drugs, cut.off=0.05, ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, ccle.sig.cnv=ccle.sig.cnv, gdsc.sig.cnv=gdsc.sig.cnv, ccle.sig.mutation=ccle.sig.mutation, gdsc.sig.mutation=gdsc.sig.mutation)

integrateEstimatesScatterplot(biomarkers=all.types.biomarkers, method="continuous_all", drugs=drugs, fdr.cut.off=0.05)

fdrBarplot(biomarkers=all.types.biomarkers, method="continuous_all", drugs=drugs)
integrate.all.validation <- integrateBiomarkersValidation(biomarkers=all.types.biomarkers, method="continuous_all", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=integrate.all.validation, method="integrate_continuous", cell="all")

##common continuous
all.types.biomarkers <- integrateDrugBasedBiomarkers(method="continuous_common", drugs=drugs, cut.off=0.05, ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, ccle.sig.cnv=common.ccle.sig.cnv, gdsc.sig.cnv=common.gdsc.sig.cnv, ccle.sig.mutation=common.ccle.sig.mutation, gdsc.sig.mutation=common.gdsc.sig.mutation)

integrateEstimatesScatterplot(biomarkers=all.types.biomarkers, method="continuous_common", drugs=drugs, fdr.cut.off=0.05)

fdrBarplot(biomarkers=all.types.biomarkers, method="continuous_common", drugs=drugs)
integrate.all.validation <- integrateBiomarkersValidation(biomarkers=all.types.biomarkers, method="continuous_common", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=integrate.all.validation, method="integrate_continuous", cell="common")


##all binary
all.types.biomarkers <- integrateDrugBasedBiomarkers(method="binary_all", drugs=drugs, cut.off=0.05, ccle.sig.rna=ccle.sig.rna.bin, gdsc.sig.rna=gdsc.sig.rna2.bin, ccle.sig.cnv=ccle.sig.cnv.bin, gdsc.sig.cnv=gdsc.sig.cnv.bin, ccle.sig.mutation=ccle.sig.mutation.bin, gdsc.sig.mutation=gdsc.sig.mutation.bin)

integrateEstimatesScatterplot(biomarkers=all.types.biomarkers, method="binary_all", drugs=drugs, fdr.cut.off=0.05)

fdrBarplot(biomarkers=all.types.biomarkers, method="binary_all", drugs=drugs)
integrate.all.validation <- integrateBiomarkersValidation(biomarkers=all.types.biomarkers, method="binary_all", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=integrate.all.validation, method="integrate_binary", cell="all")

##common binary
all.types.biomarkers <- integrateDrugBasedBiomarkers(method="binary_common", drugs=drugs, cut.off=0.05, ccle.sig.rna=common.ccle.sig.rna.bin, gdsc.sig.rna=common.gdsc.sig.rna2.bin, ccle.sig.cnv=common.ccle.sig.cnv.bin, gdsc.sig.cnv=common.gdsc.sig.cnv.bin, ccle.sig.mutation=common.ccle.sig.mutation.bin, gdsc.sig.mutation=common.gdsc.sig.mutation.bin)

integrateEstimatesScatterplot(biomarkers=all.types.biomarkers, method="binary_common", drugs=drugs, fdr.cut.off=0.05)

fdrBarplot(biomarkers=all.types.biomarkers, method="binary_common", drugs=drugs)
integrate.all.validation <- integrateBiomarkersValidation(biomarkers=all.types.biomarkers, method="binary_common", drugs=drugs, fdr.cut.off=0.05, nperm=100)
plotValidation(validation.result=integrate.all.validation, method="integrate_binary", cell="common")

