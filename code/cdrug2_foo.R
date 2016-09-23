#################################################
### Functions
#################################################

### Scatterplot with transparency
myScatterPlot <- 
  function(x, y, method=c("plain", "transparent", "smooth"), transparency=0.10, smooth.pch=".", pch=16, minp=50, col=blues9[7], smooth.col=c("white", blues9), ...) {
  require(grDevices) || stop("Library grDevices is not available!")
  method <- match.arg(method)
  if (length(col) != length(x)) {
    col <- rep(col, length.out=length(x))
  }
  ccix <- complete.cases(x, y)
  x <- x[ccix]
  y <- y[ccix]
  col <- col[ccix]
  
  if (sum(ccix) < minp) {
    ## too few points, no transparency, no smoothing
    if (sum(ccix) > 0) { rr <- plot(x=x, y=y, col=col, pch=pch, ...) } else { rr <- plot(x=x, y=y, col=col, pch=pch, ...) }
  } else {
    ## enough data points
    switch(method,
           "plain"={
             rr <- plot(x=x, y=y, col=col, pch=pch, ...)
           },
           "transparent"={
             myrgb <- sapply(col, grDevices::col2rgb, alpha=FALSE) / 255
             myrgb <- apply(myrgb, 2, function (x, transparency) {
               return (rgb(red=x[1], green=x[2], blue=x[3], alpha=transparency, maxColorValue=1))
             }, transparency=transparency)
             rr <- plot(x=x, y=y, pch=pch, col=myrgb, ...)
           },
           "smooth"={
             rr <- smoothScatter(x=x, y=y, col="lightgray", colramp=colorRampPalette(smooth.col), pch=smooth.pch, ...)
           }
    )
  }
  
  invisible(rr)
}


#################################################
## compute consistency of sensitivity measurments using several methods
computeConsistencySensitivity <- 
  function (x, y, type=c("auc", "ic50"), drugs, concentrations, cutoff, cutoff.cytotoxic) {
  
  type <- match.arg(type)
  mystats <- c("pcc.full", "pcc.sens", "scc.full", "scc.sens", "dxy.full", "dxy.sens", "cosine", "mcc", "kappa", "cramerv", "inform")
  consis <- array(NA, dim=c(length(drugs), length(mystats), 2), dimnames=list(names(drugs), mystats, c("estimate", "p")))
  for (i in names(drugs)) {
    x.i <- x[i, ]
    y.i <- y[i, ]
    if (type == "ic50") {
      if (!missing(concentrations)) {
        ## truncate ic50= and ic50=Inf by the min and max concentrations tested in each study
        ## ic50 = 0 when the drug dose-response curve starts below 50% viability
        ## ic50 = Inf when 50% viability is never reached
        x.i[x.i < concentrations[[1]][i, 1]] <- concentrations[[1]][i, 1]
        x.i[x.i > concentrations[[1]][i, 2]] <- concentrations[[1]][i, 2]
        y.i[y.i < concentrations[[2]][i, 1]] <- concentrations[[2]][i, 1]
        y.i[y.i > concentrations[[2]][i, 2]] <- concentrations[[2]][i, 2]
      }
      ## - log10 IC50 in microMolar
      x.i <- - log10(x.i / 1000)
      y.i <- - log10(y.i / 1000)
    }
    ## binarization
    cc <- ifelse(drugs[i] == 1, cutoff.cytotoxic, cutoff)
    x.i.bin <- factor(ifelse(x.i > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))
    y.i.bin <- factor(ifelse(y.i > cc, "sensitive", "resistant"), levels=c("resistant","sensitive"))
    ## consistency measures
    iix <- !(is.na(x.i.bin) | is.na(y.i.bin)) & !(x.i.bin == "resistant" & y.i.bin == "resistant")
    ccix <- complete.cases(x.i, y.i)
    ## pcc full
    res <- try(cor.test(x.i, y.i, method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.full", ] <- c(res$estimate, res$p.value)
    }
    ## pcc sensitive
    res <- try(cor.test(x.i[iix], y.i[iix], method="pearson", use="complete.obs", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "pcc.sens", ] <- c(res$estimate, res$p.value)
    }
    ## scc full
    res <- try(cor.test(x.i, y.i, method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.full", ] <- c(res$estimate, res$p.value)
    }
    ## scc sensitive
    res <- try(cor.test(x.i[iix], y.i[iix], method="spearman", use="complete.obs", alternative="greater", exact=FALSE), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "scc.sens", ] <- c(res$estimate, res$p.value)
    }
    ## dxy
    res <- try(mRMRe::correlate(X=x.i, Y=y.i, method="cindex", alternative="greater"), silent=TRUE) 
    if (class(res) != "try-error") {
      consis[i, "dxy.full", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## dxy sensitive
     res <- try(mRMRe::correlate(X=x.i[iix], Y=y.i[iix], method="cindex", alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "dxy.sens", ] <- c((res$estimate - 0.5) * 2, res$p)
    }
    ## cosine
    res <- try(PharmacoGx::cosinePerm(x=x.i[ccix] - cc, y=y.i[ccix] - cc, nperm=nperm, nthread=nbcore, alternative="greater"), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "cosine", ] <- c(res$estimate, res$p.value)
    }
    ## mcc
    #res <- try(PharmacoGx::mcc(x=x.i.bin, y=y.i.bin, nperm=nperm, alternative="greater", nthread=nbcore), silent=TRUE)
    res <- try(mcc(x=x.i.bin, y=y.i.bin, nperm=nperm, nthread=nbcore), silent=TRUE)
    if (class(res) != "try-error") {
      consis[i, "mcc", ] <- c(res$estimate, res$p.value)
    }
    ## kappa
    tt <- table(x=x.i.bin, y=y.i.bin)
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
    ## Powers, David M W (2011). "Evaluation: From Precision, Recall and F-Measure to ROC, Informedness, Markedness & Correlation" (PDF). Journal of Machine Learning Technologies 2 (1): 37–63.
    ## On the Statistical Consistency of Algorithms for Binary Classification under Class Imbalance 
    ## Authors: Aditya Menon, Harikrishna Narasimhan, Shivani Agarwal and Sanjay Chawla
    ## Conference: Proceedings of the 30th International Conference on Machine Learning (ICML-13), Year: 2013, Pages: 603-611
    ## TODO: implement permutation test
    res <- try(caret::sensitivity(data=x.i.bin, reference=y.i.bin) + caret::specificity(data=x.i.bin, reference=y.i.bin) - 1)
    if (class(res) != "try-error") {
      consis[i, "inform", "estimate"] <- res
    }
  }
  return (consis)
}


#################################################
## Determine the specificity of bimomarkers to the studies across drugs using Jaccard Index
biomarkersSpecificity <- 
  function(ccle.sig.rna, ccle.sig.cnv, ccle.sig.mutation, gdsc.sig.rna, gdsc.sig.cnv, gdsc.sig.mutation, cell, method, drugs, fdr.cut.off) {
    require(xtable)
    require(abind)
    dd <- intersectList(dimnames(ccle.sig.rna)[[3]], dimnames(ccle.sig.cnv)[[3]], dimnames(ccle.sig.mutation)[[3]])
    ccle.sigs <- NULL
    temp <- ccle.sig.rna
    rownames(temp) <- paste0(rownames(temp), "_rna")
    ccle.sigs <- abind(ccle.sigs, temp[,,dd], along=1)
    
    temp <- ccle.sig.cnv
    rownames(temp) <- paste0(rownames(temp), "_cnv")
    ccle.sigs <- abind(ccle.sigs, temp[,,dd], along=1)
    
    temp <- ccle.sig.mutation
    rownames(temp) <- paste0(rownames(temp), "_mutation")
    ccle.sigs <- abind(ccle.sigs, temp[,,dd], along=1)
    
    gdsc.sigs <- NULL
    temp <- gdsc.sig.rna
    rownames(temp) <- paste0(rownames(temp), "_rna")
    gdsc.sigs <- abind(gdsc.sigs, temp[,,dd], along=1)
    
    temp <- gdsc.sig.cnv
    rownames(temp) <- paste0(rownames(temp), "_cnv")
    gdsc.sigs <- abind(gdsc.sigs, temp[,,dd], along=1)
    
    temp <- gdsc.sig.mutation
    rownames(temp) <- paste0(rownames(temp), "_mutation")
    gdsc.sigs <- abind(gdsc.sigs, temp[,,dd], along=1)
    
    datset.specific.biomarkers <- matrix(NA, nrow=length(drugs), ncol=5)
    rownames(datset.specific.biomarkers) <- drugs
    colnames(datset.specific.biomarkers) <- c("# GDSC", "# CCLE", "% GDSC", "% CCLE", "% Both")

    for (drug in drugs)
    {
      
      top.ccle.sigs <- ccle.sigs[which((ccle.sigs[, drug, "fdr"] < fdr.cut.off) | (gdsc.sigs[, drug, "fdr"] < fdr.cut.off)), drug, ]
      if (is.null(dim(top.ccle.sigs))) {
        if(length(top.ccle.sigs) > 0){
          top.ccle.sigs <- as.matrix(t(top.ccle.sigs))
          rownames(top.ccle.sigs) <- names(which((ccle.sigs[, drug, "fdr"] < fdr.cut.off) | (gdsc.sigs[, drug, "fdr"] < fdr.cut.off)))
        }else{
          top.ccle.sigs <- matrix(NA, ncol=length(dimnames(ccle.sigs)[[3]]), dimnames=list(NA,dimnames(ccle.sigs)[[3]]))
        }
      }
      top.gdsc.sigs <- gdsc.sigs[which((ccle.sigs[, drug, "fdr"] < fdr.cut.off) | (gdsc.sigs[, drug, "fdr"] < fdr.cut.off)), drug, ]
      if (is.null(dim(top.gdsc.sigs))) {
        if(length(top.gdsc.sigs) > 0){
          top.gdsc.sigs <- as.matrix(t(top.gdsc.sigs))
          rownames(top.gdsc.sigs) <- names(which((ccle.sigs[, drug, "fdr"] < fdr.cut.off) | (gdsc.sigs[, drug, "fdr"] < fdr.cut.off)))
        }else{
          top.gdsc.sigs <- matrix(NA, ncol=length(dimnames(gdsc.sigs)[[3]]), dimnames=list(NA,dimnames(gdsc.sigs)[[3]]))
        }
      }
      rownames(top.ccle.sigs) <- paste(rownames(top.ccle.sigs), sign(top.ccle.sigs[,"estimate"]), sep="_")
      rownames(top.gdsc.sigs) <- paste(rownames(top.gdsc.sigs), sign(top.gdsc.sigs[,"estimate"]), sep="_")
      
      common <- intersect(rownames(top.ccle.sigs), rownames(top.gdsc.sigs))
      common <- length(which(top.ccle.sigs[common, "fdr"] < fdr.cut.off & top.gdsc.sigs[common, "fdr"] < fdr.cut.off))
      
      gdsc <- length(which(top.gdsc.sigs[,"fdr"] < fdr.cut.off))
      ccle <- length(which(top.ccle.sigs[,"fdr"] < fdr.cut.off))
      
      all <- gdsc + ccle - common
      
      datset.specific.biomarkers[drug, 1:5] <- c("# GDSC"=gdsc- common, 
                                              "# CCLE"=ccle- common, 
                                              "% GDSC"= round((gdsc- common)/all, digits=2) *100, 
                                              "% CCLE"=round((ccle- common)/all, digits=2) * 100,
                                              "% Both"=round(common/all, digits=2) * 100)
    }
    Jaccard <- datset.specific.biomarkers[, "% Both"]
    xtable::print.xtable(xtable::xtable(datset.specific.biomarkers, digits=0), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=file.path(saveres, sprintf("dataset_specific_biomarkers_%s_%s.tex", method, cell)), append=FALSE)
    return(Jaccard)
}

###############################################
###check known biomarkers


knownBiomarkersCheck <-
  function(ccle.sig.rna, gdsc.sig.rna, ccle.sig.mutation, gdsc.sig.mutation, gdsc.sig.fusion, ccle.sig.cnv, gdsc.sig.cnv, method, cell)
    { ### for each known biomarker, estimate gene-drug association for mutation, fusion and expression
    known.biomarkers <- read.csv(file.path("data", "known_biomarkers.csv"), stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, na.strings=c("", " ", "NA"))
    gdsc.known.biomarkers <- known.biomarkers
    gdsc.known.biomarkers[,"probe"] <- rownames(fData(GDSC@molecularProfiles$rna))[match(gdsc.known.biomarkers[, "gene"], fData(GDSC@molecularProfiles$rna)[, "Symbol"])]
    
    gdsc.known.biomarkers[,"mutation"] <- NA
    gdsc.known.biomarkers[,"fusion"] <- NA
    gdsc.known.biomarkers[,"rna"] <- NA
    for (i in 1:nrow(gdsc.known.biomarkers)) {
      gene <- gdsc.known.biomarkers[i , "gene"]
      probe <- gdsc.known.biomarkers[i , "probe"]
      if (gene %in% rownames(exprs(GDSC@molecularProfiles$mutation)) && !all(is.na(exprs(GDSC@molecularProfiles$mutation)[gene, ]))) {
        gdsc.known.biomarkers[i,"mutation"] <- 1
      }
      if (gene %in% rownames(exprs(GDSC@molecularProfiles$fusion)) && !all(is.na(exprs(GDSC@molecularProfiles$fusion)[gene, ]))) {
        gdsc.known.biomarkers[i,"fusion"] <- 1
      }
      if (probe %in% rownames(exprs(GDSC@molecularProfiles$rna)) && !all(is.na(exprs(GDSC@molecularProfiles$rna)[probe, ]))) {
        gdsc.known.biomarkers[i,"rna"] <- 1
      }
      if (gene %in% rownames(exprs(GDSC@molecularProfiles$cnv)) && !all(is.na(exprs(GDSC@molecularProfiles$cnv)[gene, ]))) {
        gdsc.known.biomarkers[i,"cnv"] <- 1
      }
    }
    for(i in 1:nrow(gdsc.known.biomarkers)) {
      if(!is.na(gdsc.known.biomarkers[i, "mutation"])) {
        gdsc.known.biomarkers[i, "mutation.pvalue"] <- gdsc.sig.mutation@.Data[gdsc.known.biomarkers[i,"gene"], gdsc.known.biomarkers[i,"drug"], "pvalue"]
        gdsc.known.biomarkers[i, "mutation.estimate"] <- gdsc.sig.mutation@.Data[gdsc.known.biomarkers[i,"gene"], gdsc.known.biomarkers[i,"drug"], "estimate"]
      }
      if(!is.na(gdsc.known.biomarkers[i, "fusion"])) {
        gdsc.known.biomarkers[i, "fusion.pvalue"] <- gdsc.sig.fusion@.Data[gdsc.known.biomarkers[i,"gene"], gdsc.known.biomarkers[i,"drug"], "pvalue"]
        gdsc.known.biomarkers[i, "fusion.estimate"] <- gdsc.sig.fusion@.Data[gdsc.known.biomarkers[i,"gene"], gdsc.known.biomarkers[i,"drug"], "estimate"]
      }
      if(!is.na(gdsc.known.biomarkers[i, "rna"])) {
        gdsc.known.biomarkers[i, "rna.pvalue"] <- gdsc.sig.rna@.Data[gdsc.known.biomarkers[i,"probe"], gdsc.known.biomarkers[i,"drug"], "pvalue"]
        gdsc.known.biomarkers[i, "rna.estimate"] <- gdsc.sig.rna@.Data[gdsc.known.biomarkers[i,"probe"], gdsc.known.biomarkers[i,"drug"], "estimate"]
      }
    }
    
    ##ccle
    ccle.known.biomarkers <- known.biomarkers
    ccle.known.biomarkers[,"probe"] <- rownames(fData(CCLE@molecularProfiles$rna))[match(gdsc.known.biomarkers[, "gene"], fData(CCLE@molecularProfiles$rna)[, "Symbol"])]
    
    ccle.known.biomarkers[,"mutation"] <- NA
    ccle.known.biomarkers[,"fusion"] <- NA
    ccle.known.biomarkers[,"rna"] <- NA
    ccle.known.biomarkers[,"cnv"] <- NA
    
    ### update CCLE PSet
    celline.bcrabl <- c("K-562", "KYO-1", "EM-3", "AR230", "KCL22", "BV-173", "CML-T1", "EM-2", "KU812", "LAMA-84", "MEG-01")
    fusion.ccle <- CCLE@molecularProfiles$mutation
    exprs(fusion.ccle) <- matrix("0", nrow=1, ncol=ncol(exprs(fusion.ccle)), dimnames=list("BCR_ABL", colnames(exprs(fusion.ccle))))
    exprs(fusion.ccle)["BCR_ABL", colnames(exprs(fusion.ccle)) %in% celline.bcrabl] <- "BCR Exon_13 to ABL Exon_2"
    myfdata <- data.frame("Symbol"=rownames(fusion.ccle), "gene_biotype"="protein_coding", stringsAsFactors=FALSE)
    rownames(myfdata) <- rownames(fusion.ccle)
    mypdata <- data.frame(cbind("batchid"=NA, "cellid"=colnames(fusion.ccle)), stringsAsFactors=FALSE)
    rownames(mypdata) <- colnames(fusion.ccle)
    fData(fusion.ccle) <- myfdata
    pData(fusion.ccle) <- mypdata
    CCLE@molecularProfiles$fusion <- fusion.ccle
    annotation(CCLE@molecularProfiles$fusion) <- "fusion"
    
    for (i in 1:nrow(ccle.known.biomarkers)) {
      gene <- ccle.known.biomarkers[i , "gene"]
      probe <- ccle.known.biomarkers[i , "probe"]
      if (gene %in% rownames(exprs(CCLE@molecularProfiles$mutation)) && !all(is.na(exprs(CCLE@molecularProfiles$mutation)[gene, ]))) {
        ccle.known.biomarkers[i,"mutation"] <- 1
      }
      if (gene %in% rownames(exprs(CCLE@molecularProfiles$fusion)) && !all(is.na(exprs(CCLE@molecularProfiles$fusion)[gene, ]))) {
        ccle.known.biomarkers[i,"fusion"] <- 1
      }
      if (probe %in% rownames(exprs(CCLE@molecularProfiles$rna)) && !all(is.na(exprs(CCLE@molecularProfiles$rna)[probe, ]))) {
        ccle.known.biomarkers[i,"rna"] <- 1
      }
    }

   
    ccle.sig.fusion <- drugSensitivitySig(pSet=CCLE, mDataType="fusion", drugs=unique(ccle.known.biomarkers[ ,"drug"]), features=unique(ccle.known.biomarkers[which(!is.na(ccle.known.biomarkers[,"fusion"])),"gene"]), sensitivity.measure="auc_published", molecular.summary.stat="or")
    for(i in 1:nrow(ccle.known.biomarkers)) {
      if(!is.na(ccle.known.biomarkers[i, "mutation"]) && ccle.known.biomarkers[i,"gene"] %in% rownames(ccle.sig.mutation@.Data)) {
        ccle.known.biomarkers[i, "mutation.pvalue"] <- ccle.sig.mutation@.Data[ccle.known.biomarkers[i,"gene"], ccle.known.biomarkers[i,"drug"], "pvalue"]
        ccle.known.biomarkers[i, "mutation.estimate"] <- ccle.sig.mutation@.Data[ccle.known.biomarkers[i,"gene"], ccle.known.biomarkers[i,"drug"], "estimate"]
      }
      if(!is.na(ccle.known.biomarkers[i, "fusion"])) {
        ccle.known.biomarkers[i, "fusion.pvalue"] <- ccle.sig.fusion@.Data[ccle.known.biomarkers[i,"gene"], ccle.known.biomarkers[i,"drug"], "pvalue"]
        ccle.known.biomarkers[i, "fusion.estimate"] <- ccle.sig.fusion@.Data[ccle.known.biomarkers[i,"gene"], ccle.known.biomarkers[i,"drug"], "estimate"]
      }
      if(!is.na(ccle.known.biomarkers[i, "rna"])) {
        ccle.known.biomarkers[i, "rna.pvalue"] <- ccle.sig.rna@.Data[ccle.known.biomarkers[i,"probe"], ccle.known.biomarkers[i,"drug"], "pvalue"]
        ccle.known.biomarkers[i, "rna.estimate"] <- ccle.sig.rna@.Data[ccle.known.biomarkers[i,"probe"], ccle.known.biomarkers[i,"drug"], "estimate"]
      }
    }
    
    cutoff <- 0.05
    xx <- NULL
    for(i in 1:nrow(ccle.known.biomarkers)) {
      ccle.min <- names(which.min(ccle.known.biomarkers[i, grep("pvalue", colnames(ccle.known.biomarkers))]))
      gdsc.min <- names(which.min(gdsc.known.biomarkers[i, grep("pvalue", colnames(gdsc.known.biomarkers))]))
      tt <- NA
      if(ccle.min == gdsc.min) {
        if(!is.na(ccle.known.biomarkers[i, ccle.min]) & !is.na(gdsc.known.biomarkers[i, ccle.min])){
          if(ccle.known.biomarkers[i, ccle.min] < cutoff & gdsc.known.biomarkers[i, ccle.min] < cutoff){
            tt <- "YES"
          }else if(ccle.known.biomarkers[i, ccle.min] >= cutoff & gdsc.known.biomarkers[i, ccle.min] >= cutoff){
            tt <- "NS"
          }else{
            tt <- "NO"
          }
        }
        rr <- c(known.biomarkers[i, "drug"],
                known.biomarkers[i, "gene"],
                gsub(".pvalue", "", ccle.min),
                gdsc.known.biomarkers[i, gsub(".pvalue", ".estimate", ccle.min)],
                gdsc.known.biomarkers[i, ccle.min],
                ccle.known.biomarkers[i, gsub(".pvalue", ".estimate", ccle.min)],
                ccle.known.biomarkers[i, ccle.min], 
                tt)
        xx <- rbind(xx, rr)
      } else{
        tt <- NA
        if(!is.na(ccle.known.biomarkers[i, gdsc.min]) & !is.na(gdsc.known.biomarkers[i, gdsc.min])){
          if(ccle.known.biomarkers[i, gdsc.min] < cutoff & gdsc.known.biomarkers[i, gdsc.min] < cutoff){
            tt <- "YES"
          }else if(ccle.known.biomarkers[i, gdsc.min] >= cutoff & gdsc.known.biomarkers[i, gdsc.min] >= cutoff){
            tt <- "NS"
          }else{
            tt <- "NO"
          }
        }
        rr <- c(known.biomarkers[i, "drug"],
                known.biomarkers[i, "gene"],
                gsub(".pvalue", "", gdsc.min),
                gdsc.known.biomarkers[i, gsub(".pvalue", ".estimate", gdsc.min)],
                gdsc.known.biomarkers[i, gdsc.min],
                ccle.known.biomarkers[i, gsub(".pvalue", ".estimate", gdsc.min)],
                ccle.known.biomarkers[i, gdsc.min], 
                tt)
        xx <- rbind(xx, rr)
        tt <- NA
        if(!is.na(ccle.known.biomarkers[i, ccle.min]) & !is.na(gdsc.known.biomarkers[i, ccle.min])){
          if(ccle.known.biomarkers[i, ccle.min] < cutoff & gdsc.known.biomarkers[i, ccle.min] < cutoff){
            tt <- "YES"
          }else if(ccle.known.biomarkers[i, ccle.min] >= cutoff & gdsc.known.biomarkers[i, ccle.min] >= cutoff){
            tt <- "NS"
          }else{
            tt <- "NO"
          }
        }
        rr <- c(known.biomarkers[i, "drug"],
                known.biomarkers[i, "gene"],
                gsub(".pvalue", "", ccle.min),
                gdsc.known.biomarkers[i, gsub(".pvalue", ".estimate", ccle.min)],
                gdsc.known.biomarkers[i, ccle.min],
                ccle.known.biomarkers[i, gsub(".pvalue", ".estimate", ccle.min)],
                ccle.known.biomarkers[i, ccle.min], 
                tt)
        xx <- rbind(xx, rr)
      }
    }
    colnames(xx) <- c("Drug", "Gene", "Type", "gdsc effect size", "gdsc pvalue", "CCLE effect size", "CCLE pvalue", "Reproducibility")
    rownames(xx) <- 1:nrow(xx)
    nilotinib <- which(xx[,"Drug"] == "Nilotinib")
    rr <- xx[nilotinib, ]
    xx <- xx[-nilotinib, ]
    xx <- rbind(rr, xx)
    xx <- as.data.frame(xx, stringsAsFactors=FALSE)
    
    xx[,"gdsc effect size"] <- as.numeric(xx[,"gdsc effect size"])
    xx[,"CCLE effect size"] <- as.numeric(xx[,"CCLE effect size"])
    xx[,"gdsc pvalue"] <- as.numeric(xx[,"gdsc pvalue"])
    xx[,"CCLE pvalue"] <- as.numeric(xx[,"CCLE pvalue"])
    xtable::print.xtable(xtable::xtable(xx, digits=c(0, 0, 0, 0, 2, -1, 2, -1, 0)), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(saveres, "known_biomarkers.tex"), append=FALSE)
    
  }

#################################################
## Create an excel file for the statistics of all the gene-drug associations in GDSC and CCLE
integrateDrugBasedBiomarkers <-
  function(method, drugs, cut.off, ccle.sig.rna, gdsc.sig.rna, ccle.sig.cnv, gdsc.sig.cnv, ccle.sig.mutation, gdsc.sig.mutation) {
    require(WriteXLS)
    all.biomarkers <- list()
    for(drug in drugs) {
      ccle.biomarkers <- ccle.sig.rna[features, drug, ]
      colnames(ccle.biomarkers) <- paste0("CCLE_", colnames(ccle.biomarkers))
      
      gdsc.biomarkers <- gdsc.sig.rna[features, drug, ]
      colnames(gdsc.biomarkers) <- paste0("GDSC_", colnames(gdsc.biomarkers))
      
      biomarkers <- cbind("Symbol"=NA, "Type"="Expression", gdsc.biomarkers, ccle.biomarkers, "Specificity"="Non significant")
      biomarkers[,"Symbol"] <- featureInfo(CCLE, "rna")[rownames(biomarkers), "Symbol"]
      biomarkers.rna <- biomarkers
      rownames(biomarkers.rna) <- sprintf("%s_rna", rownames(biomarkers))
      
      
      ccle.biomarkers <- ccle.sig.cnv[cnv.fetures, drug, ]
      colnames(ccle.biomarkers) <- paste0("CCLE_", colnames(ccle.biomarkers))
      
      gdsc.biomarkers <- gdsc.sig.cnv[cnv.fetures, drug, ]
      colnames(gdsc.biomarkers) <- paste0("GDSC_", colnames(gdsc.biomarkers))
      
      biomarkers <- cbind("Symbol"=NA, "Type"="CNV", gdsc.biomarkers, ccle.biomarkers, "Specificity"="Non significant")
      #biomarkers[,"Symbol"] <- featureInfo(CCLE, "rna")[rownames(biomarkers), "Symbol"]
      biomarkers[,"Symbol"] <- rownames(biomarkers)
      biomarkers.cnv <- biomarkers
      rownames(biomarkers.cnv) <- sprintf("%s_cnv", rownames(biomarkers))
      
      ccle.biomarkers <- ccle.sig.mutation[features.mutation, drug, ]
      colnames(ccle.biomarkers) <- paste0("CCLE_", colnames(ccle.biomarkers))
      
      gdsc.biomarkers <- gdsc.sig.mutation[features.mutation, drug, ]
      colnames(gdsc.biomarkers) <- paste0("GDSC_", colnames(gdsc.biomarkers))
      
      biomarkers <- cbind("Symbol"=NA, "Type"="Mutation", gdsc.biomarkers, ccle.biomarkers, "Specificity"="Non significant")
      #biomarkers[,"Symbol"] <- featureInfo(CCLE, "rna")[rownames(biomarkers), "Symbol"]
      biomarkers[,"Symbol"] <- rownames(biomarkers)
      biomarkers.mutation <- biomarkers
      rownames(biomarkers.mutation) <- sprintf("%s_mut", rownames(biomarkers))
      
      cc <- intersectList(colnames(biomarkers.rna), colnames(biomarkers.cnv), colnames(biomarkers.mutation))
      biomarkers <- rbind(biomarkers.rna[, cc], biomarkers.mutation[, cc], biomarkers.cnv[, cc])
      biomarkers[,"CCLE_fdr"] <- p.adjust(biomarkers[,"CCLE_pvalue"], method="fdr")
      biomarkers[,"GDSC_fdr"] <- p.adjust(biomarkers[,"GDSC_pvalue"], method="fdr")
      
      biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) < cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) < cut.off), "Specificity"] = "Both"
      biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) < cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) >= cut.off), "Specificity"] = "CCLE"
      biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) >= cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) < cut.off), "Specificity"] = "GDSC"
      biomarkers <- biomarkers[order(biomarkers[,"Specificity"]), ]
      
      all.biomarkers[[drug]] <- as.data.frame(biomarkers, stringsAsFactors=FALSE)
    }
    
    
    WriteXLS::WriteXLS(file.path(saveres, "all.biomarkers"), ExcelFileName=file.path(saveres, sprintf("all_biomarkers_%s.xlsx", method)), row.names=TRUE)
    return(all.biomarkers)
}


#################################################
##Plot biomarkers effect size for all molecular types across all drugs 
integrateEstimatesScatterplot <-
  function(biomarkers, method, drugs, fdr.cut.off) {
    pdf(file.path(saveres, sprintf("effect_size_%s.pdf", method)), height=16, width=16)
    par(mfrow=c(4, 4))
    mycol <- RColorBrewer::brewer.pal(n=8, name="Set3")
    
    
    for (drug in drugs)
    {
      biomarkers[[drug]][,3:12] <- apply(biomarkers[[drug]][,3:12], MARGIN=2, as.numeric)
      top.sig <- biomarkers[[drug]][which((biomarkers[[drug]][, "CCLE_fdr"] < fdr.cut.off) | (biomarkers[[drug]][, "GDSC_fdr"] < fdr.cut.off)), , drop=FALSE]
      features.drug <- rownames(top.sig)
      
      point.col <- vector(length=length(features.drug), mode="character")
      names(point.col) <- features.drug
      point.col[features.drug] <- "black"
      point.col[which((top.sig[, "CCLE_fdr"] < fdr.cut.off) & (top.sig[, "GDSC_fdr"] >= fdr.cut.off))] <- mycol[4]
      point.col[which((top.sig[, "CCLE_fdr"] >= fdr.cut.off) & (top.sig[, "GDSC_fdr"] < fdr.cut.off))] <- mycol[5]
      
      point.pch <- vector(length=length(features.drug), mode="numeric")
      names(point.pch) <- features.drug
      point.pch[features.drug] <- 1
      point.pch[which(top.sig[,"Type"] == "CNV")] <- 17
      point.pch[which(top.sig[,"Type"] == "Mutation")] <- 15
      
      rr <- max(abs(top.sig[,"CCLE_estimate"]), abs(top.sig[,"GDSC_estimate"]), na.rm=T)
      if(is.na(rr)| rr == Inf | rr == -Inf) {rr <- 1}
      par(mar=c(5,5,2,2))
      plot(NA, ylim=c(-rr, rr), xlim=c(-rr, rr), xlab="CCLE", ylab="GDSC", main=drug, cex.axis=1.5, cex.lab=1.5, cex.main=2)
      points(x=top.sig[which(point.col != "black"),"CCLE_estimate"], 
             y=top.sig[which(point.col != "black"),"GDSC_estimate"], 
             col=point.col[which(point.col != "black")], 
             pch=point.pch[which(point.col != "black")])
      points(x=top.sig[which(point.col == "black"),"CCLE_estimate"], 
             y=top.sig[which(point.col == "black"),"GDSC_estimate"], 
             col=point.col[which(point.col == "black")], 
             pch=point.pch[which(point.col == "black")])

      
      abline(0, 1, col="gray", lty=1, lwd=0.5 )
      abline(h=0, col="gray", lty=1, lwd=0.5)
      abline(v=0, col="gray", lty=1, lwd=0.5)

    legend("bottomright", legend=c(sprintf("Expression: %s", length(which(top.sig[,"Type"]== "Expression" & point.col == "black"))),
                                   sprintf("CNV: %s", length(which(top.sig[,"Type"]== "CNV" & point.col == "black"))), 
                                   sprintf("Mutation: %s", length(which(top.sig[,"Type"]== "Mutation" & point.col == "black")))), col="black", pch=c(1, 17, 15), bty="n" , cex=1.2)


    }
    plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
    legend("center", legend=c("Both significant", "GDSC significant","CCLE significant"), col=c("black", mycol[c(5, 4)]), pch=c(15, 15, 15), bty="n" , cex=2)
    dev.off()
  }


#################################################
## Validate bimarkers for both study and check validation ratio significance for each drug by permutation test
integrateBiomarkersValidation <-
  function(biomarkers, method=c("continuous","binary"), drugs, fdr.cut.off=0.05, nperm=100, top.ranked=0) {
    
    gdsc.biomarkers <- ccle.biomarkers <- matrix(NA, ncol=3, nrow=length(drugs))
    colnames(gdsc.biomarkers) <- colnames(ccle.biomarkers) <- c("features", "significant", "common")
    rownames(gdsc.biomarkers) <- rownames(ccle.biomarkers) <- drugs
    
    ccle.validation <- gdsc.validation <- matrix(NA, ncol=4, nrow=length(drugs))
    colnames(ccle.validation) <- colnames(gdsc.validation) <- c("features", "significant", "reported", "validated")
    rownames(ccle.validation) <- rownames(gdsc.validation) <- drugs
    
    
    for (drug in drugs)
    {
      ###ccle
      biomarkers[[drug]][,3:12] <- apply(biomarkers[[drug]][,3:12], MARGIN=2, as.numeric)
      top.ccle.sig.rna <- biomarkers[[drug]][which(biomarkers[[drug]][, "CCLE_fdr"] < fdr.cut.off), ,drop=FALSE]
      if(top.ranked != 0) {
        if(nrow(top.ccle.sig.rna) > 0) {
          top.ccle.sig.rna <- top.ccle.sig.rna[1:min(top.ranked, nrow(top.ccle.sig.rna)), , drop=FALSE]
        }
      }
      
#      top.ccle.sig.rna <- top.ccle.sig.rna[order(top.ccle.sig.rna[ , "CCLE_fdr"]), , drop=FALSE]
      reported <- rownames(top.ccle.sig.rna)
      validate <- top.ccle.sig.rna[which(top.ccle.sig.rna[,"GDSC_pvalue"] < 0.05), , drop=FALSE]
      ccle.validation[drug, ] <- c(nrow(biomarkers[[drug]]),
                                   nrow(top.ccle.sig.rna),
                                   length(reported),
                                   length(which(sign(validate[,"CCLE_estimate"]) == sign(validate[, "GDSC_estimate"]))))
      
      ccle.biomarkers[drug, c("features", "significant")] <- c(nrow(biomarkers[[drug]]), 
                                                               nrow(top.ccle.sig.rna))
      
      
      ###gdsc
      top.gdsc.sig.rna <- biomarkers[[drug]][which(biomarkers[[drug]][, "GDSC_fdr"] < fdr.cut.off), ,drop=FALSE]
      if(top.ranked != 0) {
        if(nrow(top.gdsc.sig.rna) > 0) {
          top.gdsc.sig.rna <- top.gdsc.sig.rna[1:min(top.ranked, nrow(top.gdsc.sig.rna)), , drop=FALSE]
        }
      }
      
      #      top.gdsc.sig.rna <- top.gdsc.sig.rna[order(top.gdsc.sig.rna[ , "GDSC_fdr"]), , drop=FALSE]
      reported <- rownames(top.gdsc.sig.rna)
      validate <- top.gdsc.sig.rna[which(top.gdsc.sig.rna[,"CCLE_pvalue"] < 0.05), , drop=FALSE]
      gdsc.validation[drug, ] <- c(nrow(biomarkers[[drug]]),
                                   nrow(top.gdsc.sig.rna),
                                   length(reported),
                                   length(which(sign(validate[,"GDSC_estimate"]) == sign(validate[, "CCLE_estimate"]))))
      
      gdsc.biomarkers[drug, c("features", "significant")] <- c(nrow(biomarkers[[drug]]), 
                                                               nrow(top.gdsc.sig.rna))
      
      features <- intersect(rownames(top.gdsc.sig.rna), rownames(top.ccle.sig.rna))
      ccle.biomarkers[drug, "common"] <- gdsc.biomarkers[drug, "common"] <- length(which(sign(top.gdsc.sig.rna[features, "GDSC_estimate"]) == sign(top.ccle.sig.rna[features, "CCLE_estimate"])))
    }
    
    gdsc.validation <- cbind(gdsc.validation,"p-value"=NA)
    ccle.validation <- cbind(ccle.validation,"p-value"=NA)
    
    ###permutation test
    for (drug in drugs)
    {
      ###ccle
      features <- rownames(biomarkers[[drug]])
      tt <- NULL
      for(i in 1:nperm) {
        biomarkers.temp <- features[sample(1:length(features), size=ccle.validation[drug, "reported"], replace=FALSE)]
        top.ccle.sig.rna <- biomarkers[[drug]][biomarkers.temp, , drop=FALSE]
        validate <- top.ccle.sig.rna[which(top.ccle.sig.rna[,"GDSC_pvalue"] < 0.05), , drop=FALSE]
        x <- ifelse(nrow(validate) > 0, length(which(sign(validate[,"GDSC_estimate"]) == sign(validate[, "CCLE_estimate"]))), 0)
        tt <- cbind(tt, x)
      }
      pvalue <- sum(ccle.validation[drug,"validated"] < tt)/ sum(!is.na(tt))
      ccle.validation[drug,"p-value"] <- ifelse(pvalue != 0, pvalue, 1/nperm)
      ccle.validation[drug,"p-value"] <- ifelse(ccle.validation[drug, "validated"] == 0, 1, ccle.validation[drug,"p-value"])
      
      ###gdsc
      tt <- NULL
      for(i in 1:nperm) {
        biomarkers.temp <- features[sample(1:length(features), size=gdsc.validation[drug, "reported"], replace=FALSE)]
        top.gdsc.sig.rna <- biomarkers[[drug]][biomarkers.temp, , drop=FALSE]
        validate <- top.gdsc.sig.rna[which(top.gdsc.sig.rna[,"CCLE_pvalue"] < 0.05), , drop=FALSE]
        x <- ifelse(nrow(validate) > 0, length(which(sign(validate[,"CCLE_estimate"]) == sign(validate[, "GDSC_estimate"]))), 0)
        tt <- cbind(tt, x)
      }
      pvalue <- sum(gdsc.validation[drug,"validated"] < tt)/ sum(!is.na(tt))
      gdsc.validation[drug,"p-value"] <- ifelse(pvalue != 0, pvalue, 1/nperm)
      gdsc.validation[drug,"p-value"] <- ifelse(gdsc.validation[drug, "validated"] == 0, 1, gdsc.validation[drug,"p-value"])
      
      
      
    }
    tt <- as.vector(rbind((gdsc.validation[,"validated"]/gdsc.validation[,"reported"]) * 100,
                          (ccle.validation[,"validated"]/ccle.validation[,"reported"]) * 100))
    
    names(tt) <- as.vector(rbind(rownames(gdsc.validation), sprintf("(GDSC=%s,CCLE=%s)", gdsc.validation[,"reported"], ccle.validation[,"reported"])))
    star <- as.vector(rbind(gdsc.validation[,"p-value"], ccle.validation[,"p-value"]))
    star[which(star <= 0.05)] <- "*"
    star[which(star > 0.05)] <- ""
    return(list("validation"=tt, "star"=star))
}

#################################################
## plot change of validated biomarkers ratio for various false discovery rates
fdrBarplot <-
  function(biomarkers, method, drugs, top.ranked=0) {
    fdrs <- c(0.5, 0.2, 0.1, 0.05, 0.01, 0.001)
    validation.result <- list()
    for (fdr in fdrs) {
      validation.result[[as.character(fdr)]] <- integrateBiomarkersValidation(biomarkers,
                                                                      method=method, 
                                                                      drugs=drugs, 
                                                                      fdr.cut.off=fdr, 
                                                                      nperm=100,
                                                                      top.ranked=top.ranked)
    }
    
    pdf(file.path(saveres, sprintf("validation_fdr_%s.pdf", method)), height=12.5, width=12.5)
    par(mfrow=c(4, 4))
    i = 0
    for(d in 1:length(drugs))
    {
      i = i + 1
      tt <- as.vector(sapply(validation.result, function(x){x$validation[c(i, i + 1)]}))
      names(tt) = ""
      
      star <- as.vector(sapply(validation.result, function(x){x$star[c(i, i + 1)]}))
      mycol <- RColorBrewer::brewer.pal(n=8, name="Set3")
      par(mar=c(5,5,2,2))
      bp <- barplot(tt, las=2, ylim=c(0,110), main=drugs[d], col=mycol[c(5, 4)], space=rep(c(0.7,0.2), length(tt)/2), xlab="FDR", ylab= "Validation rate", cex.lab=1.5, cex.main=2)
      mtext(sprintf("%s%%", fdrs * 100), side= 1, line=1, at=bp[1:length(bp) %% 2 == 0]-0.5, cex=0.8)
      text(x=bp, y= tt + 4, star, cex = 2)
      i = i + 1
    }
    plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
    legend("center", legend=c("GDSC","CCLE", "FDR: False Discovery Rate", "* The validation ratio is significant"), col=c(mycol[c(5, 4)], "white", "white"), pch=15, bty="n", cex=1.3)
    
    dev.off()
}


#################################################
## plot percentage of validated biomarkers across common drugs between GDSC and CCLE
plotValidation <- 
  function(validation.result, method, cell, main="")  {
    mycol <- RColorBrewer::brewer.pal(n=8, name="Set3")
    pdf(file.path(saveres, sprintf("validation_%s_%s.pdf", method, cell)), height=10, width=15)
    par(mar=c(18,7,5,5))
    nn <- names(validation.result$validation)
    names(validation.result$validation) <- ""
    bp <- barplot(validation.result$validation,
                  ylim=c(0,120), main=main, col=mycol[c(5, 4)], space=rep(c(0.7,0.2), length(validation.result$validation)/2), ylab="Validation rate", cex.lab=2.1, axes=FALSE)#, cex.names=2, cex.axis=2)
    axis(side=2, labels=c("0", "20%","40%","60%","80%","100%"), at=c(0, 20, 40, 60, 80 , 100), cex.axis=1.2, las=2)
    text(bp[1:length(bp) %% 2 == 1], par("usr")[3], labels = nn[1:length(nn) %% 2 == 1], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=2.1)
    text(bp[1:length(bp) %% 2 == 0], par("usr")[3], labels = nn[1:length(nn) %% 2 == 0], srt = 47, adj = c(1.1,1.1), xpd = TRUE, cex=1.7)
    
    
    legend("topright", legend=c("GDSC","CCLE"), col=mycol[c(5, 4)], pch=15, bty="n", cex=1.7)
    text(x=bp, y= validation.result$validation + 2, validation.result$star, cex = 2)
    dev.off()
    
}

