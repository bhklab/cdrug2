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
##Validate biomarkers across all common drugs between GDSC and CCLE for a specific FDR and a specific cut off on the number of top ranked biomarkers
biomarkers.validation <-
  function(ccle.sig.rna, gdsc.sig.rna, cell= c("all","common"), method=c("continuous","binary"), drugs, fdr.cut.off=0.05, nperm=100, top.ranked=0) {
    
    gdsc.biomarkers <- ccle.biomarkers <- matrix(NA, ncol=3, nrow=length(drugs))
    colnames(gdsc.biomarkers) <- colnames(ccle.biomarkers) <- c("features", "significant", "common")
    rownames(gdsc.biomarkers) <- rownames(ccle.biomarkers) <- drugs
    
    ccle.validation <- gdsc.validation <- matrix(NA, ncol=4, nrow=length(drugs))
    colnames(ccle.validation) <- colnames(gdsc.validation) <- c("features", "significant", "reported", "validated")
    rownames(ccle.validation) <- rownames(gdsc.validation) <- drugs
    
    
    for (drug in drugs)
    {
      ###ccle
      top.ccle.sig.rna <- ccle.sig.rna[which(ccle.sig.rna[ , drug, "fdr"] < fdr.cut.off), drug, ]
      if(is.null(dim(top.ccle.sig.rna))) {
        if(length(top.ccle.sig.rna) > 0){
          top.ccle.sig.rna <- as.matrix(t(top.ccle.sig.rna))
          rownames(top.ccle.sig.rna) <- dimnames(ccle.sig.rna)[[1]][which(ccle.sig.rna[ , drug, "fdr"] < fdr.cut.off)]
        }else{
          top.ccle.sig.rna <- matrix(NA, ncol=length(dimnames(ccle.sig.rna)[[3]]), dimnames=list(NA,dimnames(ccle.sig.rna)[[3]]))
        }
      } else{
        top.ccle.sig.rna <- top.ccle.sig.rna[order(top.ccle.sig.rna[ , "fdr"]), , drop=FALSE]
        if(top.ranked != 0) {
          if(nrow(top.ccle.sig.rna) > 0) {
            top.ccle.sig.rna <- top.ccle.sig.rna[1:min(top.ranked, nrow(top.ccle.sig.rna)), , drop=FALSE]
          }
        }
      }
      
      
      reported <- rownames(top.ccle.sig.rna)
      validate <- gdsc.sig.rna[reported, drug, ]
      if(is.null(dim(validate))) {
        if(length(validate) > 0){
          validate <- as.matrix(t(validate))
          rownames(validate) <- reported
        }else{
          validate <- matrix(NA, ncol=ncol(top.ccle.sig.rna), dimnames=list(NA,colnames(top.ccle.sig.rna)))
        }
      }
      validate <- validate[which(validate[,"pvalue"] < 0.05), , drop=FALSE]
      ccle.validation[drug, ] <- c(nrow(ccle.sig.rna[, drug, ]),
                                   nrow(top.ccle.sig.rna),
                                   length(reported),
                                   length(which(sign(validate[,"estimate"]) == sign(top.ccle.sig.rna[rownames(validate), "estimate"]))))
      
      ccle.biomarkers[drug, c("features", "significant")] <- c(nrow(ccle.sig.rna[, drug, ]), 
                                                               nrow(top.ccle.sig.rna))
      
      
      ###gdsc
      top.gdsc.sig.rna <- gdsc.sig.rna[which(gdsc.sig.rna[ , drug , "fdr"] < fdr.cut.off), drug, ]
      if(is.null(dim(top.gdsc.sig.rna))) {
        if(length(top.gdsc.sig.rna) > 0){
          top.gdsc.sig.rna <- as.matrix(t(top.gdsc.sig.rna))
          rownames(top.gdsc.sig.rna) <- dimnames(gdsc.sig.rna)[[1]][which(gdsc.sig.rna[ , drug, "fdr"] < fdr.cut.off)]
        }else{
          top.gdsc.sig.rna <- matrix(NA, ncol=length(dimnames(gdsc.sig.rna)[[3]]), dimnames=list(NA,dimnames(gdsc.sig.rna)[[3]]))
        }
      }else {
        top.gdsc.sig.rna <- top.gdsc.sig.rna[order(top.gdsc.sig.rna[ ,"fdr"]), , drop=FALSE]
        if(top.ranked != 0) {
          if(nrow(top.gdsc.sig.rna) > 0) {
            top.gdsc.sig.rna <- top.gdsc.sig.rna[1:min(top.ranked, nrow(top.gdsc.sig.rna)), , drop=FALSE]
          }
        }
      }
      
      
      reported <- rownames(top.gdsc.sig.rna)
      validate <- ccle.sig.rna[reported, drug, ]
      if(is.null(dim(validate))) {
        if(length(validate) > 0){
          validate <- as.matrix(t(validate))
          rownames(validate) <- reported
        }else{
          validate <- matrix(NA, ncol=ncol(top.gdsc.sig.rna), dimnames=list(NA,colnames(top.gdsc.sig.rna)))
        }
      }
      validate <- validate[which(validate[,"pvalue"] < 0.05), , drop=FALSE]
      gdsc.validation[drug, ] <-c(nrow(gdsc.sig.rna[, drug, ]),
                                  nrow(top.gdsc.sig.rna),
                                  length(reported),
                                  length(which(sign(validate[,"estimate"]) == sign(top.gdsc.sig.rna[rownames(validate), "estimate"]))))
      
      gdsc.biomarkers[drug, c("features", "significant")] <- c(nrow(gdsc.sig.rna[, drug, ]), 
                                                               nrow(top.gdsc.sig.rna))
      
      features <- intersect(rownames(top.gdsc.sig.rna), rownames(top.ccle.sig.rna))
      ccle.biomarkers[drug, "common"] <- gdsc.biomarkers[drug, "common"] <- length(which(sign(top.gdsc.sig.rna[features, "estimate"]) == sign(top.ccle.sig.rna[features, "estimate"])))
      
      ###convert gene-id to symbol for know biomarkers checking
      rownames(top.ccle.sig.rna) <- featureInfo(CCLE,"rna")[rownames(top.ccle.sig.rna), "Symbol"]
      rownames(top.gdsc.sig.rna) <- featureInfo(GDSC,"rna2")[rownames(top.gdsc.sig.rna), "Symbol"]
      
    }
    
    gdsc.validation <- cbind(gdsc.validation,"p-value"=NA)
    ccle.validation <- cbind(ccle.validation,"p-value"=NA)
    
    ###permutation test
    for (drug in drugs)
    {
      ###ccle
      features <- intersect (rownames(ccle.sig.rna[ , drug, ]), rownames(gdsc.sig.rna[ , drug, ]))
      tt <- NULL
      for(i in 1:nperm) {
        biomarkers <- features[sample(1:length(features), size=ccle.validation[drug, "reported"], replace=FALSE)]
        top.ccle.sig.rna <- ccle.sig.rna[biomarkers, drug, ]
        if(is.null(dim(top.ccle.sig.rna))) {
          if(length(top.ccle.sig.rna) > 0){
            top.ccle.sig.rna <- as.matrix(t(top.ccle.sig.rna))
            rownames(top.ccle.sig.rna) <- biomarkers
          }else{
            top.ccle.sig.rna <- matrix(NA, ncol=length(dimnames(ccle.sig.rna)[[3]]), dimnames=list(NA,dimnames(ccle.sig.rna)[[3]]))
          }
        }
        validate <- gdsc.sig.rna[biomarkers, drug, ]
        if(is.null(dim(validate))) {
          if(length(validate) > 0){
            validate <- as.matrix(t(validate))
            rownames(validate) <- biomarkers
          }else{
            validate <- matrix(NA, ncol=length(dimnames(ccle.sig.rna)[[3]]), dimnames=list(NA,dimnames(ccle.sig.rna)[[3]]))
          }
        }
        validate <- validate[which(validate[,"pvalue"] < 0.05), , drop=FALSE]
        x <- ifelse(nrow(validate) > 0, length(which(sign(validate[,"estimate"]) == sign(top.ccle.sig.rna[rownames(validate), "estimate"]))), 0)
        tt <- cbind(tt, x)
      }
      pvalue <- sum(ccle.validation[drug,"validated"] < tt)/ sum(!is.na(tt))
      ccle.validation[drug,"p-value"] <- ifelse(pvalue != 0, pvalue, 1/nperm)
      ccle.validation[drug,"p-value"] <- ifelse(ccle.validation[drug, "validated"] == 0, 1, ccle.validation[drug,"p-value"])
      
      ###gdsc
      tt <- NULL
      for(i in 1:nperm) {
        biomarkers <- features[sample(1:length(features), size=gdsc.validation[drug, "reported"], replace=FALSE)]
        top.gdsc.sig.rna <- gdsc.sig.rna[biomarkers, drug, ]
        if(is.null(dim(top.gdsc.sig.rna))) {
          if(length(top.gdsc.sig.rna) > 0){
            top.gdsc.sig.rna <- as.matrix(t(top.gdsc.sig.rna))
            rownames(top.gdsc.sig.rna) <- biomarkers
          }else{
            top.gdsc.sig.rna <- matrix(NA, ncol=length(dimnames(gdsc.sig.rna)[[3]]), dimnames=list(NA,dimnames(gdsc.sig.rna)[[3]]))
          }
        }
        validate <- ccle.sig.rna[biomarkers, drug, ]
        if(is.null(dim(validate))) {
          if(length(validate) > 0){
            validate <- as.matrix(t(validate))
            rownames(validate) <- biomarkers
          }else{
            validate <- matrix(NA, ncol=ncol(top.gdsc.sig.rna), dimnames=list(NA,colnames(top.gdsc.sig.rna)))
          }
        }
        validate <- validate[which(validate[,"pvalue"] < 0.05), , drop=FALSE]
        x <- ifelse(nrow(validate) > 0, length(which(sign(validate[,"estimate"]) == sign(top.gdsc.sig.rna[rownames(validate), "estimate"]))), 0)
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

#################################################
## plot change of validated biomarkers ratio for various cut offs on the number of top ranked biomarkers for a specific false discovery rate
topBarplot <-
  function(ccle.sig.rna, gdsc.sig.rna, cell, method, drugs, fdr.cut.off) {
    top.rankds <- c(1000, 500, 100, 50, 10, 1)
    validation.result <- list()
    for (top.ranked in top.rankds) {
      validation.result[[as.character(top.ranked)]] <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna, 
                                                                      gdsc.sig.rna=gdsc.sig.rna, 
                                                                      cell=cell, 
                                                                      method=method, 
                                                                      drugs=drugs, 
                                                                      fdr.cut.off=fdr.cut.off, 
                                                                      nperm=100,
                                                                      top.ranked=top.ranked)
    }
    
    pdf(file.path(saveres, sprintf("validation_fdr_%s_%s.pdf", method, cell)), height=12.5, width=12.5)
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
      bp <- barplot(tt, las=2, ylim=c(0,110), main=drugs[d], col=mycol[c(5, 4)], space=rep(c(0.7,0.2), length(tt)/2), xlab="TOP", ylab= "Validation rate", cex.lab=1.5, cex.main=2)
      mtext(sprintf("%s", top.rankds), side= 1, line=1, at=bp[1:length(bp) %% 2 == 0]-0.5, cex=0.8)
      text(x=bp, y= tt + 4, star, cex = 2)
      i = i + 1
    }
    plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
    legend("center", legend=c("GDSC as discovery set","CCLE as discovery set", "TOP: Top ranked biomarkers", "'*' The validation rate is significant"), col=c(mycol[c(5, 4)], "white", "white"), pch=15, bty="n", cex=1.3)
    
    dev.off()
}


#################################################
## plot change of validated biomarkers ratio for various false discovery rates
fdrBarplot <-
  function(ccle.sig.rna, gdsc.sig.rna, cell, method, drugs, top.ranked=0) {
    fdrs <- c(0.5, 0.2, 0.1, 0.05, 0.01, 0.001)
    validation.result <- list()
    for (fdr in fdrs) {
      validation.result[[as.character(fdr)]] <- biomarkers.validation(ccle.sig.rna=ccle.sig.rna, 
                                                                      gdsc.sig.rna=gdsc.sig.rna, 
                                                                      cell=cell, 
                                                                      method=method, 
                                                                      drugs=drugs, 
                                                                      fdr.cut.off=fdr, 
                                                                      nperm=100,
                                                                      top.ranked=top.ranked)
    }
    
    pdf(file.path(saveres, sprintf("validation_fdr_%s_%s.pdf", method, cell)), height=12.5, width=12.5)
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
##Plot expression based biomarkers effect size across all drugs 
estimatesScatterplot <-
  function(ccle.sig.rna, gdsc.sig.rna, cell, method, drugs, features, fdr.cut.off) {
    pdf(file.path(saveres, sprintf("effect_size_%s_%s.pdf", method, cell)), height=16, width=16)
    par(mfrow=c(4, 4))
    mycol <- RColorBrewer::brewer.pal(n=8, name="Set3")
    
    
    for (drug in drugs)
    {
      ###ccle
      top.ccle.sig.rna <- ccle.sig.rna[which((ccle.sig.rna[features, drug, "fdr"] < fdr.cut.off) | (gdsc.sig.rna[features, drug, "fdr"] < fdr.cut.off)), drug, ]
      if (is.null(dim(top.ccle.sig.rna))) {
        if(length(top.ccle.sig.rna) > 0){
          top.ccle.sig.rna <- as.matrix(t(top.ccle.sig.rna))
          rownames(top.ccle.sig.rna) <- names(which((ccle.sig.rna[features, drug, "fdr"] < fdr.cut.off) | (gdsc.sig.rna[features, drug, "fdr"] < fdr.cut.off)))
        }else{
          top.ccle.sig.rna <- matrix(NA, ncol=length(dimnames(ccle.sig.rna)[[3]]), dimnames=list(NA,dimnames(ccle.sig.rna)[[3]]))
        }
      }
      top.gdsc.sig.rna <- gdsc.sig.rna[which((ccle.sig.rna[features, drug, "fdr"] < fdr.cut.off) | (gdsc.sig.rna[features, drug, "fdr"] < fdr.cut.off)), drug, ]
      if (is.null(dim(top.gdsc.sig.rna))) {
        if(length(top.gdsc.sig.rna) > 0){
          top.gdsc.sig.rna <- as.matrix(t(top.gdsc.sig.rna))
          rownames(top.gdsc.sig.rna) <- names(which((ccle.sig.rna[features, drug, "fdr"] < fdr.cut.off) | (gdsc.sig.rna[features, drug, "fdr"] < fdr.cut.off)))
        }else{
          top.gdsc.sig.rna <- matrix(NA, ncol=length(dimnames(gdsc.sig.rna)[[3]]), dimnames=list(NA,dimnames(gdsc.sig.rna)[[3]]))
        }
      }
      features.drug <- rownames(top.ccle.sig.rna)
      
      point.col <- vector(length=length(features.drug), mode="character")
      names(point.col) <- features.drug
      point.col[features.drug] <- "black"
      point.col[which((top.ccle.sig.rna[features.drug, "fdr"] < fdr.cut.off) & (top.gdsc.sig.rna[features.drug, "fdr"] >= fdr.cut.off))] <- mycol[4]
      point.col[which((top.ccle.sig.rna[features.drug, "fdr"] >= fdr.cut.off) & (top.gdsc.sig.rna[features.drug, "fdr"] < fdr.cut.off))] <- mycol[5]
      ll <- min(top.ccle.sig.rna[,"estimate"], top.gdsc.sig.rna[,"estimate"], na.rm=T)
      if(is.na(ll) | ll == Inf | ll == -Inf) {ll <- -1}
      rr <- max(top.ccle.sig.rna[,"estimate"], top.gdsc.sig.rna[,"estimate"], na.rm=T)
      if(is.na(rr)| rr == Inf | rr == -Inf) {rr <- 1}
      rr <- max(rr , abs(ll))
      par(mar=c(5,5,2,2))
      plot(x=top.ccle.sig.rna[,"estimate"], y=top.gdsc.sig.rna[,"estimate"], ylim=c(-rr, rr), xlim=c(-rr, rr), xlab="CCLE", ylab="GDSC", main=drug, col=point.col, pch=20, cex.axis=1.5, cex.lab=1.5, cex.main=2)
      #plot(x=top.ccle.sig.rna[,"estimate"], y=top.gdsc.sig.rna[,"estimate"], ylim=c(ll, rr), xlim=c(ll, rr), xlab="CCLE", ylab="GDSC", main=drug, col=point.col, pch=20, cex.axis=1.5, cex.lab=1.5, cex.main=2)
      
      abline(0, 1, col="gray", lty=1, lwd=0.5 )
      abline(h=0, col="gray", lty=1, lwd=0.5)
      abline(v=0, col="gray", lty=1, lwd=0.5)
      
    }
    plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
    legend("center", legend=c("Both significant", "GDSC significant","CCLE significant"), col=c("black", mycol[c(5, 4)]), pch=15, bty="n" , cex=2)
    dev.off()
}


#################################################
## Determine the specificity of bimomarkers to the studies across drugs using Jaccard Index
biomarkersSpecificity <- 
  function(ccle.sig.rna, gdsc.sig.rna, cell, method, drugs, features, fdr.cut.off) {
    require(xtable)
    datset.specific.biomarkers <- matrix(NA, nrow=length(drugs), ncol=5)
    rownames(datset.specific.biomarkers) <- drugs
    colnames(datset.specific.biomarkers) <- c("# GDSC", "# CCLE", "% GDSC", "% CCLE", "% Both")

    for (drug in drugs)
    {
      
      top.ccle.sig.rna <- ccle.sig.rna[which((ccle.sig.rna[features, drug, "fdr"] < fdr.cut.off) | (gdsc.sig.rna[features, drug, "fdr"] < fdr.cut.off)), drug, ]
      if (is.null(dim(top.ccle.sig.rna))) {
        if(length(top.ccle.sig.rna) > 0){
          top.ccle.sig.rna <- as.matrix(t(top.ccle.sig.rna))
          rownames(top.ccle.sig.rna) <- names(which((ccle.sig.rna[features, drug, "fdr"] < fdr.cut.off) | (gdsc.sig.rna[features, drug, "fdr"] < fdr.cut.off)))
        }else{
          top.ccle.sig.rna <- matrix(NA, ncol=length(dimnames(ccle.sig.rna)[[3]]), dimnames=list(NA,dimnames(ccle.sig.rna)[[3]]))
        }
      }
      top.gdsc.sig.rna <- gdsc.sig.rna[which((ccle.sig.rna[features, drug, "fdr"] < fdr.cut.off) | (gdsc.sig.rna[features, drug, "fdr"] < fdr.cut.off)), drug, ]
      if (is.null(dim(top.gdsc.sig.rna))) {
        if(length(top.gdsc.sig.rna) > 0){
          top.gdsc.sig.rna <- as.matrix(t(top.gdsc.sig.rna))
          rownames(top.gdsc.sig.rna) <- names(which((ccle.sig.rna[features, drug, "fdr"] < fdr.cut.off) | (gdsc.sig.rna[features, drug, "fdr"] < fdr.cut.off)))
        }else{
          top.gdsc.sig.rna <- matrix(NA, ncol=length(dimnames(gdsc.sig.rna)[[3]]), dimnames=list(NA,dimnames(gdsc.sig.rna)[[3]]))
        }
      }
      rownames(top.ccle.sig.rna) <- paste(rownames(top.ccle.sig.rna), sign(top.ccle.sig.rna[,"estimate"]), sep="_")
      rownames(top.gdsc.sig.rna) <- paste(rownames(top.gdsc.sig.rna), sign(top.gdsc.sig.rna[,"estimate"]), sep="_")
      
      common <- intersect(rownames(top.ccle.sig.rna), rownames(top.gdsc.sig.rna))
      common <- length(which(top.ccle.sig.rna[common, "fdr"] < fdr.cut.off & top.gdsc.sig.rna[common, "fdr"] < fdr.cut.off))
      
      gdsc <- length(which(top.gdsc.sig.rna[,"fdr"] < fdr.cut.off))
      ccle <- length(which(top.ccle.sig.rna[,"fdr"] < fdr.cut.off))
      
      all <- gdsc + ccle - common
      
      datset.specific.biomarkers[drug, 1:5] <- c("# GDSC"=gdsc- common, 
                                              "# CCLE"=ccle- common, 
                                              "% GDSC"= round((gdsc- common)/all, digits=2) *100, 
                                              "% CCLE"=round((ccle- common)/all, digits=2) * 100,
                                              "% Both"=round(common/all, digits=2) * 100)
    }
    Jaccard <- datset.specific.biomarkers[, "% Both"]
    xtable::print.xtable(xtable::xtable(datset.specific.biomarkers, digits=0), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=sprintf("dataset_specific_biomarkers_%s_%s.tex", method, cell), append=FALSE)
    return(Jaccard)
}

###############################################
###check known biomarkers
knownBiomarkersCheck <-
  function(ccle.sig.rna, gdsc.sig.rna, ccle.sig.mutation, gdsc.sig.mutation, gdsc.sig.fusion, ccle.sig.cnv, gdsc.sig.cnv, method, cell)
  {
    
    known.biomarkers <- read.csv(file="known.biomarkers.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, na.strings=c("", " ", "NA"))
    known.biomarkers <- cbind(known.biomarkers, "GDSC effect size"=NA, "GDSC pvalue"=NA, "GDSC FDR"=NA, "CCLE effect size"=NA, "CCLE pvalue"=NA, "CCLE FDR"=NA, "Reproducible"=NA)
    known.biomarkers <- known.biomarkers[which(!is.na(known.biomarkers[ ,"type"])),]
    
    for(i in 1:nrow(known.biomarkers)) {
      if(!is.na(known.biomarkers[i ,"type"])) {
        if(known.biomarkers[i ,"type"] == "expression") {
          feature <- rownames(featureInfo(CCLE, "rna"))[which(featureInfo(CCLE, "rna")$Symbol == known.biomarkers[i ,"gene"])]
          known.biomarkers[i ,c("CCLE effect size", "CCLE pvalue", "CCLE FDR")] <- ccle.sig.rna[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
          known.biomarkers[i ,c("GDSC effect size", "GDSC pvalue", "GDSC FDR")] <- gdsc.sig.rna[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
          known.biomarkers[i, "Reproducible"] <- ifelse(known.biomarkers[i ,"CCLE pvalue"] < 0.05 & 
                                                          known.biomarkers[i, "GDSC pvalue"] < 0.05 & 
                                                          sign(known.biomarkers[i ,"CCLE effect size"]) == sign(known.biomarkers[i ,"GDSC effect size"]), "YES", "NO")
        }else if(known.biomarkers[i ,"type"] == "mutation") {
          feature <- known.biomarkers[i ,"gene"]
          known.biomarkers[i ,c("CCLE effect size", "CCLE pvalue", "CCLE FDR")] <- ccle.sig.mutation[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
          known.biomarkers[i ,c("GDSC effect size", "GDSC pvalue", "GDSC FDR")] <- gdsc.sig.mutation[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
          known.biomarkers[i, "Reproducible"] <- ifelse(known.biomarkers[i ,"CCLE pvalue"] < 0.05 & 
                                                          known.biomarkers[i, "GDSC pvalue"] < 0.05 & 
                                                          sign(known.biomarkers[i ,"CCLE effect size"]) == sign(known.biomarkers[i ,"GDSC effect size"]), "YES", "NO")
        }else if(known.biomarkers[i ,"type"] == "fusion") {
          feature <- known.biomarkers[i ,"gene"]
          known.biomarkers[i ,c("GDSC effect size", "GDSC pvalue", "GDSC FDR")] <- gdsc.sig.fusion[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
          #known.biomarkers[i, "Reproducible"] <- "NO"
        }else if(known.biomarkers[i ,"type"] == "amplification") {
          feature <- rownames(featureInfo(CCLE, "cnv"))[which(featureInfo(CCLE, "cnv")$Symbol == known.biomarkers[i ,"gene"])]
          known.biomarkers[i ,c("CCLE effect size", "CCLE pvalue", "CCLE FDR")] <- ccle.sig.cnv[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
          known.biomarkers[i ,c("GDSC effect size", "GDSC pvalue", "GDSC FDR")] <- gdsc.sig.cnv[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
          known.biomarkers[i, "Reproducible"] <- ifelse(known.biomarkers[i ,"CCLE pvalue"] < 0.05 & 
                                                          known.biomarkers[i, "GDSC pvalue"] < 0.05 & 
                                                          sign(known.biomarkers[i ,"CCLE effect size"]) == sign(known.biomarkers[i ,"GDSC effect size"]), "YES", "NO")
        }
      }
    }
    colnames(known.biomarkers)[1:3] <- capitalize(colnames(known.biomarkers)[1:3])
    xtable::print.xtable(xtable::xtable(known.biomarkers[, c(1:3, 7:8, 10:11, 13)], digits=c(0, 0, 0, 0, 2, -1, 2, -1, 0)), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=sprintf("known_biomarkers_%s_%s.tex", method, cell), append=FALSE)
    
    
}


#################################################
## Create an excel file for the statistics of expression based gene-drug associations in GDSC and CCLE
drugBasedBiomarkers <-
  function(ccle.sig.rna, gdsc.sig.rna, cell, method, drugs, features, cut.off) {
    require(WriteXLS)
    all.biomarkers <- list()
    for(drug in drugs) {
      ccle.biomarkers <- ccle.sig.rna[features, drug, ]
      colnames(ccle.biomarkers) <- paste0("CCLE_", colnames(ccle.biomarkers))
      
      gdsc.biomarkers <- gdsc.sig.rna[features, drug, ]
      colnames(gdsc.biomarkers) <- paste0("GDSC_", colnames(gdsc.biomarkers))
      
      biomarkers <- cbind("Symbol"=NA, gdsc.biomarkers, ccle.biomarkers, "Specificity"="Non significant")
      biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) < cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) < cut.off), "Specificity"] = "Both"
      biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) < cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) >= cut.off), "Specificity"] = "CCLE"
      biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) >= cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) < cut.off), "Specificity"] = "GDSC"
      #biomarkers[,"Symbol"] <- featureInfo(CCLE, "rna")[rownames(biomarkers), "Symbol"]
      biomarkers[,"Symbol"] <- rownames(biomarkers)
      all.biomarkers[[drug]] <- as.data.frame(biomarkers, stringsAsFactors=FALSE)
    }
    
    WriteXLS::WriteXLS("all.biomarkers", ExcelFileName=sprintf("all_biomarkers_%s_%s.xlsx", method, cell), row.names=TRUE)
    
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
      
      biomarkers <- rbind(biomarkers.rna, biomarkers.mutation, biomarkers.cnv)
      biomarkers[,"CCLE_fdr"] <- p.adjust(biomarkers[,"CCLE_pvalue"], method="fdr")
      biomarkers[,"GDSC_fdr"] <- p.adjust(biomarkers[,"GDSC_pvalue"], method="fdr")
      
      biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) < cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) < cut.off), "Specificity"] = "Both"
      biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) < cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) >= cut.off), "Specificity"] = "CCLE"
      biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) >= cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) < cut.off), "Specificity"] = "GDSC"
      biomarkers <- biomarkers[order(biomarkers[,"Specificity"]), ]
      
      all.biomarkers[[drug]] <- as.data.frame(biomarkers, stringsAsFactors=FALSE)
    }
    
    
    WriteXLS::WriteXLS("all.biomarkers", ExcelFileName=sprintf("all_biomarkers_%s.xlsx", method), row.names=TRUE)
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

