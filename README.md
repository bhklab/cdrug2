# Code for "Revisiting inconsistency in large pharmacogenomic studies"


## Full reproducibility of the analysis results
 
We will describe how to fully reproduce the figures and tables reported in the main manuscript. We automated the analysis pipeline so that minimal manual interaction is required to reproduce our results. To do this, one must simply:

  * Set up the software environment
  * Run the R scripts
  * Generate the Supplementary Information



The code and associated files are publicly available on GitHub: https://github.com/bhklab/cdrug2

## Set up the software environment

We developed and tested our analysis pipeline using R running on Linux and Mac OSX platforms.

To mimic our software environment the following R packages should be installed.


  * R version 3.3.1 (2016-06-21), x86_64-apple-darwin13.4.0
  * Base packages: base, datasets, graphics, grDevices, grid, methods, parallel, stats, utils
  * Other packages: Biobase 2.28.0, BiocGenerics 0.14.0, biomaRt 2.24.1, devtools 1.12.0, Formula 1.2-1, futile.logger 1.4.3, genefu 1.18.0, ggplot2 2.1.0, Hmisc 3.17-4, lattice 0.20-33, mclust 5.2, PharmacoGx 1.3.2, pROC 1.8, prodlim 1.5.7, psych 1.6.6, RColorBrewer 1.1-2, rJava 0.9-8, survcomp 1.18.0, survival 2.39-5, VennDiagram 1.6.17, xlsx 0.5.7, xlsxjars 0.6.1, xtable 1.8-2
  * Loaded via a namespace (and not attached): acepack 1.3-3.3, amap 0.8-14, AnnotationDbi 1.30.1, BiocInstaller 1.18.5, bitops 1.0-6, bootstrap 2015.2, caTools 1.17.1, celestial 1.3, chron 2.3-47, cluster 2.0.4, colorspace 1.2-6, curl 1.2, data.table 1.9.6, DBI 0.5, digest 0.6.10, downloader 0.4, foreign 0.8-66, futile.options 1.0.0, gdata 2.17.0, GenomeInfoDb 1.4.3, git2r 0.15.0, gplots 3.0.1, gridExtra 2.2.1, gtable 0.2.0, gtools 3.5.0, httr 1.2.1, igraph 1.0.1, IRanges 2.2.9, KernSmooth 2.23-15, lambda.r 1.1.9, latticeExtra 0.6-28, lava 1.4.4, limma 3.24.15, lsa 0.73.1, magicaxis 2.0.0, magrittr 1.5, mapproj 1.2-4, maps 3.1.1, marray 1.46.0, MASS 7.3-45, Matrix 1.2-7, memoise 1.0.0, mnormt 1.5-4, munsell 0.4.3, nnet 7.3-12, piano 1.8.2, plotrix 3.6-3, plyr 1.8.4, R6 2.1.3, RANN 2.5, Rcpp 0.12.6, RCurl 1.95-4.8, relations 0.6-6, rmeta 2.16, rpart 4.1-10, RSQLite 1.0.0, S4Vectors 0.6.6, scales 0.4.0, sets 1.0-16, slam 0.1-37, sm 2.2-5.4, SnowballC 0.5.1, splines 3.3.1, stats4 3.3.1, SuppDists 1.1-9.2, survivalROC 1.0.3, tools 3.3.1, withr 1.0.2, XML 3.98-1.4

All these packages are available on CRAN (http://cran.r-project.org) or Bioconductor (http://www.bioconductor.org).

Run the following commands in a R session to install all the required packages:

```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("VennDiagram", "Hmisc", "xtable", "RColorBrewer", "pROC", "Biobase", "genefu"
	"PharmacoGx", "xlsx"))
```

Note that *PharmacoGx* requires that several packages are installed. However, all dependencies are available from CRAN or Bioconductor.

Once the packages are installed, clone the *cdrug2* GitHub repository (https://github.com/bhklab/cdrug2) This should create a directory on the file system containing the following files:

  * cdrug2_analysis.R Script generating all the figures and tables reported in the manuscript.
  * cdrug2_foo.R Additional functions implemented specifically for the analysis and results visualization.


All the files required to run the automated analysis pipeline are now in place. It is worth noting that 500 MB storage for downloading CGP and CCLE PSets which is done Automatically when user runs \texttt{cdrug2\_analysis.R} script. 

Run the R scripts
Open a terminal window and go to the \texttt{cdrug} directory. You can easily run the analysis pipeline either in batch mode or in a R session. 

To run the full analysis pipeline in an R session, simply type the following command:

```R
nbcore <- 4
### to allocate four CPU cores for instance.
source("code/cdrug2_analysis_.R")
```

Key messages will be displayed to monitor the progress of the analysis.

The analysis pipeline was developed so that all intermediate analysis results are saved in the directory output. Therefore, in case of interruption, the pipeline will restart where it stopped.

## Generate the Supplementary Information

After completion of the analysis pipeline a directory output will be created to contain all the intermediate results, tables and figures reported in the main manuscript and this Supplementary Information.
