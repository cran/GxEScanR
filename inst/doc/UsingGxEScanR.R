## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(GxEScanR)
library(BinaryDosage)

## ----genetic-data-convert, eval = FALSE---------------------------------------
# vcffile <- system.file("extdata", "gendata.vcf.gz", package = "GxEScanR")
# 
# exampledir <- tempdir()
# bdosefile <- file.path(exampledir, "gendata.bdose")
# 
# BinaryDosage::vcftobd(vcffile = vcffile, bdose_file = bdosefile)
# bdinfo <- BinaryDosage::getbdinfo(bdosefile)

## ----genetic-data-------------------------------------------------------------
bdosefile <- system.file("extdata", "gendata.bdose", package = "GxEScanR")
bdinfo <- BinaryDosage::getbdinfo(bdosefile)
exampledir <- tempdir()

## ----subject-data-------------------------------------------------------------
subjectfile <- system.file("extdata", "subdata.rds", package = "GxEScanR")
subjectdata <- readRDS(subjectfile)
head(subjectdata)

## ----linear-model-------------------------------------------------------------
# Remove y_logistic from the subject data
lineardata <- subjectdata[, c(1, 2, 4, 5)]

# Keep only subjects with complete data
lineardata <- lineardata[complete.cases(lineardata), ]

# Keep only subjects with genetic data
lineardata <- lineardata[lineardata$subid %in% bdinfo$samples$sid, ]

# Fit the linear model
linearmodel <- glm(y_linear ~ x2 + x1, data = lineardata)

## ----logistic-model-----------------------------------------------------------
# Remove y_linear from the subject data
logisticdata <- subjectdata[, c(1, 3, 4, 5)]

# Keep only subjects with complete data
logisticdata <- logisticdata[complete.cases(logisticdata), ]

# Keep only subjects with genetic data
logisticdata <- logisticdata[logisticdata$subid %in% bdinfo$samples$sid, ]

# Fit the logistic model
logisticmodel <- glm(y_logistic ~ x2 + x1, family = binomial, data = logisticdata)

## ----allocate-memory----------------------------------------------------------
# Continuous outcome
linearmem <- gweis.mem(gemdl = linearmodel,
                       subids = lineardata$subid,
                       tests = c("bg_ge", "bg_gxe", "bgxe", "joint"))

# Dichotomous outcome
logisticmem <- gweis.mem(gemdl = logisticmodel,
                         subids = logisticdata$subid,
                         tests = c("bg_ge", "bg_gxe", "bgxe", "joint",
                                   "bg_eg", "bg_case", "bg_ctrl"))

## ----run-gweis----------------------------------------------------------------
snpindex <- 1:nrow(bdinfo$snps)

# Continuous outcome
linearresults <- file.path(exampledir, "linear.txt")
rungweis(gweismem = linearmem,
         bdinfo = bdinfo,
         snps = snpindex,
         outfilename = linearresults)

# Dichotomous outcome
logisticresults <- file.path(exampledir, "logistic.txt")
rungweis(gweismem = logisticmem,
         bdinfo = bdinfo,
         snps = snpindex,
         outfilename = logisticresults)

## ----read-results-------------------------------------------------------------
lineardf <- read.table(linearresults, header = TRUE, sep = "\t")
logisticdf <- read.table(logisticresults, header = TRUE, sep = "\t")

## ----snp-info-linear----------------------------------------------------------
knitr::kable(lineardf[1:3, 1:6], caption = "SNP information — continuous outcome")

## ----snp-info-logistic--------------------------------------------------------
knitr::kable(logisticdf[1:3, 1:8], caption = "SNP information — dichotomous outcome")

## ----model1-output------------------------------------------------------------
knitr::kable(lineardf[1:3, 7:8], caption = "Model 1 results")

## ----model2-output------------------------------------------------------------
knitr::kable(lineardf[1:3, 9:13], caption = "Model 2 results")

## ----model3-5-output----------------------------------------------------------
knitr::kable(logisticdf[1:3, 16:21], caption = "Models 3-5 results")

