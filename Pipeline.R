## load the package
library(GEOquery)
library(arrayQualityMetrics)
library(impute) 
library(limma)  

##step1 -- get data
# read in data from GEO database
gse <- getGEO("GSE8664",GSEMatrix=TRUE)[[1]] ## use Series Matirx format -- ExpressionSet class
class(gse) ## list 
show(gse)  ## show 
head(exprs(gse))  ## get value -- judge log2 transformed or not
head(featureNames(gse)) ## feature
sampleNames(gse)  ## sample name
phenoData(gse)    

## step2 -- quality assessment of raw data
arrayQualityMetrics(expressionset = gse,
                    outdir = "fig",
                    force = TRUE,
                    do.logtransform = F) ## set the out directory 


##!! the results show that the 15th sample may be a outlier -- do not use it for further analysis

## step3 -- get data matrix, select probeset with ORf
ORF                   <- featureData(gse)@data$ORF 
use_probe             <- which(is.na(ORF)==F & match(ORF,"",nomatch = 0)==0) ## only use probeset with ORF 
data.matrix           <- exprs(gse)[use_probe,] ## only use probes with ORF and delete the 15th sample
rownames(data.matrix) <- ORF[use_probe]
head(data.matrix) ## see data -- matrix format
rm(ORF,use_probe)
Med <- median(data.matrix,na.rm = T) 
if(Med > 16) data.matrix <- log2(data.matrix) ## judge the data -- log2 transformed or not
rm(Med)

## step4 -- fill in the NA value & normalization
na.length <- length(which(is.na(data.matrix)==T))
if(na.length > 0) data.matrix <- impute.knn(data.matrix)$data ## use impute.knn for normalization
rm(na.length)
data.matrix <- normalizeBetweenArrays(data.matrix) ## Normalizes expression intensities 
## so that the intensities or log-ratios 
##	have similar distributions across a set of arrays.
## combine the duplicate probeset
tmp <- aggregate(data.matrix,list(rownames(data.matrix)),median)
head(tmp) ## first column is the ORF name
data.matrix <- as.matrix(tmp[,-1])
rownames(data.matrix) <- tmp[,1] 
head(data.matrix) ##
dim(data.matrix)  ##  number of samples and genes
rm(tmp)

## quality assessment after normalization
gse1 <- gse ## better not overwrite the original file 
exprs(gse1) <- data.matrix
arrayQualityMetrics(expressionset = gse1,
                    outdir = "fig_norm",
                    force = TRUE,
                    do.logtransform = F) ## set the out directory 

