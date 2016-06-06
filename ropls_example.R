## Sample code from
## https://www.bioconductor.org/packages/3.3/bioc/vignettes/ropls/inst/doc/ropls.pdf

library(ropls)
data("sacurine")
names(sacurine)
attach(sacurine)
## The sacurine data contains list of 3 data sets;
## (1) dataMatrix       --> numeric data; kept original scale, with sample id as the row name & variates id as the col name
## (2) sampleMetadata   --> usually categorical type data of the samples are seperated into this data; used for color
## (3) variableMetadata  --> used for classification of variates
## The raw data is expected to be preprocessed into these three kinds of format before opls fitting 
strF(sacurine)
## run the 1st try by default.
## key settings by default;
###   (1) 'x'     ;Given input 'dataMatrix' will be taken as x (Numeric data frame or matrix, NAs are allowed) 
###   (2) 'y'     ;Because no 'y' is given, thus will be automatically be treated as 'PCA'
###   (3) 'predI' ;When PCA it will perform auto fit by default with max 10 componets and stop at the variance is less than the mean variance of all components.
###   (4) 'algoC' ;use 'svd' by default in case 'no missing value in x' for PCA case, otherwise (contains missing value), use 'nipals' instead.
###   (5) 'scaleC';the default setting will use mean-centering and unit variance scaling ('standard')   
sacurine.pca <- opls(dataMatrix)

## Use the value of the output to check some propertis of the fitted model
sacurine.pca$typeC                ### "PCA" 
sacurine.pca$descriptionMC        ### "183" samples with "109" x_variables and no near_zero_excluded_x_varialbes and no missing value
sacurine.pca$modelDF              ### automatically fit 8 PCs 
sacurine.pca$pcaVarVn             ### variance of each predic PCs.
##
genderFc <- sampleMetadata[,"gender"]
plot(sacurine.pca, typeVc = "x-score", parAsColFcVn = genderFc,parEllipsesL = TRUE)




## iris data can also be seperated into dataMatrix and sampleMetaData
iris_dm <- iris[,1:4]
iris_nMeta <- iris[,5]
## if 'predI' is not specified, then only one component is extracted by default because the variance of component 2 is less than the mean variance of all components
iris_pca <- opls(iris_dm,predI = 2)

## try to use the full explanation of the toall variance
iris_pca_fu <- opls(iris_dm, predI = 4)
iris_pca_fu$modelDF 
#R2X R2X(cum) Iter.
#p1 0.730    0.730     0
#p2 0.229    0.959     0
#p3 0.037    0.996     0
#p4 0.005    1.000     0




## Read the EDC sample data 

setwd("d:/userdata/R/data")

edc_sample <- read.csv("simca-sample.csv") 
## Seperate the edc_sammple into data matrix & sample meta
#edc_sample <- edc_sample[edc_sample$PICO_TFT.EQP=="PIPR06",]
edc_sample_md <- edc_sample[,1:19]                   ##The 1st~19 col assign to the sample meta

## Convert the Loss rate in % to numeric value
edc_sample_md$液晶氣泡LOSS <- as.numeric(sub("%","",edc_sample_md$液晶氣泡LOSS))/100
## Function to transform the Loss Rate
RankLoss <- function (lr) {
  if (lr >= 0.0 & lr <= 0.1) return('G')
  if (lr > 0.1 & lr <= 0.3) return('W')
  if (lr > 0.3) return('B')
}
edc_sample_md$RankLoss <- sapply(edc_sample_md$液晶氣泡LOSS,RankLoss)

## Convert the metadata into factor type
mdcol <- colnames(edc_sample_md)
for (i in 1:length(mdcol)) {
  if (class(edc_sample_md[,i]) != "factor") {
    edc_sample_md[,i] <- as.factor(edc_sample_md[,i])
    print(paste(mdcol[i],class(edc_sample_md[,i]),sep = " Convert to "))
  }
}


edc_sample_dm <- edc_sample[,20:ncol(edc_sample)]    ## data matrix contain only numeric type
row.names(edc_sample_dm) <- edc_sample[,1]
row.names(edc_sample_md) <- edc_sample[,1]

## Auto fit - default 
edc_sample_pca_svd <- opls(edc_sample_dm,predI = 4)

edc_sample_pca_nipals <- opls(edc_sample_dm,algoC = "nipals",predI = 4)

i <- 15
lableX <- colnames(edc_sample_md)
paste(lableX[i],levels(edc_sample_md[,i]),sep="")
plot(edc_sample_pca, typeVc = "x-score", parAsColFcVn = edc_sample_md[,lableX[i]], parEllipsesL = TRUE)
#plot(edc_sample_pca, typeVc =   , parAsColFcVn = edc_sample_md[,lableX[i]], parEllipsesL = TRUE)



for (i in 1:length(lableX)) {
  #possibleErr <- tryCatch
   #(
      plot(edc_sample_pca, typeVc = "x-score", parAsColFcVn = edc_sample_md[,lableX[i]], parEllipsesL = TRUE)
    #  error=function(e) e
   # )
 # if(inherits(possibleErr,"error")) next,
}













