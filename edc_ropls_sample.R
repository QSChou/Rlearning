library(ropls)

## UDF Function for preprocessing the data 
### Transform the Loss Rate into 'G'/'W'/'B' groups
RankLoss <- function (lr) {
  if (lr >= 0.0 & lr <= 0.1) return('G')
  if (lr > 0.1 & lr < 0.5) return('W')
  if (lr >= 0.5) return('B')
}

unit_test <- "T" #; unit_test <- "F"
## UDF:plot_loading
plot_loading <- function(fit_m,pcn) {
  if (missing(pcn)) pcn <- 1 
  else if(pcn <1 || pcn > dim(getLoadingMN(fit_m))[2]) return("pcn is out of range!")
  par(mfrow=c(1,1))
  par(mai = c(3,1,1,1) )
  barplot(sort(getLoadingMN(fit_m)[,pcn],decreasing = FALSE),las=2,col="green3")
  print("Max Loading in order (X)")
  print(head( sort(abs(getLoadingMN(fit_m)[,1]),decreasing = TRUE),3))
}
## UDF:plot score loading with colored by sample meta (Yi)
colrByYi <- function(fit_pca,sampleMeta,yi) { 
  lableX <- colnames(sampleMeta)
  plot(fit_pca, typeVc = "x-score", parAsColFcVn = sampleMeta[,lableX[yi]], parEllipsesL = TRUE)
}
#### Test Code: for plot_loading & colrByYi
if (unit_test == "T") {
  print("Unit Test for plot_loading; use iris as the test data~")
  attach(iris)
  iris_pca <- opls(iris[,1:4],predI =4)
  plot_loading(iris_pca,1)
  colrByYi(iris_pca,subset(iris,select = 5),1)
  detach(iris)
  remove(iris_pca)
}
## UDF:getPCMaxLoadV; arg pcn & maxN are optional, get the PC1 & 1st variable by default 
getPCMaxLoadV <- function(fit_m,pcn,maxN) {
  if (missing(pcn)) pcn <- 1 
  else if(pcn <1 || pcn > dim(getLoadingMN(fit_m))[2]) return("pcn is out of range!")
  if (missing(maxN)) maxN <- 1 
  else if(maxN <1 || maxN > dim(getLoadingMN(fit_m))[1]) return("maxN is out of range!")
  return(sort(abs(getLoadingMN(fit_m)[,pcn]),decreasing = TRUE)[1:maxN])
}
#### Test Code: for getPCMaxLoadV
#### if you are doing below test, try to examine the print result(PC1 & PC2() and check the location in loading plot
if (unit_test == "T") {
  print("Unit Test for getPCMaxLoadV; use iris as the test data~")
  attach(iris)
  iris_pca <- opls(iris[,1:4],predI =4)
  print("getPCMaxLoadV unit test by all default")
  getPCMaxLoadV(iris_pca,1,4)
  getPCMaxLoadV(iris_pca,2,4)
  detach(iris)
  remove(iris_pca)
}

### Do scaling on the data by center mean and unit variance
scale_d <- function(xMN) {
  xMeanVn <- apply(xMN, 2, function(colVn) mean(colVn, na.rm = TRUE))
  xSdVn <-  apply(xMN, 2, function(colVn) sd(colVn, na.rm = TRUE))
  xMN <- scale(xMN, center = xMeanVn, scale = xSdVn)
  return(xMN)
}
#### Test Code: for scale_d
if (unit_test == "T") {
  print("Unit Test for scale_d; use iris as the test data~")
  attach(iris)
  scale_iris <- scale_d(iris[,1:4])
  detach(iris)
  remove(scale_iris)
}

### Due to data matrix has not the col of the gls id, cannot use merge function to do join select; it's a bit of lousy to do
### two steps subset for extracting data matrix data by the criterion in the sample metadata.
### Step 1: sel y's id (glass id) according to given levels of the factor in Y (sample meta)
### Step 2: sel x's variables subset by the id in step1
selybyfaclevel <- function(smd,fac,lvl) {  ##smd -> data frame; fac -> col index of the factor; lvl -> lvl wish to extract
  return(subset(smd,smd[,fac] %in% lvl,select = 1))
}
selxbyid <- function(dm1,dm2rnv,vn) {
  res <- subset(dm1,row.names(dm1) %in% dm2rnv,select = vn) #Because data matrix has not the id col, use row.names to do filtering
  return(res)
}
### Test code for selybyfaclevel & selxbyid
if (unit_test == "T") {
  print("Unit Test for selybyfaclevel & selxbyid; use iris as the test data~")
  attach(iris)
  dataMatrix <- iris[,1:4]
  sampleMeta <- cbind(row.names(dataMatrix),as.factor(iris[,5]))
  res1 <- selybyfaclevel(sampleMeta,2,"2")
  res2 <- selxbyid(dataMatrix,res1,1)
  print(res1) ##Select result for id
  print(res2) ##Select result for X col
  detach(iris)
  remove(res1);remove(res2)
}

## plsda-fit with returning maxload variables' data case when exceeding the r2 threshold 
plsda_fit <- function(dataMatrix,sampleMeta,sampleyi,r2,maxLoadVN,sampleyj) {  #sampleyi is used for pls-da & sampleyj used for highlighting 
  possibleErr <- tryCatch(
    plsda<- opls(dataMatrix,sampleMeta[,sampleyi],orthoI = 0), ## orthoI = 0 to Do PLS-DA --
    error= function(e) e,
    finally = print(paste(sampleyi,colnames(sampleMeta)[sampleyi],sep = " var:"))
  )
  if(inherits(possibleErr,"error")) return("error")
  else if(getSummaryDF(plsda)[1] > r2 || getSummaryDF(plsda)[2] > r2) {  ## If R2X or R2Y > predefined R2 'threshold'
    plot_loading(plsda)
    if (missing(maxLoadVN)) maxLoadVN <- 1 #maxLoadVN define the number of variables return in order of loading.
    if (missing(sampleyj)) sampleyj <- sampleyi #If not the sampleyj assigned, used the sampleyi as the default.
    maxLoadV <- getPCMaxLoadV(plsda,1,maxLoadVN)
    print(paste("Max Load Var of PC1: ",maxLoadV, sep = ""))
    print(paste("plot the raw data of maxLoadV by different levels of ",levels(sampleMeta[,sampleyi]),sep = ":" ))
    YN <- colnames(sampleMeta)[sampleyi]
    YM <- colnames(sampleMeta)[sampleyj]
    XNi <- which(colnames(dataMatrix) %in% names(maxLoadV))
    print(paste("Class Y",YN,sep = ":")); print(paste("maxLoadV of X",names(maxLoadV), sep = ":"))
    allG <- cbind(sampleMeta[,sampleyi],subset(dataMatrix,select = XNi),sampleMeta[,sampleyj])
    return(list(plsda,allG))   
  } else return(plsda)
}
### graphD spec: row.name(=rowkey)/column1(=fixedY)/column2(=Y for DA)/column3~n(=X selected)
## scatter_plot will plot the input data by fixed format; the input data has the structure col1 -> factor; col2~colk -> numeric
scatter_plot <- function(graphD,markY) {
  graphD <-as.data.frame(graphD)
  if(is.factor(graphD[,1]) == "FALSE") graphD[,1] <- as.factor(graphD[,1])
  par(mfrow=c(1,1))
  par(mai = c(1,1,1,1) )
  colvl <- c("darkseagreen","darkolivegreen4","gold","deepskyblue2","darkseagreen2","gold2","darkseagreen4","gold4","deepskyblue",
             "darkolivegreen","deepskyblue4","darkolivegreen2")
  for (i in 2:dim(graphD)[2]-1) {
    print(i)
    u <- mean(graphD[,i])
    sdev <- sd(graphD[,i])
    row.names(graphD) <- 1:nrow(graphD) ## Avoidiing X error
    plot(row.names(graphD),graphD[,i],las=1,cex.lab=0.1)
    abline(h=u,col="gray60"); abline(h=u+3*sdev,col="gray60");abline(h=u-3*sdev,col="gray60")
    title(paste(colnames(graphD)[i],colnames(graphD)[1],sep = " by "))
    for (j in 1:length(levels(graphD[,1]))) {
      lvl <- levels(graphD[,1])[j]
      points(subset(row.names(graphD),graphD[,1] == lvl),
             subset(graphD[,i],graphD[,1] == lvl),
             col = colvl[j], pch = lvl)
    }
    
    lastcol <- dim(graphD)[2]
    points(subset(row.names(graphD),graphD[,lastcol] == markY),
           subset(graphD[,i],graphD[,lastcol] == markY),
           col = "red2", pch = markY)
  }
}

### Test code for plsda_fit & scatter_plot for maxLoadV variables 
if (unit_test == "T") {
  print("Unit Test for plsda_fit; use iris as the test data~")
  attach(iris)
  dataMatrix <- iris[,1:4]
  sampleMeta <- cbind(row.names(dataMatrix),as.factor(iris[,5]))
  #iris_plsda <- plsda_fit(dataMatrix,sampleMeta,2,0.5)
  iris_plsda <- plsda_fit(dataMatrix,sampleMeta,2,0.5,4)
  if(length(iris_plsda) == 2) scatter_plot(iris_plsda[2],"2")
  detach(iris)
  remove(dataMatrix);remove(sampleMeta);remove(iris_plsda)
}

#plsda_fit1 <- plsda_fit(edc_sample_dm_gb,gls_GB,14,0.1)

## dev.off() to reset the plot device
dev.off()
## Read the EDC sample data and then pre-processing the data  
setwd("~/data")
#edc_sample <- read.csv("simca-sample.csv") 
edc_sample <- read.csv("simca-sample.csv",fileEncoding = "latin1") 
## Seperate the edc_sammple into data matrix & sample meta
edc_sample_smd <- edc_sample[,1:19]                   ##The 1st~19 col assign to the sample meta

## Convert the Loss rate in % to numeric value
edc_sample_smd$LQBubbleLOSS <- as.numeric(sub("%","",edc_sample_smd$LQBubbleLOSS))/100

edc_sample_smd$RankLoss <- sapply(edc_sample_smd$LQBubbleLOSS,RankLoss)

## Convert the metadata into factor type
mdcol <- colnames(edc_sample_smd)
for (i in 1:length(mdcol)) {
  if (class(edc_sample_smd[,i]) != "factor") {
    edc_sample_smd[,i] <- as.factor(edc_sample_smd[,i])
    print(paste(mdcol[i],class(edc_sample_smd[,i]),sep = " Convert to "))
  }
}

edc_sample_dm <- edc_sample[,20:ncol(edc_sample)]    ## data matrix contain only numeric type
row.names(edc_sample_dm) <- edc_sample[,1]        ## suffix 'dm' denotes data matrix 
row.names(edc_sample_smd) <- edc_sample[,1]       ## suffix 'smd' denotes sample meta data
# -----Here's data ready --------------------------------------------------------------------------- 

## Auto fit - 1st try; all by default 
system.time(edc_sample_pca <- opls(edc_sample_dm))  ## only x is given, will fit PCA by degault 
#system.time(edc_sample_pca_svd <- opls(edc_sample_dm,predI = 4))  ## pca fitting, if not the 'algoC' assigned, use 'svd' as default algorithm
#edc_sample_pca_nipals <- opls(edc_sample_dm,algoC = "nipals",predI = 4) ## pca fitting, use 'nipals' as the algorithm

for (i in 2:dim(edc_sample_smd)[2]) {
  possibleErr <- tryCatch
  {
    print(i)
    Yi <- subset(edc_sample_smd,select = i)
    #if(is.finite(as.matrix(Yi))) 
    colrByYi(edc_sample_pca,subset(edc_sample_smd,select = i),1)
    error=function(e) e
  }
   if(inherits(possibleErr,"error")) next
}

edc_sample_pls <- opls(edc_sample_dm,edc_sample_smd[,20])
par(mfrow = c(1,1))
barplot(sort(getLoadingMN(edc_sample_opls)[,1],decreasing = TRUE))

## Seperate the data into only <= 10% and >= 50
###gls_GB <- subset(edc_sample_smd,RankLoss == c("G","B")) ##Without Re-factor
gls_GB <- as.data.frame(lapply(subset(edc_sample_smd,RankLoss == c("G","B")),function(x) if(is.factor(x)) factor(x) else x))
edc_sample_dm_gb <- subset(edc_sample_dm,row.names(edc_sample_dm) %in% gls_GB$"i..glass")

## 1st try by full set of X variables
edc_sample_pca_fv <- opls(edc_sample_dm_gb) # 1st try with full set of varialbes X

## simulate the PCA by eliminating the variable in order by the variance (min-> max)

edc_sample_dm_gb_sc <- scale_d(edc_sample_dm_gb)
edc_sample_dm_gb_sc_var <- sort(apply(edc_sample_dm_gb_sc,2,var),decreasing = FALSE)
#edc_sample_dm_gb_var <- sort(apply(edc_sample_dm_gb,2,var),decreasing = FALSE)
unlist <- names(edc_sample_dm_gb_sc_var)
candidate <- unlist
MaxR2 <- getSummaryDF(edc_sample_pca_fv)[1]

if (length(getPcaVarVn(edc_sample_pca_fv))>=2 & length(unlist) >=2) {
  i <- 1
  pcn <- length(getPcaVarVn(edc_sample_pca_fv))
  while (i <=length(unlist) && pcn >=2) {
    newlist <- unlist[i:length(unlist)]
    dt <- subset(edc_sample_dm_gb,select = newlist )
    edc_sample_pca <- opls(dt)
    i <- i+1
    pcn <- length(getPcaVarVn(edc_sample_pca))
    print(paste(i,pcn,sep = " run, PC Number: "))
    plot(edc_sample_pca, typeVc = "x-score", parAsColFcVn = gls_GB[,2], parEllipsesL = TRUE)
    barplot(sort(getLoadingMN(edc_sample_pca)[,1],decreasing = TRUE))
    print("Max Loading in order (X)")
    print(head( sort(abs(getLoadingMN(edc_sample_pca)[,1]),decreasing = TRUE),3))
    print("Min Loading in order (X)")
    print(head( sort(abs(getLoadingMN(edc_sample_pca)[,1]),decreasing = FALSE),3))
    
    if (getSummaryDF(edc_sample_pca)[1] > MaxR2) {
      MaxR2 <- getSummaryDF(edc_sample_pca)[1]
      candidate <- newlist
    }
  }
  candidate1 <- candidate  ## The candidate1 will be the reduction by single variable's variance
}

# simulate the PCA by eliminating the varialbe X in order by the PC1 loading (min->max)
pc1_x_loading_incr <- sort(abs(getLoadingMN(edc_sample_pca_fv)[,1]),decreasing = FALSE)
pc2_x_loading_incr <- sort(abs(getLoadingMN(edc_sample_pca_fv)[,2]),decreasing = FALSE)
unlist <- names(pc1_x_loading_incr)
MaxR2 <- getSummaryDF(edc_sample_pca_fv)[1]
if (length(getPcaVarVn(edc_sample_pca_fv))>=2 & length(unlist) >=2) {
  i <- 1
  pcn <- length(getPcaVarVn(edc_sample_pca_fv))
  while (i <=length(unlist) && pcn >=2) {
    newlist <- unlist[i:length(unlist)]
    dt <- subset(edc_sample_dm_gb,select = newlist )
    edc_sample_pca <- opls(dt)
    i <- i+1
    pcn <- length(getPcaVarVn(edc_sample_pca))
    print(paste(i,pcn,sep = " run, PC Number: "))
    plot(edc_sample_pca, typeVc = "x-score", parAsColFcVn = gls_GB[,2], parEllipsesL = TRUE)
    barplot(sort(getLoadingMN(edc_sample_pca)[,1],decreasing = TRUE))
    print("Max Loading in order (X)")
    print(head( sort(abs(getLoadingMN(edc_sample_pca)[,1]),decreasing = TRUE),3))
    print("Min Loading in order (X)")
    print(head( sort(abs(getLoadingMN(edc_sample_pca)[,1]),decreasing = FALSE),3))
    
    if (getSummaryDF(edc_sample_pca)[1] > MaxR2) {
      MaxR2 <- getSummaryDF(edc_sample_pca)[1]
      candidate <- newlist
    }
  }
  candidate2 <- candidate #The candidate2 will be the reduction by loading of PC1
}

## Get the MaxR2 model and color by different sample meta
dt <- subset(edc_sample_dm_gb,select = candidate1 )
edc_sample_pca <- opls(dt)
for (i in 1:dim(gls_GB)[2]) {
    possibleErr <- tryCatch(
    plot(edc_sample_pca, typeVc = "x-score", parAsColFcVn = gls_GB[,i], parEllipsesL = TRUE ),
    error= function(e) e, finally = print("Great!"))
  if(inherits(possibleErr,"error")) next
}

barplot(sort(getLoadingMN(edc_sample_pca)[,1],decreasing = TRUE))
print("Max Loading in order (X)")
print(head( sort(abs(getLoadingMN(edc_sample_pca)[,1]),decreasing = TRUE),3))
print("Min Loading in order (X)")
print(head( sort(abs(getLoadingMN(edc_sample_pca)[,1]),decreasing = FALSE),3))

## Get the MaxR2 model and color by different sample meta
dt <- subset(edc_sample_dm_gb,select = candidate2 )
edc_sample_pca <- opls(dt)
for (i in 1:dim(gls_GB)[2]) {
  possibleErr <- tryCatch(
    plot(edc_sample_pca, typeVc = "x-score", parAsColFcVn = gls_GB[,i], parEllipsesL = TRUE ),
    error= function(e) e, finally = print("Great!"))
  if(inherits(possibleErr,"error")) next
}



plot_loading(edc_sample_pca)
## Above 2 cases try to get the model with MaxR2 by comparing different X variables;
## Approach 1 will eliminating the variable with minimal variance step by step
## Approach 2 will eliminating the variable with minimal loading in PC1 from the previous fit
## The test show that get the candiate with MaxR2 is not necessary the best guess!
## Need to Consider different pumutation & index for evaluating the model quality..

dt <- subset(edc_sample_dm_gb) ## ,select = newlist )
edc_sample_opls <- opls(dt,gls_GB[,20],orthoI = NA) ## Do OPLS-DA fitting
plot_loading(edc_sample_opls)

selrawx <- selxbyid(edc_sample_dm_gb,selybyfaclevel(gls_GB,14,2)[,1],14)

## Using PLS-DA for factors with levels number > 2
row.names(gls_GB) <- gls_GB$i..glass

## fit y(14)
yi <- 2; th <- 0.25; yj <-20
maxit <- dim(gls_GB)[2]
#yi <- 15; maxit <- yi

colvl <- c("darkseagreen","darkolivegreen4","gold","deepskyblue2","darkseagreen2","gold2","darkseagreen4","gold4","deepskyblue",
           "darkolivegreen","deepskyblue4","darkolivegreen2")
while(yi <= maxit) {
  rpls <- plsda_fit(edc_sample_dm_gb,gls_GB,yi,th,1,yj) ## yi -> col i of y; th -> R2x,R2y threshold; 1 -> 1 var with max loading, yj -> 2nd col j of y for highlight
  if (length(rpls) >1) { 
    graphD <-as.data.frame(rpls[2]) 
    scatter_plot(graphD,"B")
  }
  yi <- yi+1
}
        


