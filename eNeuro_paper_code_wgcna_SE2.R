

# eNeuro paper code 


##########################################################################################################
## Eric Dammer - adapted code for WGCNA from Neelroop Parikshak, Vivek Swarup, and Divya Nandakumar
## SeyfriedLab&ProteomicsCorePipeline.R
##
## Applied to ROSMAP BA37 BULK Data from Herskowitz Lab at UAB - Evan Liu
##
## Goal: This code will carry out coexpression network analysis and systems biology for Protein Abundance
##       Also tested with alternate code for cleaning sparse RNA-Seq FPKM or TPM data (commented out)
##########################################################################################################

## CODE BLOCKS:
# Set up environment, load data,
# Clean up abundance matrix of rows/gene products with too many 0 FPKM/TPM or NA protein abundance values
# Perform TAMPOR robust median polish removal of unwanted batch covariance (preserving wanted biological variance)
# Outlier Removal using WGCNA z.k connectivity detection per each case-sample
# Nonparametric Bootstrap Regression of unwanted variance (covariance with confounding traits)

# WGCNA blockwiseModules network construction
# GlobalNetworkPlots and kMEtable Output
# FET to cell type specific symbol lists, modified to optionally adjust for cross-species symbol lookup inefficiency/loss
# ANOVA / DiffEx table generation -- list of ANOVAout dataframes for different subgroup comparisons, if necessary; otherwise, the list contains one data.frame
# Generate Volcano plots with WGCNA module color overlay (multiple options) 
# iGRAPHs (Multiple Toggle Options, e.g. BioGRID interactome overlap) // CONNECTIVITY PLOT
# One-Step GO-ELITE -- configurable user parameter section at top of code block
# Speakeasy2 code

## Set up environment
# (Always run the below code after loading a saved.image .Rdata)
#############################################################################################
## Set folders and libraries
rm(list=ls())
options(stringsAsFactors=FALSE)
rootdir <- "C:/Users/ehobby/Documents/EH_Emory41_Update/"
#"C:/Users/liuey/Documents/ROSMAP_WPCNA/" # This is the folder containing all of the analysis scripts, input, and output for this project
#functiondir <- "Code"
datadir <- "Input"
outputfigs <- "Figures"
outputtabs <- "Tables"

setwd(rootdir)

#functiondir = paste0(rootdir,functiondir,"/")
datadir = paste0(rootdir,datadir,"/")
outputfigs = paste0(rootdir,outputfigs,"/")
outputtabs = paste0(rootdir,outputtabs,"/")
dir.create(file.path(outputfigs)) #harmless warning given if already exists
dir.create(file.path(outputtabs))



library(WGCNA) # Network analysis package
library(NMF) # this package has a great annotated heatmap function - aheatmap
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(Cairo) # nicer graphics, anti-aliased, etc. --text from windows output PDFs using CairoPDF() function may not load in Illustrator, though -- so also use the pdf() standard output function when generating PDF figures
##Only for macs:
#CairoFonts(regular="Arial:style=Regular",bold="Arial:style=Bold",italic="Arial:style=Italic",bolditalic="Arial:style=Bold Italic,BoldItalic",symbol="Symbol")

#note: other libraries are loaded as needed in certain blocks but also listed below for completeness
library(reshape2)     #TAMPOR block
library(gridExtra)    #TAMPOR block
library(ggpubr)       #TAMPOR block
library("doParallel") #Bootstrap Regression block and GlobalNetworkPlots
library("biomaRt")    #Fisher Exact Test/list overlap block (enables cross-species lookup)
library("NMF")        #WGCNA GlobalNetworkPlots - eigengene heatmap and clustering
library("plotly")     #Volcano Figure Generation
library("stringr")    #GO-Elite block
library(boot)         #bootstrap regression block
library(tidyverse)
library(gplots) #for col2hex() fn (module plots at end)


##Check and change your input filenames below
inputAbundanceFile="post_TAMPOR_cleanDat.csv"     #abundance.CSV: rows are genes/proteins, columns are samples
inputTraitsFile="post_TAMPOR_traits.csv"                        #traits.CSV: rows are samples, columns are traits (as many as possible should be numerically coded)
#IMPORTANT:                                                    Make sure your inputTraitsFile has a "Group" column calling out (expected) different subsets of case-samples. Groups can be numerically coded if a single severity scale applies, but each group (number) should be used for at least 3 samples.

##Part of most output filenames specific to this project
projectFilesOutputTag="EH_41BULK_update"


## Load and clean the data
################################################################################################

cleanDat <- read.csv(file=paste0(datadir,inputAbundanceFile),header=TRUE,row.names=1)
range(cleanDat, na.rm=T) #get a sense of whether these are all positive >0 #will not be all >0 if already TAMPORed
# cleanDat<-log2(cleanDat) #only if input protein abundance data not already log2 transformed
# colnames(cleanDat)<-gsub("\\.","-",colnames(cleanDat))
GI <- rownames(cleanDat)


## Match up metadata and plot it for visualization - we want to know what variables are correlated with other variables
numericMeta <- read.csv(paste0(datadir,inputTraitsFile),header=TRUE,row.names=1)
rownames(numericMeta)==colnames(cleanDat) #sanity check -- are sample names in same order?
colnames(cleanDat) <- rownames(numericMeta)

cleanDat<-cleanDat[,match(rownames(numericMeta),colnames(cleanDat))] #cull to keep only samples we have traits for, and match the column (sample) order of cleanDat to the row order of numericMeta
#numericMeta <- numericMeta[match(colnames(cleanDat),rownames(numericMeta)),] #use this line instead of above if you have more samples in your traits file than you do in abundance data; trait sample (row) order will be matched to column names of cleanDat
rownames(numericMeta)==colnames(cleanDat) #sanity check -- are sample names in same order?




## If working on log2(protein abundance or ratio) with NA missing values; Enforce <50% missingness (1 less than half of cleanDat columns (or round down half if odd number of columns))
LThalfSamples<-length(colnames(cleanDat))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanDat)) %% 2)==1) { 0.5 } else { 1.0 }

#remove rows with >=50% missing values (only if there are some rows to be removed)
IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
temp2<-as.data.frame(cleanDat[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples),])
#handle condition if temp2 is for one row of cleandat (a vector instead of a DF)
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(cleanDat)[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples)]
}

if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples); rowsRemoved<-rownames(cleanDat)[IndexHighMissing]; cleanDat<-cleanDat[-IndexHighMissing,]; }

dim(cleanDat)
#8212   41

#write cleanDat abundance data to CSV file
write.csv(cleanDat,file=paste0(outputtabs,"cleanDat.log2-MissingData50pctControlled",projectFilesOutputTag,".csv"))

#save session data structures to Rdata for reloading intermittently
save.image(paste0("saved.image.",projectFilesOutputTag,".Rdata"))



##########################################################################################################
#=============================#
#  Check and Remove Outliers  #
#=============================#

if(!exists("numericMeta")) numericMeta<-traits
library(WGCNA)

sdout=3 #Z.k SD fold for outlier threshold
outliers.noOLremoval<-outliers.All<-vector()
cleanDat.noOLremoval<-cleanDat

for (repeated in 1:5) {
  normadj <- (0.5+0.5*bicor(cleanDat,use="pairwise.complete.obs")^2)
  
  ## Calculate connectivity
  netsummary <- fundamentalNetworkConcepts(normadj)
  ku <- netsummary$Connectivity
  z.ku <- ku-(mean(ku))/sqrt(var(ku))
  ## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
  outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
  print(paste0("There are ",sum(outliers)," outlier samples based on a bicor distance sample network connectivity standard deviation above ",sdout,".  [Round ",repeated,"]"))
  targets.All=numericMeta
  
  cleanDat <- cleanDat[,!outliers] 
  numericMeta <- targets <- targets.All[!outliers,]
  outliers.All<-c(outliers.All,outliers)
} #repeat 5 times

#All outliers removed
print(paste0("There are ",sum(outliers.All)," total outlier samples removed in ",repeated," iterations:"))
names(which(outliers.All))
outliersRemoved<-names(which(outliers.All))
#Note outliers as comment below, copied from R session.


## Enforce <50% missingness (1 less than half of cleanDat columns (or round down half if odd number of columns))
LThalfSamples<-length(colnames(cleanDat))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanDat)) %% 2)==1) { 0.5 } else { 1.0 }

## If operating on log2(FPKM) data, remove rows with >=50% originally 0 FPKM values (only if there are some rows to be removed)
#IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
#temp2<-data.frame(ThrowOut=apply(cleanDat,1,function(x) length(x[x==log2(0+0.05)])>LThalfSamples))
#cleanDat<-cleanDat[!temp2$ThrowOut,]
#dim(cleanDat) #still have x genes, now for y total samples

## If working on log2(protein abundance or ratio) with NA missing values; Enforce <50% missingness (1 less than half of cleanDat columns (or round down half if odd number of columns))
#remove rows with >=50% missing values (only if there are some rows to be removed)
IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
temp2<-as.data.frame(cleanDat[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples),])
#handle condition if temp2 is for one row of cleandat (a vector instead of a data frame)
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(cleanDat)[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples)]
}

if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples); rowsRemoved<-rownames(cleanDat)[IndexHighMissing]; cleanDat<-cleanDat[-IndexHighMissing,]; }

dim(cleanDat)


#Remove "Exclude" cases before regression
#Exclude may not be a homogeneous condition so remove
######################################
traits <- numericMeta

excludeVec <- grepl("Exclude", traits$Group)
cleanDat <- cleanDat[,!excludeVec]
traits <- traits[!excludeVec,]

rownames(traits) == colnames(cleanDat)
dim(cleanDat)
dim(traits)

na.coltest <- function(x) {
  w <- sapply(x, function(x)all(is.na(x)))
  vectorNA <- as.vector(w)
  return(vectorNA)
}

colnames(cleanDat) == rownames(traits)
bool_vectorNA <- na.coltest(cleanDat) #removes allNA columns
cleanDat <- cleanDat[,!bool_vectorNA]  
traits <- traits[!bool_vectorNA,] 
dim(cleanDat)
dim(traits)

LThalfSamples<-length(colnames(cleanDat))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanDat)) %% 2)==1) { 0.5 } else { 1.0 }

#remove rows with >=50% missing values (only if there are some rows to be removed)
IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
temp2<-as.data.frame(cleanDat[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples),])
#handle condition if temp2 is for one row of cleandat (a vector instead of a data frame)
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(cleanDat)[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples)]
}

if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples); rowsRemoved<-rownames(cleanDat)[IndexHighMissing]; cleanDat<-cleanDat[-IndexHighMissing,]; }

dim(cleanDat)




## Bootstrap regression set up EH edits

# # fix batch.channel - if batch channel is b01.126.AD we need it to be b01
# # Split the 'batch.channel' column into parts
# split_info <- do.call(rbind, strsplit(as.character(numericMeta$batch.channel), "\\."))
# 
# # Assign to new variables
# batch <- as.factor(split_info[, 1])
# channel <- as.factor(split_info[, 2])
# #diagnosis <- as.factor(split_info[, 3])
# 
# #sanity check
# table(batch)
# table(channel)
# #table(diagnosis)
# 
# #add back to numericMeta
# numericMeta$batch <- batch
# numericMeta$channel <- channel
# #numericMeta$diagnosis <- diagnosis


#######################
#Parallel Bootstrap Regression

cleanDat.unreg<-cleanDat
numericMeta <- traits
dim(numericMeta)

library("doParallel")
#when Eric is running at Emory (requires RSA public key-ssh, & manual run of shell script from command prompt to start server backend):  
# parallelThreads=30
# clusterLocal <- makeCluster(c(rep("haplotein.biochem.emory.edu",parallelThreads)), type = "SOCK", port=10191, user="edammer", rscript="/usr/bin/Rscript",rscript_args="OUT=/dev/null SNOWLIB=/usr/lib64/R/library",manual=FALSE)
# ##OR to run parallel processing threads for regression: :
parallelThreads=20 #max is number of processes that can run on your computer at one time
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")

registerDoParallel(clusterLocal)

##OR for no parallel processing skip all doParallel functions (much slower):
# parallelThreads=8


library(boot)
boot <- TRUE
numboot <- 1000
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows bootstrap function to select samples
  fit <- lm(formula, data=d)
  return(coef(fit))
}  

condition <- as.factor(numericMeta$Group)

# colnames(numericMeta)[grep("age", colnames(numericMeta))] = "Age"
# colnames(numericMeta)[grep("msex", colnames(numericMeta))] = "SEX"
# colnames(numericMeta)[grep("pmi", colnames(numericMeta))] = "PMI"

#regvars <- as.data.frame(cbind(condition, batch,PMI))

condition <- as.factor(numericMeta$Group)
batch     <- as.factor(numericMeta$Batch)   # preserves "b01", "b02", etc.
PMI       <- as.numeric(numericMeta$PMI)

regvars <- as.data.frame(cbind(condition, batch,PMI))

## Run the regression
normExpr.reg <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(cleanDat))
rownames(normExpr.reg) <- rownames(cleanDat)
colnames(normExpr.reg) <- colnames(cleanDat)
coefmat <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(regvars)+2) ## change this to ncol(regvars)+2 when condition has 2 levels if BOOT=TRUE, +1 if BOOT=FALSE

#another RNG seed set for reproducibility
set.seed(8675309);

if (parallelThreads > 1) {
  
  if (boot==TRUE) { #ORDINARY NONPARAMETRIC BOOTSTRAP
    set.seed(8675309)
    cat('[bootstrap-PARALLEL] Working on ORDINARY NONPARAMETRIC BOOTSTRAP regression with ', parallelThreads, ' threads over ', nrow(cleanDat), ' iterations.\n Estimated time to complete:', round(120/parallelThreads*nrow(cleanDat)/2736,1), ' minutes.\n') #intermediate progress printouts would not be visible in parallel mode
    coefmat <- foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% {
      set.seed(8675309)
      options(stringsAsFactors=FALSE)
      library(boot)
      thisexp <- as.numeric(cleanDat[i,])
      bs.results <- boot(data=data.frame(thisexp,regvars), statistic=bs,
                         R=numboot, formula=thisexp~ condition +batch+PMI)  ## run 1000 resamplings
      ## get the median - we can sometimes get NA values here... so let's exclude these - old code #bs.stats <- apply(bs.results$t,2,median) 
      bs.stats <- rep(NA,ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
      for (n in 1:ncol(bs.results$t)) {
        bs.stats[n] <- median(na.omit(bs.results$t[,n]))
      }
      bs.stats
      #cat('[bootstrap] Done for Protein ',i,'\n') #will not be visible
    }
    #    normExpr.reg <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(cleanDat))
    normExpr.reg <-foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% { cleanDat[i,]- coefmat[i,3]*regvars[,"batch"] - coefmat[i,4]*regvars[,"PMI"] }
  } else { #linear model regression; faster but incomplete regression of Age, Sex, PMI effects, SO NOT USED WITH boot=TRUE (requires changing coefmat matrix ncol to 1 less above)
    coefmat<-coefmat[,-ncol(coefmat)] #handles different column requirement for lm regression method
    for (i in 1:nrow(cleanDat)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(cleanDat[i,])~condition +age+sex+PMI,data=regvars) #+PMI
      ##datpred <- predict(object=lmmod1,newdata=regvars)
      coef <- coef(lmmod1)
      coefmat[i,] <- coef
      normExpr.reg[i,] <- coef[1] + coef[2]*regvars[,"condition"] + lmmod1$residuals ## The full data - the undesired covariates
      ## Also equivalent to <- thisexp - coef*var expression above
      #cat('Done for Protein ',i,'\n')
    }
  } #end parallel option -- Average run time estimate printed in console in minutes is calculated based on benchmark of a 2+ GHz intel Xeon with 30 threads 
} else {
  if (boot==TRUE) { #ORDINARY NONPARAMETRIC BOOTSTRAP
    for (i in 1:nrow(cleanDat)) {
      if (i%%1000 == 0) {print(i)}
      thisexp <- as.numeric(cleanDat[i,])
      bs.results <- boot(data=data.frame(thisexp,regvars), statistic=bs,
                         R=numboot, formula=thisexp~ condition +age+sex+PMI) ## run 1000 resamplings
      #                       R=numboot, formula=thisexp~ condition +age+sex+PMI) ## run 1000 resamplings
      ## get the median - we can sometimes get NA values here... so let's exclude these - old code #bs.stats <- apply(bs.results$t,2,median) 
      bs.stats <- rep(NA,ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
      for (n in 1:ncol(bs.results$t)) {
        bs.stats[n] <- median(na.omit(bs.results$t[,n]))
      }
      coefmat[i,] <- bs.stats
      normExpr.reg[i,] <- thisexp - bs.stats[3]*regvars[,"batch"]- bs.stats[4]*regvars[,"PMI"]
      cat('[bootstrap] Done for Protein ',i,'\n')
    }
  } else { #linear model regression; faster but NOT USED WITH boot=TRUE (requires changing coefmat matrix ncol to 1 less above)
    # coefmat<-coefmat[,-ncol(coefmat)] #handles different column requirement for lm regression method
    for (i in 1:nrow(cleanDat)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(cleanDat[i,])~condition +age+sex+PMI,data=regvars)
      ##datpred <- predict(object=lmmod1,newdata=regvars)
      coef <- coef(lmmod1)
      coefmat[i,] <- coef
      normExpr.reg[i,] <- coef[1] + coef[2]*regvars[,"condition"] + lmmod1$residuals ## The full data - the undesired covariates
      ## Also equivalent to <- thisexp - coef*var expression above
      cat('Done for Protein ',i,'\n')
    }
  } #end nonparallel option
}



#quantiles let us check if our regression worked

## Sanity Check -- Did regression do something unexpected to our abundance data?         EH this should be the same
quantile(cleanDat[,1],c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE)

#Emorys
# 0%       2.5%        25%        50%        75%      97.5%       100% 
# -3.1232752 -0.7816937 -0.1445080  0.0000000  0.1286251  0.5280247  3.8222000
#EH
# 0%       2.5%        25%        50%        75%      97.5%       100% 
# -3.1232752 -0.7816937 -0.1445080  0.0000000  0.1286251  0.5280247  3.8222000 

quantile(normExpr.reg[,1],c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE)                 #This should be slightly different 

#EH
# 0%         2.5%          25%          50%          75%        97.5%         100% 
# -3.066539173 -0.825317463 -0.160689774 -0.001135936  0.138690652  0.588191765  3.634146427 



##Overwrite cleanDat with regressed data
##(DO NOT RERUN OUT OF CONTEXT)
############################################
cleanDatReg<-normExpr.reg
rownames(cleanDat)<-rownames(cleanDat.unreg)
############################################

save.image(paste0("regressed.batchandPMI.saved.image.",projectFilesOutputTag,".Rdata"))  #overwrites
#importantly, now contains final cleanDat, cleanDat.unreg, final numericMeta, outliersRemoved

write.table(cleanDatReg,file=paste0(outputtabs,"/cleanDatReg.final_",projectFilesOutputTag,"_",nrow(cleanDatReg),"x",ncol(cleanDatReg),"_BootAgeSexPMIregr_GroupProtected.txt"),sep="\t") #check and apply changes to static tail of filename if necessary

colnames(cleanDatReg) == rownames(numericMeta)
bool_vectorNA <- na.coltest(cleanDatReg) # check again for NA columns
cleanDatReg <- cleanDatReg[,!bool_vectorNA]  
numericMeta <- numericMeta[!bool_vectorNA,] 
traits <- numericMeta

dim(cleanDatReg)
dim(numericMeta)

LThalfSamples<-length(colnames(cleanDatReg))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanDatReg)) %% 2)==1) { 0.5 } else { 1.0 }

#remove rows with >=50% missing values (only if there are some rows to be removed)
IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
temp2<-as.data.frame(cleanDatReg[which(rowSums(as.matrix(is.na(cleanDatReg)))>LThalfSamples),])
#handle condition if temp2 is for one row of cleandat (a vector instead of a data frame)
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(cleanDatReg)[which(rowSums(as.matrix(is.na(cleanDatReg)))>LThalfSamples)]
}

if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(cleanDatReg)))>LThalfSamples); rowsRemoved<-rownames(cleanDatReg)[IndexHighMissing]; cleanDatReg<-cleanDatReg[-IndexHighMissing,]; }

dim(cleanDatReg)


## WGCNA blockwiseModules (for signed bicor coexpression network)
##############################################################################

##Above code may have borrowed server for parallel processing of bootstrap regression -- return parallel processing to local workstation
library("doParallel")
library("snow")
# stopCluster(clusterLocal)
parallelThreads=20 #set to # of threads on your computer
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)
allowWGCNAThreads() #speeds the pickSoftThreshold function

powers <- seq(4,12,by=1)  #initial power check -- try to get SFT.R.sq to go > 0.80
sft <- pickSoftThreshold(t(cleanDatReg),blockSize=nrow(cleanDatReg)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")


#EH 04.24.2025
# Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
# 1     4    0.697 -5.18          0.920   722.0     700.0   1070
# 2     5    0.799 -4.05          0.950   443.0     420.0    800
# 3     6    0.857 -3.43          0.963   282.0     260.0    626
# 4     7    0.884 -3.00          0.962   185.0     166.0    506
# 5     8    0.907 -2.73          0.964   126.0     108.0    419
# 6     9    0.907 -2.58          0.959    87.6      71.4    353
# 7    10    0.904 -2.46          0.946    62.6      48.3    302
# 8    11    0.911 -2.33          0.946    45.8      33.3    261
# 9    12    0.918 -2.23          0.952    34.2      23.2    228

# Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1     4    0.747 -17.40          0.966  566.00    560.00  716.0
# 2     5    0.789 -11.80          0.973  312.00    306.00  452.0
# 3     6    0.822  -8.69          0.974  176.00    171.00  298.0
# 4     7    0.843  -6.99          0.973  101.00     96.90  204.0
# 5     8    0.864  -5.71          0.973   59.50     55.90  144.0
# 6     9    0.897  -4.73          0.978   35.80     32.80  105.0
# 7    10    0.926  -4.03          0.985   22.00     19.60   77.8
# 8    11    0.948  -3.53          0.988   13.80     11.80   59.1
# 9    12    0.960  -3.24          0.993    8.86      7.22   47.5


#plot initial SFT.R.sq vs. power curve
tableSFT<-sft[[2]]
plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R?")

# Plot the results
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(tableSFT[,1], -sign(tableSFT[,3])*tableSFT[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(tableSFT[,1], -sign(tableSFT[,3])*tableSFT[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(tableSFT[,1], tableSFT[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(tableSFT[,1], tableSFT[,5], labels=powers, cex=cex1,col="red")


allowWGCNAThreads() #speeds the pickSoftThreshold function

powers <- seq(5,12,by=0.5) #finer grained check over a honed range of power values.  No need to ever go above 30 (problem with data if so); higher power gives lower connectivity (k) and therefore more uncorrelated gene products in the network (grey with no coex module)
sft <- pickSoftThreshold(t(cleanDatReg),blockSize=nrow(cleanDatReg)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")

#EH 04.24.2025 same on 5.3.25 - chose power 8
# Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
# 1    5.0    0.799 -4.05          0.950   443.0     420.0    800
# 2    5.5    0.829 -3.71          0.957   351.0     329.0    704
# 3    6.0    0.857 -3.43          0.963   282.0     260.0    626  
# 4    6.5    0.869 -3.19          0.961   227.0     207.0    561
# 5    7.0    0.884 -3.00          0.962   185.0     166.0    506
# 6    7.5    0.896 -2.86          0.963   152.0     133.0    459
# 7    8.0    0.907 -2.73          0.964   126.0     108.0    419
# 8    8.5    0.898 -2.67          0.955   105.0      87.5    384
# 9    9.0    0.907 -2.58          0.959    87.6      71.4    353
# 10   9.5    0.912 -2.50          0.958    73.9      58.5    326   << power
# 11  10.0    0.904 -2.46          0.946    62.6      48.3    302
# 12  10.5    0.911 -2.38          0.949    53.4      40.0    280
# 13  11.0    0.911 -2.33          0.946    45.8      33.3    261
# 14  11.5    0.912 -2.29          0.947    39.5      27.7    243
# 15  12.0    0.918 -2.23          0.952    34.2      23.2    228

#Evans 
# Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1    5.0    0.789 -11.80          0.973  312.00    306.00  452.0
# 2    5.5    0.807 -10.00          0.975  233.00    228.00  365.0
# 3    6.0    0.822  -8.69          0.974  176.00    171.00  298.0
# 4    6.5    0.841  -7.66          0.976  133.00    128.00  246.0
# 5    7.0    0.843  -6.99          0.973  101.00     96.90  204.0
# 6    7.5    0.853  -6.29          0.972   77.40     73.50  171.0  << power
# 7    8.0    0.864  -5.71          0.973   59.50     55.90  144.0
# 8    8.5    0.873  -5.25          0.971   46.10     42.80  122.0
# 9    9.0    0.897  -4.73          0.978   35.80     32.80  105.0
# 10   9.5    0.907  -4.38          0.980   28.00     25.30   90.0
# 11  10.0    0.926  -4.03          0.985   22.00     19.60   77.8
# 12  10.5    0.940  -3.74          0.988   17.40     15.20   67.6
# 13  11.0    0.948  -3.53          0.988   13.80     11.80   59.1
# 14  11.5    0.956  -3.37          0.991   11.00      9.22   52.7
# 15  12.0    0.960  -3.24          0.993    8.86      7.22   47.5



#plot fine-grained results, looking for first power where SFT.R.sq has approached an asymptote
tableSFT<-sft[[2]]
plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R?")
#Notes on this data: looks choppy and SFT R^2 not improving in this range ... asymptote reached.

#choose power 10 elbow of SFT R? curve approaching asymptote near or ideally above 0.80
power=8


## Run an automated network analysis (ds=4 and mergeCutHeight=0.07, more liberal)
# choose parameters deepSplit and mergeCutHeight to get respectively more modules and more stringency sending more low connectivity genes to grey (not in modules).
allowWGCNAThreads(nThreads = 16)

net <- blockwiseModules(t(cleanDatReg),power=power,deepSplit=1,minModuleSize=30,
                        mergeCutHeight=0.07,TOMdenom="mean", #detectCutHeight=0.9999,                        #TOMdenom="mean" may get more small modules here.
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=nrow(cleanDatReg)+1000,reassignThresh=0.05)       #maxBlockSize always more than the number of rows in cleanDat
#blockwiseModules can take 30 min+ for large numbers of gene products/proteins (10000s of rows); much quicker for smaller proteomic data sets

net <- net.ds2.V2

# net.ds2.V2 is ds2 run on dec 9th, i think the other ds2 something is wrong

# net.ds1 <- net
# net.ds2.V2 <- net
# net.ds3 <- net
# net.ds4 <- net

# net.ds1.pwr6 <- net
# net.ds2.pwr6 <- net
# net.ds3.pwr6 <- net
# net.ds4.pwr6 <- net

# net.ds1.pwr10 <- net
# net.ds2.pwr10 <- net
# net.ds3.pwr10 <- net
# net.ds4.pwr10 <- net


nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))

##copy R session output;

# EH 5.03.2025 power 8
# Mnum           Color Size
# turquoise         M1       turquoise 1389
# blue              M2            blue  432
# brown             M3           brown  426
# yellow            M4          yellow  314
# green             M5           green  303
# red               M6             red  279
# black             M7           black  264
# pink              M8            pink  255
# magenta           M9         magenta  253
# purple           M10          purple  245
# greenyellow      M11     greenyellow  236
# tan              M12             tan  220
# salmon           M13          salmon  200
# cyan             M14            cyan  141
# midnightblue     M15    midnightblue  131
# lightcyan        M16       lightcyan  128
# grey60           M17          grey60  116
# lightgreen       M18      lightgreen  108
# lightyellow      M19     lightyellow  101
# royalblue        M20       royalblue   93
# darkred          M21         darkred   87
# darkgreen        M22       darkgreen   86
# darkturquoise    M23   darkturquoise   84
# darkgrey         M24        darkgrey   81
# orange           M25          orange   78
# darkorange       M26      darkorange   77
# white            M27           white   74
# skyblue          M28         skyblue   70
# saddlebrown      M29     saddlebrown   64
# steelblue        M30       steelblue   62
# paleturquoise    M31   paleturquoise   60
# violet           M32          violet   57
# darkolivegreen   M33  darkolivegreen   53
# darkmagenta      M34     darkmagenta   51
# sienna3          M35         sienna3   51
# yellowgreen      M36     yellowgreen   50
# skyblue3         M37        skyblue3   50
# plum1            M38           plum1   50
# orangered4       M39      orangered4   49
# mediumpurple3    M40   mediumpurple3   46
# lightsteelblue1  M41 lightsteelblue1   44
# lightcyan1       M42      lightcyan1   43
# ivory            M43           ivory   37




#clean dat 8218 41

## saved image of R session after running and finalizing blockwiseModules() function WGCNA output (now includes net data structure)
save.image(paste0("modules.saved.image.",projectFilesOutputTag,".Rdata"))  #overwrites

numericMeta1 <- numericMeta

limma::plotMDS(cleanDatReg, col=numericMeta$Color, main="After regressing out Batch/PMI and TAMPOR")

# Generate a boxplot for protein of interest.  (preferentially uses numericMeta to pair to cleanDat, but numericMeta exists only after removing outliers)
protein <- "UNC5B"

idx <- grepl(protein,rownames(cleanDatReg))
data_sub <- as.data.frame(cleanDatReg[idx,])
# data_sub <- data_sub[2,]
data_sub <- as.data.frame(t(data_sub))
# traits$Batch <- sampleIndex$batch[match(rownames(traits), rownames(numericMeta1))]

vectorGISdatasub <- grepl("GIS", rownames(data_sub))
colnames(data_sub) <- "Ratio"
data_sub$Group <- numericMeta$Group[match(rownames(data_sub),rownames(numericMeta1))]
# data_sub$Batch <- sampleIndex$batch[match(rownames(data_sub),rownames(numericMeta1))]

data_sub <- as.data.frame(data_sub[!vectorGISdatasub,])

plot <- ggplot(data_sub, aes(group = data_sub$Group, y = Ratio, fill = data_sub$Group)) + geom_boxplot(outlier.colour = "black", outlier.shape = 20, outlier.size = 1) + scale_fill_ghibli_d("MarnieLight1", -1)
plot + ggtitle(protein) +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 11, face = "bold"),
    axis.title.x = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y = element_text(color = "black", size = 11, face = "bold")
  ) 
## Output GlobalNetworkPlots and kMEtable
####################################################################################################################

numericMeta1 <- numericMeta

##############################
#EH added new spine traits and took out the incorrect ones
#EH removes unwanted columns from traits - numericMeta1
numericMeta1 <- numericMeta1[, -c(25:110)]

#only pick the spine traits that you are interested in merging
newTraits <- newTraits[, -c(2,3,5,8,10,11,12,17:22,31:34)]

#JH said get rid of these
# Total length
# Total spines
# X10
# Surface area
# Neck diameter
# Both ratios
# All percents
# Num dendrites

#E05-130   0.8498237   1.2989037   0.4129278

# Read with custom column names
newTraits <- read.csv(
  file = file.path(datadir, "CW_spines_PFC.csv"),
  header       = TRUE,
  fileEncoding = "UTF-8-BOM",
  check.names  = TRUE
)

library(dplyr)

newTraits <- newTraits %>%
  rename(SampleID = X)

newTraits <- newTraits[, -c(2:19,22,26,28:30,35:40,49:52)]

colnames(newTraits) <- make.names(colnames(newTraits), unique = TRUE)

#sanity check
str(newTraits)
head(colnames(newTraits), 10)

# Ensure SampleID formats match
numericMeta1$SampleID <- as.character(numericMeta1$SampleID)
newTraits$SampleID <- as.character(newTraits$SampleID)

# Match new spine traits to numericMeta1 by SampleID
matched_spineTraits <- newTraits[match(numericMeta1$SampleID, newTraits$SampleID), ]

# Check: all SampleIDs aligned?
if (!all(numericMeta1$SampleID == matched_spineTraits$SampleID)) {
  warning("⚠️ SampleID mismatch! Double check before merging.")
} else {
  message("✅ SampleIDs aligned safe to merge.")
}

# Drop SampleID from the trait data before merging
matched_spineTraits <- matched_spineTraits[, -which(colnames(matched_spineTraits) == "SampleID")]

# Merge into numericMeta1
numericMeta1 <- cbind(numericMeta1, matched_spineTraits)

#get rid of missing spine data
numericMeta1[numericMeta1 == "#DIV/0!"] <- NA



# EH i did this before merging 5.3.2025
# fix names in numericMeta1 if needed
# 1) The 36 old names
oldNames <- c(
  "Total.length..um.",   "Total.spines",             "Spine.Density.per.1um",
  "Spine.Density.per.10um","Backbone.Length.µm.",    "Volume.µm..",
  "Surface.Area.µm..",   "Head.Diameter.µm.",       "Neck.Diameter",
  "Head.Diameter.Neck.Diameter.µm.","Backbone.Length.Head.Diameter.µm",
  "Thin.spine.density",   "Stubby.spine.density",    "Mushroom.spine.density",
  "Filopodia.spine.density","Thin.spine..",          "Stubby.spine..",
  "Mushroom.spine..",     "Filopodia.spine..",       "Braak",
  "Dendrites",            "Length.of.Thin",          "Length.of.stubby",
  "Length.of.Mushroom",   "Length.of.Filopodia",     "Head.D...thin",
  "Head.D...stubby",      "Head.D...mushroom",       "Head.D...filopodia",
  "Neck.Diameter...thin", "Neck.Diameter...stubby",  "Neck.Diameter...Mush",
  "Neck.Diameter...f",     "Volume...T",              "Volume...S",
  "Volume...M",           "Volume...F"
)

# 2) Your desired new names, same length/order (fill in with whatever you like):
newNames <- c(
  "Total.Length",      "Total.Spines",            "Spine.Density",
  "X10","Backbone.Length",      "Volume",
  "Surface.Area",     "Head.Diameter",        "Neck.Diameter",
  "Head.Vs.Neck.Ratio",  "Backbone.To.Head.Ratio", "Thin.Spine.Density",
  "Stubby.Spine.Density",  "Mushroom.Spine.Density",   "Filopodia.Spine.Density",
  "Thin.Spine.Percent",       "Stubby.Spine.Percent",        "Mushroom.Spine.Percent",
  "Filopodia.Spine.Percent",  "Braak",             "Num.Dendrites",
  "Length.of.Thin.Spines",          "Length.of.Stubby.Spines",           "Length.of.Mushroom.Spines",
  "Length.of.Filopodia.Spines",     "Thin.Head.Diameter",           "Stubby.Head.Diameter",
  "Mush.Head.Diameter",    "Filopodia.Head.Diameter",      "Thin.Neck.Diameter",
  "Stubby.Neck.Diameter",      "Mush.Neck.Diameter",       "Filopodia.Neck.Diameter",
  "Volume.Thin",          "Volume.Stubby",           "Volume.Mushroom",
  "Volume.Filopodia"
)

# 3) Match and rename in numericMeta1
idx <- match(oldNames, names(newTraits))
names(newTraits)[idx] <- newNames

# 4) Verify
all(newNames %in% names(newTraits))  # should return TRUE


saveRDS(numericMeta1, file = "numericMeta1.rds")


################

FileBaseName=paste0(projectFilesOutputTag,power,"_MergeHeight_")
CairoPDF(file=paste0(rootdir,"01.30.2026.Global_plots_41Bulk.pdf"),width=16,height=12)


# # Open a larger plotting window
# dev.new(width = 16, height = 12)  # adjust width and height as needed

## Plot dendrogram with module colors and trait correlations
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDatReg), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta1)

numericIndices<-unique(c( which(!is.na(apply(numericMeta1,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta1,2,function(x) sum(as.numeric(x),na.rm=T)))==0) ))
#Warnings OK; This determines which traits are numeric and if forced to numeric values, non-NA values do not sum to 0

geneSignificance <- cor(sapply(numericMeta1[,numericIndices],as.numeric),t(cleanDatReg),use="pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta1)[numericIndices]
geneSigColors <- t(numbers2colors(t(geneSignificance),signed=TRUE,lim=c(-1,1),naColor="black"))
rownames(geneSigColors) <- colnames(numericMeta1)[numericIndices]


plotDendroAndColors(dendro=net$dendrograms[[1]],
                    colors=t(rbind(net$colors,geneSigColors)),
                    cex.dendroLabels=1.2,addGuide=FALSE,
                    dendroLabels=FALSE,
                    groupLabels = rep("", nrow(geneSigColors) + 1))



#groupLabels=c("Module Colors",colnames(numericMeta1)[numericIndices]))

## Plot eigengene dendrogram/heatmap - using bicor
tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL

plotEigengeneNetworks(tmpMEs, "Eigengene Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))



# #ANOVA
# numericMeta1$AD      <- 0
# numericMeta1$AsymAD  <- 0
# numericMeta1$Control <- 0
# 
# numericMeta1$AD     [ numericMeta1$Group == "AD"      ] <- 1
# numericMeta1$AsymAD [ numericMeta1$Group == "AsymAD"     ] <- 1
# numericMeta1$Control[ numericMeta1$Group == "Control" ] <- 1
# 
# Grouping <- character(nrow(numericMeta1))
# 
# Grouping[ numericMeta1$AD      == 1 ] <- "AD"
# Grouping[ numericMeta1$AsymAD  == 1 ] <- "AsymAD"
# Grouping[ numericMeta1$Control == 1 ] <- "Control"
# 
# Grouping <- factor(Grouping, levels = c("Control","AsymAD","AD"))


## new 
######################
## Find differences between Groups (as defined in Traits input file); Finalize Grouping of Samples for ANOVA

#Set a vector of strings that represent each sample in order, calling out each sample as a member of named groups (used by GlobalNetworkPlot boxplots, and later, ANOVA DiffEx)

# Create binary indicator columns
numericMeta1$AD      <- as.numeric(numericMeta1$Group == "AD")
numericMeta1$AsymAD  <- as.numeric(numericMeta1$Group == "AsymAD")
numericMeta1$Control <- as.numeric(numericMeta1$Group == "CT")  # CT = Control in the raw data

# Create the 'Grouping' factor with readable labels
Grouping <- character(nrow(numericMeta1))
Grouping[numericMeta1$AD == 1]      <- "AD"
Grouping[numericMeta1$AsymAD == 1]  <- "AsymAD"
Grouping[numericMeta1$Control == 1] <- "Control"

# Convert to factor with desired level order
Grouping <- factor(Grouping, levels = c("Control", "AsymAD", "AD"))



# This gets ANOVA (Kruskal-Wallis) nonparametric p-values for groupwise comparison of interest.
# look at numericMeta (traits data) and choose traits to use for linear model-determination of p value
head(numericMeta1)

# 1. Build your covariate data frame properly
regvars2 <- data.frame(
  AD      = as.factor(numericMeta1$AD),
  AsymAD  = as.numeric(numericMeta1$AsymAD),
  Control = as.numeric(numericMeta1$Control)
)

# 2. Compute one-way group p-values for each module eigengene
pvec <- sapply(seq_len(ncol(MEs)), function(i) {
  fit <- lm(MEs[, i] ~ AD, data = regvars2)
  f   <- summary(fit)$fstatistic
  pf(f[1], f[2], f[3], lower.tail = FALSE)
})
names(pvec) <- colnames(MEs)

# Inspect your results
head(pvec)
head(sort(pvec))

## should match these
#d2  5.03.2025
# darkgreen darkturquoise       magenta  midnightblue        salmon     lightcyan 
# 0.601497159   0.166752398   0.006907908   0.009415648   0.061151141   0.004737429 

#d4 5.03.2025
# lightgreen    turquoise         cyan midnightblue      magenta       salmon 
# 0.94697770   0.38526605   0.01072184   0.06083762   0.47142577   0.91066440 

# head(sort(pvec))
# red    darkgreen     darkgrey  greenyellow   lightcyan1  yellowgreen 
# 5.039833e-05 2.708332e-04 3.420530e-04 1.174826e-03 1.217107e-03 3.280813e-03 


# OLD 
#d2  5.01.2025
# midnightblue        black    darkgreen         pink    royalblue         cyan 
# 0.0002823662 0.7497950572 0.0331756359 0.5619077440 0.0533683103 0.0380083421 

#d2 upd
# green   lightgreen        brown       yellow midnightblue         cyan 
# 0.9872114945 0.4083586816 0.0001702289 0.4890324201 0.0298237801 0.2372072603 

#d2
# green          tan   lightgreen         blue        brown    royalblue 
# 0.4860135460 0.1885660078 0.2883007302 0.0426817366 0.0001788748 0.5958606283 

#d4
#  turquoise greenyellow      grey60     magenta         tan  orangered4 
#  0.315514402 0.005006144 0.071038995 0.579607308 0.024336486 0.002676361 

ApoE<- numericMeta1$ApoE
ApoE[numericMeta1$ApoE==22]<-"2/2"
ApoE[numericMeta1$ApoE==23]<-"2/3"
ApoE[numericMeta1$ApoE==24]<-"2/4"
ApoE[numericMeta1$ApoE==33]<-"3/3"
ApoE[numericMeta1$ApoE==34]<-"3/4"
ApoE[numericMeta1$ApoE==44]<-"4/4"

# drop empty apoe
numericMeta1$ApoE[numericMeta1$ApoE == ""] <- NA
numericMeta1$ApoE <- factor(numericMeta1$ApoE)

regvars3$ApoE <- numericMeta1$ApoE


## Find differences between APOE Risk Groups
regvars3 <- data.frame(as.factor(numericMeta1[,"ApoE"]),as.numeric(numericMeta1[,"AD"]),as.factor(numericMeta1[,"AsymAD"]),as.numeric(numericMeta1[,"Control"]))
colnames(regvars3) <- c("ApoE","AD","AsymAD","Control") ## data frame with covaraites in case we want to try multivariate regression
#aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
lm1 <- lm(data.matrix(MEs)~ApoE,data=regvars3)

pvec.Apoe <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.Apoe[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.Apoe) <- colnames(MEs)

head(pvec.Apoe)
head(sort(pvec.Apoe))

# 01.30.26
# head(pvec.Apoe)
# lightgreen    turquoise         cyan midnightblue      magenta       salmon 
# 0.07444765   0.14021615   0.50427188   0.61755890   0.43258865   0.24780579 

# head(sort(pvec.Apoe))
# lightyellow      violet   royalblue  lightgreen   turquoise      yellow 
# 0.02821897  0.05202414  0.05870586  0.07444765  0.14021615  0.21311128 
# 


## OLD
#d2
#  turquoise greenyellow      grey60     magenta         tan  orangered4 
#  0.17238305  0.32075923  0.50401965  0.46350051  0.36929353  0.04307129 

#d4
# turquoise  greenyellow       grey60      magenta          tan   orangered4 
# 8.404509e-05 3.851860e-05 2.389466e-02 4.179531e-01 7.645752e-02 1.152232e-02 



## this code does the exact same as above as confirmed by modules bellow
# 1) Recode APOE using the vector you already created
ApoE_char <- as.character(numericMeta1$ApoE)        # e.g. "E3/4", "E2/4"
# strip the leading "E" so you get "3/4", "2/4", etc.
ApoE_clean <- sub("^E", "", ApoE_char)

# 2) Build a factor in the exact order you care about:
ApoE_factor <- factor(
  ApoE_clean,
  levels = c("2/2","2/3","2/4","3/3","3/4","4/4")
)

# 3) Build a data.frame for the model — include any covariates if you like:
regvars_ApoE <- data.frame(
  ApoE = ApoE_factor,
  Age  = as.numeric(numericMeta1$Age),      # if you want to adjust for age
  Sex  = factor(numericMeta1$Sex),          # or keep as numeric 0/1
  PMI  = as.numeric(numericMeta1$PMI)
)

# 4) Fit one LM per module eigengene
pvec.Apoe <- sapply(seq_len(ncol(MEs)), function(i) {
  fit <- lm(MEs[, i] ~ ApoE, data = regvars_ApoE)
  f   <- summary(fit)$fstatistic
  pf(f[1], f[2], f[3], lower.tail = FALSE)
})
names(pvec.Apoe) <- colnames(MEs)

# 5) Quick check
print(head(pvec.Apoe))
head(sort(pvec.Apoe))


# print(head(pvec.Apoe))
# lightgreen    turquoise         cyan midnightblue      magenta       salmon 
# 0.07444765   0.14021615   0.50427188   0.61755890   0.43258865   0.24780579 
# 
# 
# head(sort(pvec.Apoe))
# lightyellow      violet   royalblue  lightgreen   turquoise      yellow 
# 0.02821897  0.05202414  0.05870586  0.07444765  0.14021615  0.21311128 




######################
## Get sigend kME values
kMEdat <- signedKME(t(cleanDatReg), tmpMEs, corFnc="bicor")


######################
## Plot eigengene-trait correlations - using p value of bicor for heatmap scale
library(RColorBrewer)
MEcors <- bicorAndPvalue(MEs,numericMeta1[,numericIndices])
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p


textMatrix = apply(moduleTraitCor,2,function(x) signif(x, 2))
#textMatrix = paste(signif(moduleTraitCor, 2), " (",
#  signif(moduleTraitPvalue, 1), ")", sep = "");
#dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3));

## Display the correlation values within a heatmap plot
cexy <- if(nModules>75) { 0.8 } else { 1 }
colvec <- rep("white",1500)
colvec[1:500] <- colorRampPalette(rev(brewer.pal(8,"BuPu")[2:8]))(500)
colvec[501:1000]<-colorRampPalette(c("white",brewer.pal(8,"BuPu")[2]))(3)[2] #interpolated color for 0.05-0.1 p
labeledHeatmap(Matrix = apply(moduleTraitPvalue,2,as.numeric),
               xLabels = colnames(numericMeta1)[numericIndices],
               yLabels = paste0("ME",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colvec,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y= cexy,
               zlim = c(0,0.15),
               main = paste("Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Student correlation p value"),
               cex.main=0.8)


######################
## Plot eigengene-trait heatmap custom - using bicor color scale

numericMetaCustom<-numericMeta1[,numericIndices]
MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

moduleTraitPvalue<-signif(moduleTraitPvalue, 1)
moduleTraitPvalue[moduleTraitPvalue > as.numeric(0.05)]<-as.character("")

textMatrix = moduleTraitPvalue; #paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#textMatrix = gsub("()", "", textMatrix,fixed=TRUE)

labelMat<-matrix(nrow=(length(names(MEs))), ncol=2,data=c(rep(1:(length(names(MEs)))),labels2colors(1:(length(names(MEs))))))
labelMat<-labelMat[match(names(MEs),labelMat[,2]),]
for (i in 1:(length(names(MEs)))) { labelMat[i,1]<-paste("M",labelMat[i,1],sep="") }
for (i in 1:length(names(MEs))) { labelMat[i,2]<-paste("ME",labelMat[i,2],sep="") }

#rowMin(moduleTraitPvalue) # if we want to resort rows by min P value in the row
xlabAngle <- if(nModules>75) { 90 } else { 45 }

par(mar=c(16, 12, 3, 3) )
par(mfrow=c(1,1))

bw<-colorRampPalette(c("#0058CC", "white"))
wr<-colorRampPalette(c("white", "#CC3300"))

colvec<-c(bw(50),wr(50))

labeledHeatmap(Matrix = t(moduleTraitCor)[,],
               yLabels = colnames(numericMetaCustom),
               xLabels = labelMat[,2],
               xSymbols = labelMat[,1],
               xColorLabels=TRUE,
               colors = colvec,
               textMatrix = t(textMatrix)[,],
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x = cexy,
               xLabelsAngle = xlabAngle,
               verticalSeparator.x=c(rep(c(1:length(colnames(MEs))),as.numeric(ncol(MEs)))),
               verticalSeparator.col = 1,
               verticalSeparator.lty = 1,
               verticalSeparator.lwd = 1,
               verticalSeparator.ext = 0,
               horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
               horizontalSeparator.col = 1,
               horizontalSeparator.lty = 1,
               horizontalSeparator.lwd = 1,
               horizontalSeparator.ext = 0,
               zlim = c(-1,1),
               main = "Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
               cex.main=0.8)


# # turn your grouping column into an ordered factor
# Group <- factor(numericMeta1$Group,
#                 levels = c("CT", "AsymAD", "AD"))
# 
# # then build your metadata frame
# metdat <- data.frame(
#   Group  = Group,
#   Age    = as.numeric(numericMeta1$Age),
#   Gender = Gender
# )
# 
# # quick check
# str(metdat)
# 


## Plot annotated heatmap - annotate all the metadata, plot the eigengenes!
# This is where we will first use the Grouping vector of string group descriptions we set above.
toplot <- MEs

colnames(toplot) <- colnames(MEs)
rownames(toplot) <- rownames(MEs)
toplot <- t(toplot)

pvec <- pvec[match(names(pvec),rownames(toplot))]
#rownames(toplot) <- paste(rownames(toplot),"\np = ",signif(pvec,2),sep="")
rownames(toplot) <- paste(orderedModules[match(colnames(MEs),orderedModules[,2]),1]," ",rownames(toplot),"p=",signif(pvec,2),sep="")

# add any traits of interest you want to be in the legend
Gender=as.numeric(numericMeta1$Sex)
Gender[Gender==0]<-"Female"
Gender[Gender==1]<-"Male"
metdat=data.frame(Group=Grouping,Age=as.numeric(numericMeta1$Age), Gender=Gender)

# set colors for the traits in the legend

dev.new(width = 10, height = 10)           #opens in bigger ploting window

heatmapLegendColors=list('Group'=c("dodgerblue","goldenrod","seagreen3"), #,"hotpink","purple"),
                         'Age'=c("white","darkgreen"), #young to old
                         'Gender'=c("pink","dodgerblue"), #F, M
                         'Modules'=sort(colnames(MEs)))

library(NMF)
par(mfrow=c(1,1))
aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES IN ORIGINAL, e.g. BATCH OR REGION ORDER",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=NA) ## Do not cluster columns - keep given order

aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES CLUSTERED",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE,Colv=TRUE) ## Cluster columns


dev.off()

# library(dplyr)
# 
# numericMeta1 <- numericMeta1 %>%
#   rename(
#     FontralDP = FrontalDP,
#   )



######################################
## Change the below code in the for loop using the following session output
library(gplots) #for col2hex() fn
library(beeswarm)

## Get module-trait bicor correlations (append to verboseScatterplot title below)
#numericMetaCustom<-numericMeta[,numericIndices]
MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

#These are your numerically coded traits:
colnames(numericMeta1)[numericIndices] #choose traits for correlation scatterplots (verboseScatterplot functions below)

#These are your ANOVA sample groups and the number of samples in each
table(Grouping) #alphabetically ordered, you choose the order of groups in the boxplot function by typing them in

## Make changes after checking output on console for the above 2 lines
par(mfrow=c(4,6))
par(mar=c(4.5,6,4.5,1.5))

for (i in 1:(nrow(toplot))) {  # grey already excluded, no -1
  titlecolor<-if(signif(pvec,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(Grouping,c("Control","AsymAD","AD")),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\np = ",signif(pvec,2)[i]),xlab=NULL,las=2,col.main=titlecolor)  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(Grouping,c("Control","AsymAD","AD")),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism ; #bg=goldenrod #DDA43B "#DDA43B99"
  
  # titlecolor<-if(signif(pvec.Apoe,2)[i] <0.05) { "red" } else { "black" }
  # boxplot(toplot[i,]~factor(ApoE,c("2/3","3/3","3/4","4/4")),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nK-W p = ",signif(pvec.Apoe,2)[i]),xlab="APOE genotype",las=2,col.main=titlecolor)  #no outliers: ,outline=FALSE)
  # transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  # beeswarm(toplot[i,]~factor(ApoE,c("2/3","3/3","3/4","4/4")),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism ; #bg=goldenrod #DDA43B "#DDA43B99"
  
  verboseScatterplot(x=numericMeta1[,"CERAD"],y=toplot[i,],xlab="CERAD Score",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"CERAD"],2),", p=",signif(moduleTraitPvalue[i,"CERAD"],2),"\n"),col.main=if(moduleTraitPvalue[i,"CERAD"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta1[,"BRAAK"],y=toplot[i,],xlab="Braak Score",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"BRAAK"],2),", p=",signif(moduleTraitPvalue[i,"BRAAK"],2),"\n"),col.main=if(moduleTraitPvalue[i,"BRAAK"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta1[,"FrontalNP"],y=toplot[i,],xlab="FrontalNP",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FrontalNP"],2),", p=",signif(moduleTraitPvalue[i,"FrontalNP"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FrontalNP"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta1[,"FontralDP"],y=toplot[i,],xlab="FrontalDP",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FontralDP"],2),", p=",signif(moduleTraitPvalue[i,"FontralDP"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FontralDP"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta1[,"FrontalNFT"],y=toplot[i,],xlab="FrontalNFT",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FrontalNFT"],2),", p=",signif(moduleTraitPvalue[i,"FrontalNFT"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FrontalNFT"]<0.05) { "red" } else { "black" })
  
  verboseScatterplot(x=numericMeta1[,"MMSE"],y=toplot[i,],xlab="MMSE",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"MMSE"],2),", p=",signif(moduleTraitPvalue[i,"MMSE"],2),"\n"),col.main=if(moduleTraitPvalue[i,"MMSE"]<0.05) { "red" } else { "black" })
  #verboseScatterplot(x=numericMeta[,"ABC"],y=toplot[i,],xlab="ABC",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"ABC"],2),", p=",signif(moduleTraitPvalue[i,"ABC"],2),"\n"),col.main=if(moduleTraitPvalue[i,"ABC"]<0.05) { "red" } else { "black" })
  
}

dev.off()









# Open PDF device
pdf("01.30.26.updated_module_plots.pdf", width = 14, height = 10)  # adjust width/height as needed

# Set plotting layout
par(mfrow=c(4,6))
par(mar=c(4.5,6,4.5,1.5))

# Your existing for loop
for (i in 1:(nrow(toplot))) {  # grey already excluded, no -1
  titlecolor <- if(signif(pvec,2)[i] < 0.05) { "red" } else { "black" }
  
  # Boxplot + beeswarm
  boxplot(
    toplot[i,] ~ factor(Grouping, c("Control","AsymAD","AD")),
    col = colnames(MEs)[i],
    ylab = "Eigenprotein Value",
    main = paste0(orderedModules[match(colnames(MEs)[i], orderedModules[,2]),1], " ", colnames(MEs)[i], "\np = ", signif(pvec,2)[i]),
    xlab = NULL,
    las = 2,
    col.main = titlecolor
  )
  
  transcol = paste0(col2hex(colnames(MEs)[i]), "99")
  beeswarm(
    toplot[i,] ~ factor(Grouping, c("Control","AsymAD","AD")),
    method = "swarm",
    add = TRUE,
    corralWidth = 0.5,
    vertical = TRUE,
    pch = 21,
    bg = transcol,
    col = "black",
    cex = 0.8,
    corral = "gutter"
  )
  
  # Verbose scatterplots (your existing calls)
  verboseScatterplot(x=numericMeta1[,"CERAD"],y=toplot[i,],xlab="CERAD Score",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"CERAD"],2),", p=",signif(moduleTraitPvalue[i,"CERAD"],2),"\n"),col.main=if(moduleTraitPvalue[i,"CERAD"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta1[,"BRAAK"],y=toplot[i,],xlab="Braak Score",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"BRAAK"],2),", p=",signif(moduleTraitPvalue[i,"BRAAK"],2),"\n"),col.main=if(moduleTraitPvalue[i,"BRAAK"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta1[,"FrontalNP"],y=toplot[i,],xlab="FrontalNP",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FrontalNP"],2),", p=",signif(moduleTraitPvalue[i,"FrontalNP"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FrontalNP"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta1[,"FontralDP"],y=toplot[i,],xlab="FrontalDP",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FontralDP"],2),", p=",signif(moduleTraitPvalue[i,"FontralDP"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FontralDP"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta1[,"FrontalNFT"],y=toplot[i,],xlab="FrontalNFT",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FrontalNFT"],2),", p=",signif(moduleTraitPvalue[i,"FrontalNFT"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FrontalNFT"]<0.05) { "red" } else { "black" })
  
  verboseScatterplot(x=numericMeta1[,"MMSE"],y=toplot[i,],xlab="MMSE",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"MMSE"],2),", p=",signif(moduleTraitPvalue[i,"MMSE"],2),"\n"),col.main=if(moduleTraitPvalue[i,"MMSE"]<0.05) { "red" } else { "black" })
  
  # Repeat for the other numeric traits (BRAAK, FrontalNP, etc.)
  # ... your existing verboseScatterplot calls ...
}

# Close PDF device
dev.off()



# clean spine names

numericMeta1 <- numericMeta1 %>%
  rename(
    FrontalDP = FontralDP,
    total_length_um = Total.length..um.,
    total_spines = Total.spines,
    spine_density = Spine.Density.per.10um,
    spine_length = Backbone.Length.µm.,
    volume = Volume.µm..,
    head_diameter = Head.Diameter.µm.,
    thin_spine_density = Thin.spine.density,
    stubby_spine_density = Stubby.spine.density,
    mushroom_spine_density = Mushroom.spine.density,
    filopodia_spine_density = Filopodia.spine.density,
    thin_spine_length = Length.of.Thin,
    stubby_spine_length = Length.of.stubby,
    mushroom_spine_length = Length.of.Mushroom,
    filopodia_spine_length = Length.of.Filopodia,
    thin_head_diameter = Head.D...thin,
    stubby_head_diameter = Head.D...stubby,
    mushroom_head_diameter = Head.D...mushroom,
    filopodia_head_diameter = Head.D...filopodia,
    thin_volume = Volume...T,
    stubby_volume = Volume...S,
    mushroom_volume = Volume...M,
    filopodia_volume = Volume...F
  )







numericIndices <- match(allTraits, colnames(numericMeta1))



library(gplots) #for col2hex() fn
library(WGCNA)
library(Cairo)
library(beeswarm)

# Define all spine traits to loop over
allTraits <- c(
  "total_length_um", "total_spines", "spine_density", "spine_length",
  "volume", "head_diameter",
  "thin_spine_density", "stubby_spine_density", "mushroom_spine_density", "filopodia_spine_density",
  "thin_spine_length", "stubby_spine_length", "mushroom_spine_length", "filopodia_spine_length",
  "thin_head_diameter", "stubby_head_diameter", "mushroom_head_diameter", "filopodia_head_diameter",
  "thin_volume", "stubby_volume", "mushroom_volume", "filopodia_volume"
)

#sanity check
setdiff(allTraits, colnames(numericMeta1))

# Open PDF device
# CairoPDF(file = paste0(rootdir, "01.30.26_spine_plots.pdf"), width = 16, height = 12)

# Set up layout and margins
par(mfrow = c(4, 6))
par(mar = c(4.5, 6, 4.5, 1.5))

# Loop through each module
for (i in seq_len(nrow(toplot))) {
  y <- toplot[i, ]
  moduleCol <- colnames(MEs)[i]
  moduleName <- rownames(toplot)[i]
  
  # Panel 1: Boxplot + Beeswarm by group
  boxplot(y ~ factor(Grouping, c("Control", "AsymAD", "AD")),
          col = moduleCol, ylab = "Eigengene Value",
          main = moduleName, xlab = NULL, las = 2)
  beeswarm(y ~ factor(Grouping, c("Control", "AsymAD", "AD")),
           add = TRUE, method = "swarm", pch = 21,
           bg = paste0(col2hex(moduleCol), "99"), corral = "gutter")
  
  # Panels 2–N: Trait correlations
  for (trait in allTraits) {
    if (!trait %in% colnames(numericMeta1)) next
    
    x_raw <- numericMeta1[[trait]]
    x <- suppressWarnings(as.numeric(as.character(x_raw)))
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 2) next
    
    verboseScatterplot(
      x = x[ok],
      y = y[ok],
      xlab = trait,
      ylab = "Eigenprotein",
      abline = TRUE,
      cex.axis = 1, cex.lab = 1, cex = 1,
      col = "black", bg = moduleCol, pch = 21,
      main = paste0(
        "bicor=", if (is.numeric(moduleTraitCor[i, trait]) && is.finite(moduleTraitCor[i, trait])) signif(moduleTraitCor[i, trait], 2) else "NA",
        ", p=", if (is.numeric(moduleTraitPvalue[i, trait]) && is.finite(moduleTraitPvalue[i, trait])) signif(moduleTraitPvalue[i, trait], 2) else "NA"
      ),
      col.main = if (is.numeric(moduleTraitPvalue[i, trait]) && moduleTraitPvalue[i, trait] < 0.05) "red" else "black"
    )
  }
  
  frame()  # starts a new page of 24 panels if layout is full
}

# Close the PDF device
dev.off()





for (i in 1:(nrow(toplot))) {  # grey already excluded, no -1
  verboseScatterplot(
    x = numericMeta1$spine_density, 
    y = toplot[i,],
    xlab = "Spine Density 10um",
    ylab = "Eigenprotein",
    abline = TRUE,
    cex.axis = 1, cex.lab = 1, cex = 1,
    col = "black",
    bg = colnames(MEs)[i],
    pch = 21,
    main = paste0(
      "bicor=", signif(moduleTraitCor[i, "spine_density"], 2),
      ", p=", signif(moduleTraitPvalue[i, "spine_density"], 2), "\n"
    ),
    col.main = if (moduleTraitPvalue[i, "spine_density"] < 0.05) { "red" } else { "black" }
  )
}




library(RColorBrewer)
library(Cairo)
library(beeswarm)

# Define all spine traits to loop over
allTraits <- c(
  "total_length_um", "total_spines", "spine_density", "spine_length",
  "volume", "head_diameter",
  "thin_spine_density", "stubby_spine_density", "mushroom_spine_density", "filopodia_spine_density",
  "thin_spine_length", "stubby_spine_length", "mushroom_spine_length", "filopodia_spine_length",
  "thin_head_diameter", "stubby_head_diameter", "mushroom_head_diameter", "filopodia_head_diameter",
  "thin_volume", "stubby_volume", "mushroom_volume", "filopodia_volume"
)

# Sanity check
setdiff(allTraits, colnames(numericMeta1))  # should be empty

# Open PDF device
CairoPDF(file = paste0(rootdir, "01.30.26_spine_plots_fixed.pdf"), width = 16, height = 12)

# Set up layout and margins
par(mfrow = c(4, 6))
par(mar = c(4.5, 6, 4.5, 1.5))

# Loop through each module
for (i in seq_len(nrow(toplot))) {
  
  # Extract module vector
  y <- toplot[i, ]
  
  # Extract module color (column of MEs)
  moduleCol <- colnames(MEs)[i]
  
  # Extract module name (just the color, remove any "M# " prefix or "p=" suffix)
  moduleName <- rownames(toplot)[i]
  moduleColor <- sub(".* ", "", moduleName)       # remove M# prefix
  moduleColor <- sub("p=.*$", "", moduleColor)    # remove p-value suffix
  
  # Panel 1: Boxplot + Beeswarm by group
  boxplot(y ~ factor(Grouping, c("Control", "AsymAD", "AD")),
          col = moduleCol, ylab = "Eigengene Value",
          main = moduleName, xlab = NULL, las = 2)
  beeswarm(y ~ factor(Grouping, c("Control", "AsymAD", "AD")),
           add = TRUE, method = "swarm", pch = 21,
           bg = paste0(col2hex(moduleCol), "99"), corral = "gutter")
  
  # Panels 2–N: Trait correlations
  for (trait in allTraits) {
    
    # Skip if trait doesn't exist
    if (!trait %in% colnames(numericMeta1)) next
    
    # Get numeric trait vector
    x_raw <- numericMeta1[[trait]]
    x <- suppressWarnings(as.numeric(as.character(x_raw)))
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 2) next
    
    # Extract correlation and p-value
    bicor_val <- moduleTraitCor[moduleColor, trait]
    pval_val  <- moduleTraitPvalue[moduleColor, trait]
    
    # Plot
    verboseScatterplot(
      x = x[ok],
      y = y[ok],
      xlab = trait,
      ylab = "Eigenprotein",
      abline = TRUE,
      cex.axis = 1, cex.lab = 1, cex = 1,
      col = "black", bg = moduleCol, pch = 21,
      main = paste0(
        "bicor=", if (is.numeric(bicor_val) && is.finite(bicor_val)) signif(bicor_val, 2) else "NA",
        ", p=", if (is.numeric(pval_val) && is.finite(pval_val)) signif(pval_val, 2) else "NA"
      ),
      col.main = if (is.numeric(pval_val) && pval_val < 0.05) "red" else "black"
    )
  }
  
  frame()  # start new page if layout full
}

# Close the PDF device
dev.off()




## EH tells us the module’s overall expression is associated with disease status
# Sort all modules by ascending p-value
sortedPV <- sort(pvec)

# View the five smallest
head(sortedPV, 42)

#d2    5.3.2025
# brown    royalblue    lightcyan      magenta midnightblue 
# 0.0001046431 0.0003628619 0.0047374295 0.0069079084 0.0094156482 

#d4    5.3.2025
# red    darkgreen     darkgrey  greenyellow   lightcyan1 
# 5.039833e-05 2.708332e-04 3.420530e-04 1.174826e-03 1.217107e-03 


#d4
# black         cyan     darkgrey       yellow   orangered4 
# 0.0001534039 0.0001678707 0.0010807062 0.0022883984 0.0026763605 


########################################
#Write Module Membership/kME table
orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:nrow(cleanDatReg)),rownames(cleanDatReg),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
write.table(kMEtable,file=paste0(outputtabs,"/01.30.26.power8.ModuleAssignments-",FileBaseName,".txt"),sep="\t",row.names=FALSE)
#(load above file in excel and apply green-yellow-red conditional formatting heatmap to the columns with kME values); then save as excel.

## saved image of R session
save.image(paste0("01.30.26.module_membership_saved.image.",projectFilesOutputTag,".Rdata"))  #overwrites


# GO ELITE #
######################## EDIT THESE VARIABLES (USER PARAMETERS SET IN GLOBAL ENVIRONMENT) ############################################
#inputFile <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
#inputFile <- "ModuleAssignments_Jingting32TH_BOOTaspRegr_power8_MergeHeight0.07_PAMstageTRUE_ds2.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
#INPUT CSV FILE - in the filePath folder.
#Can be formatted as Kme table from WGCNA pipeline, or
#can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
#in this case, the longest list is used as background or the "universe" for the FET contingencies
#  For simple columnwise list input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.

filePath <- "C:/Users/ehobby/Documents/EH_Emory41_Update/GO_Emory41"   #gsub("//","/",outputfigs)
#Folder that (may) contain the input file specified above, and which will contain the outFilename project Folder.

outFilename <- "GO_Emory41_Results_d4"
#SUBFOLDER WITH THIS NAME WILL BE CREATED, and .PDF + .csv file using the same name will be created within this folder.

outputGOeliteInputs=FALSE
#If TRUE, GO Elite background file and module or list-specific input files will be created in the outFilename subfolder.

maxBarsPerOntology=5
#Ontologies per ontology type, used for generating the PDF report; does not limit tabled output

GMTdatabaseFile="C:/Users/ehobby/Documents/EH_Emory41_Update/GO_Emory41/Human_GO_AllPathways_noPFOCR_with_GO_iea_September_16_2024_symbol.gmt"   # e.g. "Human_GO_AllPathways_with_GO_iea_June_01_2022_symbol.gmt"
# Current month release will be downloaded if file does not exist.
# **Specify a nonexistent file to always download the current database to this folder.**
# Database .GMT file will be saved to the specified folder with its original date-specific name.
#path/to/filename of ontology database for the appropriate species (no conversion is performed)
#BaderLab website links to their current monthly update of ontologies for Human, Mouse, and Rat, minimally
#http://download.baderlab.org/EM_Genesets/current_release/
#For more information, see documentation:  http://baderlab.org/GeneSets

panelDimensions=c(3,2)    #dimensions of the individual parblots within a page of the main barplot PDF output
pageDimensions=c(8.5,11)  #main barplot PDF output page dimensions, in inches

color=c("darkseagreen3","lightsteelblue1","lightpink4","goldenrod","darkorange","gold")

# color <- c("#D0E3CA", "#A4D4A0", "#76C37A", "#4FA554", "#347C3A", "#1B4E23")

#colors respectively for ontology Types:
#"Biological Process","Molecular Function","Cellular Component","Reactome","WikiPathways","MSig.C2"
#must be valid R colors

modulesInMemory=TRUE
#uses cleanDat, net, and kMEdat from WGCNA systems biology pipeline already run, and these variables must be in memory
#inputFile will be ignored
ANOVAgroups=FALSE
#if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
#inputFile will be ignored

############ MUST HAVE AT LEAST 2 THREADS ENABLED TO RUN ############################################################################

parallelThreads=20

removeRedundantGOterms=TRUE
#if true, the 3 GO ontology types are collapased into a minimal set of less redundant terms using the below OBO file
GO.OBOfile<-"C:/Users/ehobby/Documents/EH_Emory41_Update/GO_Emory41/go.obo"
#only used and needed if above flag to remove redundant GO terms is TRUE.
#Download from http://current.geneontology.org/ontology/go.obo will commence into the specified folder if the specified filename does not exist.
#Does not appear to be species specific, stores all GO term relations and is periodically updated.

cocluster=TRUE
#If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)

######################## END OF PARAMETER VARIABLES ###################################################################################

# colnames(MEs)[19] <- "lightyellow"
# colnames(MEs)[27] <- "white"

library(piano)
source("C:/Users/ehobby/Documents/EH_Emory41_Update/GO_Emory41/GOparallel-FET.R")
GOparallel()  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.
# priority is given to modulesInMemory


## saved image of R session
save.image(paste0("d4.power8.final.saved.image.",projectFilesOutputTag,".Rdata"))  #overwrites


#############################################

## speakeasy 2 code

#load Rdata file from WGCNA aka d4.power8.final.saved.image.EH_41BULK

library(speakeasyR)

library(WGCNA) # Network analysis package
library(NMF) # this package has a great annotated heatmap function - aheatmap
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(Cairo) # nicer graphics, anti-aliased, etc. --text from windows output PDFs using CairoPDF() function may not load in Illustrator, though -- so also use the pdf() standard output function when generating PDF figures
##Only for macs:
#CairoFonts(regular="Arial:style=Regular",bold="Arial:style=Bold",italic="Arial:style=Italic",bolditalic="Arial:style=Bold Italic,BoldItalic",symbol="Symbol")


adj_cleanDatReg <- adjacency(t(cleanDatReg), type="signed", power=1) #calculates signed adjacency

#Used one level of subclustering b/c gaiteri paper shows 1 level of subclustering could seperate "large communities into smaller communities"
#SE2 provides subclustering where the individual communities of the initial clustering will in turn be clustered into smaller communities. 
#This behavior can be turned on by setting the subclusters to parameter to a value greater than 1. 
#(The min_clust parameter determines the smallest community to consider for subclustering, if a community has fewer than min_clust nodes, it will not be subclustered further.)
set.seed(111)
se_mod <- speakeasyR::cluster(adj_cleanDatReg, seed = 111, subcluster = 2, min_clust = 100, verbose = TRUE, is_directed = TRUE) 
ordering <- speakeasyR::order_nodes(adj_cleanDatReg, se_mod)


# confirm levels
dim(se_mod)                 # should be 3 x N
length(unique(se_mod[1, ]))
length(unique(se_mod[2, ]))
length(unique(se_mod[3, ]))


level <- 2
level_order <- ordering[level,]
level_memb <- se_mod[level, level_order]
color <- labels2colors(level_memb)


#save clustering heatmap as PNG, PDF literally does not load
heatmap(adj_cleanDatReg[level_order, level_order], scale = "none", Rowv = NA, Colv = NA, 
        RowSideColors = color, xlab = "SE2 Module")


# labeling the se2 modules with the protein names from cleanDatReg

SE2modColor <- labels2colors(se_mod[level, ])
names(SE2modColor) <- rownames(cleanDatReg)

stopifnot(identical(colnames(adj_cleanDatReg), rownames(adj_cleanDatReg)))


table(SE2modColor)
# black         blue        brown         cyan        green  greenyellow      magenta midnightblue 
#     2         1254            3         1069          987            2          645         1034 
#  pink       purple          red       salmon          tan    turquoise       yellow 
#   657          565            4         1164            1          822            3 

# modules with binned grey module
# EH

# blue         cyan        green         grey      magenta    midnightblue         pink 
# 1254         1069          987           15          645         1034            657 
# purple       salmon    turquoise 
# 565         1164          822 



##assign color to each protein species similar to net$colors in WGCNA and run plots below 
##if <50 species in a cluster (module) then re-allocate to grey

allocate_grey_SE2 <- names(table(SE2modColor)[sapply(table(SE2modColor), FUN = function(x)x<50)])
SE2modColor[SE2modColor %in% allocate_grey_SE2] <- "grey"
orderedModules <-  matrix(nrow=length(unique(SE2modColor)), ncol=2)


colnames(orderedModules) <- c("Mnum", "Color")
orderedModules[,1] <- paste0("M", 1:length(unique(SE2modColor)))
orderedModules[,2] <- c("blue", "green", "turquoise", "purple", "magenta", "pink", "midnightblue", "cyan", "salmon", "grey")

orderedModules[match("grey", orderedModules[,2]),1] <- NA

numericMeta2 <- readRDS("numericMetawSpines.rds")
numericMeta <- numericMeta1

numericMeta <- numericMeta1_CW_spines_01.31.26

## 02.01.26
## spine data fix

numericMeta1 <- readRDS("numericMeta1_CW_spines_01.31.26.rds")
numericMeta <- numericMeta1

saveRDS(cleanDatReg, file = "cleanDatReg")



######################################################################################################
## Output GlobalNetworkPlots and kMEtable
####################################################################################################################


projectFilesOutputTag = "EH_41BULK"
FileBaseName=paste0(projectFilesOutputTag,"_SE2")


out_dir <- "C:/Users/ehobby/Documents/SE2_rerun_02.01.26"
out_file <- file.path(out_dir, "02.16.26_GlobalNetworkPlots-41Bulk.pdf")

CairoPDF(file = out_file, width = 16, height = 12)

#plot(1:10, 1:10, main = "TEST")

## Plot dendrogram with module colors and trait correlations
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDatReg), colors = SE2modColor)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)

#numericIndices<-unique(c( which(!is.na(apply(numericMeta,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta,2,function(x) sum(as.numeric(x),na.rm=T)))==0) ))
#Warnings OK; This determines which traits are numeric and if forced to numeric values, non-NA values do not sum to 0
#custom column selection in graphs
numericIndices<-c(8:10, 14:21, 25:46)

geneSignificance <- cor(sapply(numericMeta[,numericIndices],as.numeric),t(cleanDatReg),use="pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
geneSigColors <- t(numbers2colors(t(geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]


######################
## Find differences between Groups (as defined in Traits input file); Finalize Grouping of Samples for ANOVA

#Set a vector of strings that represent each sample in order, calling out each sample as a member of named groups (used by GlobalNetworkPlot boxplots, and later, ANOVA DiffEx)
Grouping<-numericMeta$Group  #typically there is a column "Group" loaded as a column in the traits.csv file

ApoE<- numericMeta$ApoE

ABC<- numericMeta$ABC
#ABC[numericMeta$ABC==0]<-"None"
#ABC[numericMeta$ABC==1]<-"Low"
#ABC[numericMeta$ABC==2]<-"Intermediate"
#ABC[numericMeta$ABC==3]<-"High"


# This gets ANOVA (Kruskal-Wallis) nonparametric p-values for groupwise comparison of interest.
# look at numericMeta (traits data) and choose traits to use for linear model-determination of p value
head(numericMeta)


# # Change below line to point to a factored trait, which will define groups for ANOVA
# regvars <- data.frame(as.factor( Grouping ), as.numeric(numericMeta$Age), as.numeric(numericMeta$Sex))
# colnames(regvars) <- c("Grouping","Age","Sex") ## data frame with covaraites incase we want to try multivariate regression
# ##aov1 <- aov(data.matrix(MEs)~AD,data=regvars) ## ANOVA framework yields same results
# lm1 <- lm(data.matrix(MEs)~Grouping,data=regvars) #sex and age effects are removed by the linear model
# 
# pvec <- rep(NA,ncol(MEs))
# for (i in 1:ncol(MEs)) {
#   f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
#   pvec[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
# }
# names(pvec) <- colnames(MEs)
# 
# ## Find differences between APOE Risk Groups
# regvars <- data.frame(as.factor(numericMeta[,"ApoE"]),as.numeric(numericMeta[,"Age"]),as.factor(numericMeta[,"Sex"]),as.numeric(numericMeta[,"PMI"]))
# colnames(regvars) <- c("Grouping","Age","batch","PMI") ## data frame with covaraites in case we want to try multivariate regression
# #aov1 <- aov(data.matrix(MEs)~Group,data=regvars) ## ANOVA framework yields same results
# lm1 <- lm(data.matrix(MEs)~Grouping,data=regvars)
# 
# pvec.Apoe <- rep(NA,ncol(MEs))
# for (i in 1:ncol(MEs)) {
#   f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
#   pvec.Apoe[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
# }
# names(pvec.Apoe) <- colnames(MEs)


## EH editts - changed covariates to Batch and PMI

# Create the covariate data frame
regvars <- data.frame(
  Grouping = as.factor(numericMeta$Group),
  Batch = as.factor(numericMeta$Batch),  # treat Batch as a factor
  PMI = as.numeric(numericMeta$PMI)      # treat PMI as numeric
)
colnames(regvars) <- c("Grouping","Batch","PMI") ## data frame with covaraites incase we want to try multivariate regression

# Fit the multivariate linear model (one model across all SE2 modules)
lm1 <- lm(data.matrix(MEs) ~ Grouping + Batch + PMI, data = regvars)

# Extract p-values for each module
pvec <- rep(NA, ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic
  pvec[i] <- pf(f[1], f[2], f[3], lower.tail = FALSE)
}
names(pvec) <- colnames(MEs)

# (Optional) Adjust for multiple comparisons
pvec_fdr <- p.adjust(pvec, method = "fdr")


## Find differences between APOE Risk Groups
regvars2 <- data.frame(
  ApoE = as.factor(numericMeta$ApoE),
  Batch = as.factor(numericMeta$Batch),  # treat Batch as a factor
  PMI = as.numeric(numericMeta$PMI)      # treat PMI as numeric
)
colnames(regvars2) <- c("Grouping","Batch","PMI") ## data frame with covaraites incase we want to try multivariate regression

# Fit the multivariate linear model (one model across all SE2 modules)
lm1 <- lm(data.matrix(MEs) ~ Grouping + Batch + PMI, data = regvars2)

pvec.Apoe <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec.Apoe[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec.Apoe) <- colnames(MEs)



######################
## Get sigend kME values
kMEdat <- signedKME(t(cleanDatReg), MEList$eigengenes, corFnc="bicor")


######################
## Plot eigengene-trait correlations - using p value of bicor for heatmap scale
library(RColorBrewer)
MEcors <- bicorAndPvalue(MEs,numericMeta[,numericIndices])
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p


textMatrix = apply(moduleTraitCor,2,function(x) signif(x, 2))
#textMatrix = paste(signif(moduleTraitCor, 2), " (",
#  signif(moduleTraitPvalue, 1), ")", sep = "");
#dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3));

## Display the correlation values within a heatmap plot
cexy <- if(nModules>75) { 0.8 } else { 1 }
colvec <- rep("white",1500)
colvec[1:500] <- colorRampPalette(rev(brewer.pal(8,"BuPu")[2:8]))(500)
colvec[501:1000]<-colorRampPalette(c("white",brewer.pal(8,"BuPu")[2]))(3)[2] #interpolated color for 0.05-0.1 p
labeledHeatmap(Matrix = apply(moduleTraitPvalue,2,as.numeric),
               xLabels = colnames(numericMeta)[numericIndices],
               yLabels = paste0("ME",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colvec,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y= cexy,
               zlim = c(0,0.15),
               main = paste("Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Student correlation p value"),
               cex.main=0.8)


######################
## Plot eigengene-trait heatmap custom - using bicor color scale

numericMetaCustom<-numericMeta[,numericIndices]
MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

moduleTraitPvalue<-signif(moduleTraitPvalue, 1)
moduleTraitPvalue[moduleTraitPvalue > as.numeric(0.05)]<-as.character("")

textMatrix = moduleTraitPvalue; #paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#textMatrix = gsub("()", "", textMatrix,fixed=TRUE)

labelMat<-matrix(nrow=(length(names(MEs))), ncol=2,data=c(rep(1:(length(names(MEs)))),names(MEs)))
labelMat<-labelMat[match(names(MEs),labelMat[,2]),]
for (i in 1:(length(names(MEs)))) { labelMat[i,1]<-paste("M",labelMat[i,1],sep="") }
for (i in 1:length(names(MEs))) { labelMat[i,2]<-paste("ME",labelMat[i,2],sep="") }

#rowMin(moduleTraitPvalue) # if we want to resort rows by min P value in the row
xlabAngle <- if(nModules>75) { 90 } else { 45 }


par(mar=c(16, 12, 3, 3) )
par(mfrow=c(1,1))

bw<-colorRampPalette(c("#0058CC", "white"))
wr<-colorRampPalette(c("white", "#CC3300"))

colvec<-c(bw(50),wr(50))

labeledHeatmap(Matrix = t(moduleTraitCor)[,],
               yLabels = colnames(numericMetaCustom),
               xLabels = labelMat[,2],
               xSymbols = labelMat[,1],
               xColorLabels=TRUE,
               colors = colvec,
               textMatrix = t(textMatrix)[,],
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x = cexy,
               xLabelsAngle = xlabAngle,
               verticalSeparator.x=c(rep(c(1:length(colnames(MEs))),as.numeric(ncol(MEs)))),
               verticalSeparator.col = 1,
               verticalSeparator.lty = 1,
               verticalSeparator.lwd = 1,
               verticalSeparator.ext = 0,
               horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
               horizontalSeparator.col = 1,
               horizontalSeparator.lty = 1,
               horizontalSeparator.lwd = 1,
               horizontalSeparator.ext = 0,
               zlim = c(-1,1),
               main = "Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
               cex.main=0.8)


## Plot annotated heatmap - annotate all the metadata, plot the eigengenes!
# This is where we will first use the Grouping vector of string group descriptions we set above.
toplot <- MEs

colnames(toplot) <- colnames(MEs)
rownames(toplot) <- rownames(MEs)
toplot <- t(toplot)

pvec <- pvec[match(names(pvec),rownames(toplot))]
#rownames(toplot) <- paste(rownames(toplot),"\np = ",signif(pvec,2),sep="")
rownames(toplot) <- paste(orderedModules[match(colnames(MEs),orderedModules[,2]),1]," ",rownames(toplot),"  |  K-W p=",signif(pvec,2),sep="")

# add any traits of interest you want to be in the legend
Gender=as.numeric(numericMeta$Sex)
Gender[Gender==0]<-"Female"
Gender[Gender==1]<-"Male"
metdat=data.frame(Group=Grouping,Age=as.numeric(numericMeta$Age), Gender=Gender)

# set colors for the traits in the legend
heatmapLegendColors=list('Group'=c("dodgerblue","goldenrod","seagreen3","hotpink","purple"),
                         'Age'=c("white","darkgreen"), #young to old
                         'Gender'=c("pink","dodgerblue"), #F, M
                         'Modules'=sort(colnames(MEs)))

library(NMF)
par(mfrow=c(1,1))
aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES IN ORIGINAL, e.g. BATCH OR REGION ORDER",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=NA) ## Do not cluster columns - keep given order

aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES CLUSTERED",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE,Colv=TRUE) ## Cluster columns





######################################
## Change the below code in the for loop using the following session output
library(gplots) #for col2hex() fn
library(beeswarm)

## Get module-trait bicor correlations (append to verboseScatterplot title below)
numericMetaCustom<-numericMeta[,numericIndices]
MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p

#These are your numerically coded traits:
colnames(numericMeta)[numericIndices] #choose traits for correlation scatterplots (verboseScatterplot functions below)

#These are your ANOVA sample groups and the number of samples in each
table(Grouping) #alphabetically ordered, you choose the order of groups in the boxplot function by typing them in

## Make changes after checking output on console for the above 2 lines
par(mfrow=c(4,6))
par(mar=c(4.5,6,4.5,1.5))

for (i in 1:(nrow(toplot))) {  # grey already excluded, no -1
  titlecolor<-if(signif(pvec,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(Grouping,c("CT","AsymAD","AD")),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nK-W p = ",signif(pvec,2)[i]),xlab=NULL,las=2,col.main=titlecolor)  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(Grouping,c("CT","AsymAD","AD")),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism ; #bg=goldenrod #DDA43B "#DDA43B99"
  
  titlecolor<-if(signif(pvec.Apoe,2)[i] <0.05) { "red" } else { "black" }
  boxplot(toplot[i,]~factor(ApoE,c("E2/3","E3/3","E3/4","E4/4")),col=colnames(MEs)[i],ylab="Eigenprotein Value",main=paste0(orderedModules[match(colnames(MEs)[i],orderedModules[,2]),1]," ",colnames(MEs)[i],"\nK-W p = ",signif(pvec.Apoe,2)[i]),xlab="APOE genotype",las=2,col.main=titlecolor)  #no outliers: ,outline=FALSE)
  transcol=paste0(col2hex(colnames(MEs)[i]),"99")
  beeswarm(toplot[i,]~factor(ApoE,c("E2/3","E3/3","E3/4","E4/4")),method="swarm",add=TRUE,corralWidth=0.5,vertical=TRUE,pch=21,bg=transcol,col="black",cex=0.8,corral="gutter") #more like prism ; #bg=goldenrod #DDA43B "#DDA43B99"
  
  verboseScatterplot(x=numericMeta[,"CERAD"],y=toplot[i,],xlab="CERAD Score",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"CERAD"],2),", p=",signif(moduleTraitPvalue[i,"CERAD"],2),"\n"),col.main=if(moduleTraitPvalue[i,"CERAD"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta[,"BRAAK"],y=toplot[i,],xlab="Braak Score",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"BRAAK"],2),", p=",signif(moduleTraitPvalue[i,"BRAAK"],2),"\n"),col.main=if(moduleTraitPvalue[i,"BRAAK"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta[,"FrontalNP"],y=toplot[i,],xlab="FrontalNP",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FrontalNP"],2),", p=",signif(moduleTraitPvalue[i,"FrontalNP"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FrontalNP"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta[,"FrontalDP"],y=toplot[i,],xlab="FrontalDP",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FrontalDP"],2),", p=",signif(moduleTraitPvalue[i,"FrontalDP"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FrontalDP"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta[,"FrontalNFT"],y=toplot[i,],xlab="FrontalNFT",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"FrontalNFT"],2),", p=",signif(moduleTraitPvalue[i,"FrontalNFT"],2),"\n"),col.main=if(moduleTraitPvalue[i,"FrontalNFT"]<0.05) { "red" } else { "black" })
  
  verboseScatterplot(x=numericMeta[,"MMSE"],y=toplot[i,],xlab="MMSE",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"MMSE"],2),", p=",signif(moduleTraitPvalue[i,"MMSE"],2),"\n"),col.main=if(moduleTraitPvalue[i,"MMSE"]<0.05) { "red" } else { "black" })
  
}

##outputs sample-by-sample eigenprotein barplots (not useful for large number of samples)
#while(!par('page')) plot.new()
#for (i in 1:nrow(toplot)) {
# barplot(height=rev(toplot[i,]),width=5,col=colnames(MEs)[i],xlab=paste(colnames(MEs)[i]," Eigenprotein Relative Expression"),main=rownames(toplot)[i],ylab=NULL,las=2,space=0.4,horiz=TRUE) #las=2 for rotated 90° X-axis labels  main=rownames(toplot)[i]
## text(bargr,par("usr")[3] - 0.025, srt=45, adj =1, labels= c(colnames(toplot)),xpd=TRUE,font=2) # bargr <- barplot(... above; gives rotated 45° x-axis labels but overwrites on top of existing ones
#}

dev.off()

# #### All spine corr plots
# 
# CairoPDF(file=paste0(rootdir,"41Bulk-module-spine-plots.pdf"),width=16,height=12)
# 
# 
## Make changes after checking output on console for the above 2 lines
par(mfrow=c(4,6))
par(mar=c(4.5,6,4.5,1.5))

for (i in 1:(nrow(toplot))) {  # grey already excluded, no -1
  verboseScatterplot(x=numericMeta$Spine.Density.10?m,y=toplot[i,],xlab="Spine Density 10um",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Spine.Density.10?m"],2),", p=",signif(moduleTraitPvalue[i,"Spine.Density.10?m"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Spine.Density.10?m"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Spine.Length..?m.,y=toplot[i,],xlab="Spine.Length..?m.",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Spine.Length..?m."],2),", p=",signif(moduleTraitPvalue[i,"Spine.Length..?m."],2),"\n"),col.main=if(moduleTraitPvalue[i,"Spine.Length..?m."]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Head.Diameter.?m.,y=toplot[i,],xlab="Head.Diameter.?m.",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Head.Diameter.?m."],2),", p=",signif(moduleTraitPvalue[i,"Head.Diameter.?m."],2),"\n"),col.main=if(moduleTraitPvalue[i,"Head.Diameter.?m."]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Neck.Diameter,y=toplot[i,],xlab="Neck.Diameter",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Neck.Diameter"],2),", p=",signif(moduleTraitPvalue[i,"Neck.Diameter"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Neck.Diameter"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Thin.SD,y=toplot[i,],xlab="Thin.SD",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Thin.SD"],2),", p=",signif(moduleTraitPvalue[i,"Thin.SD"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Thin.SD"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Stubby.SD,y=toplot[i,],xlab="Stubby.SD",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Stubby.SD"],2),", p=",signif(moduleTraitPvalue[i,"Stubby.SD"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Stubby.SD"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Mushroom.SD,y=toplot[i,],xlab="Mushroom.SD",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Mushroom.SD"],2),", p=",signif(moduleTraitPvalue[i,"Mushroom.SD"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Mushroom.SD"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Filopodia.D,y=toplot[i,],xlab="Filopodia.D",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Filopodia.D"],2),", p=",signif(moduleTraitPvalue[i,"Filopodia.D"],2),"\n"),col.main=if(moduleTraitPvalue[i,"Filopodia.D"]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Thin.spine..,y=toplot[i,],xlab="Thin.spine..",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Thin.spine.."],2),", p=",signif(moduleTraitPvalue[i,"Thin.spine.."],2),"\n"),col.main=if(moduleTraitPvalue[i,"Thin.spine.."]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Stubby.spine..,y=toplot[i,],xlab="Stubby.spine..",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Stubby.spine.."],2),", p=",signif(moduleTraitPvalue[i,"Stubby.spine.."],2),"\n"),col.main=if(moduleTraitPvalue[i,"Stubby.spine.."]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Mushroom.spine..,y=toplot[i,],xlab="Mushroom.spine..",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Mushroom.spine.."],2),", p=",signif(moduleTraitPvalue[i,"Mushroom.spine.."],2),"\n"),col.main=if(moduleTraitPvalue[i,"Mushroom.spine.."]<0.05) { "red" } else { "black" })
  verboseScatterplot(x=numericMeta$Filopodia.spine..,y=toplot[i,],xlab="Filopodia.spine..",ylab="Eigenprotein",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col="black",bg=colnames(MEs)[i],pch=21,main=paste0("bicor=",signif(moduleTraitCor[i,"Filopodia.spine.."],2),", p=",signif(moduleTraitPvalue[i,"Filopodia.spine.."],2),"\n"),col.main=if(moduleTraitPvalue[i,"Filopodia.spine.."]<0.05) { "red" } else { "black" })
}

dev.off()

########################################
#Write Module Membership/kME table
orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(SE2modColor,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:nrow(cleanDatReg)),rownames(cleanDatReg),SE2modColor,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
write.table(kMEtable,file=paste0(rootdir,"/02.17.26_ModuleAssignments-41Bulk_SE2.txt"),sep="\t",row.names=FALSE)
#(load above file in excel and apply green-yellow-red conditional formatting heatmap to the columns with kME values); then save as excel.

## saved image of R session
save.image(paste0("02.01.26_saved.image.",projectFilesOutputTag,".Rdata"))  #overwrites





#######WRAPPER CALL GO-ELITE GENE ONTOLOGIES
#######***NOTE FROM EVAN LIU - I HAD TO EDIT LINE 350 TO MAKE THE LIST COLOR FROM SPEAKEASY2 AND NOT WGCNA

######################## EDIT THESE VARIABLES (USER PARAMETERS SET IN GLOBAL ENVIRONMENT) ############################################

colnames(kMEtable)[3] <- "net.colors"
write.csv(kMEtable, file=paste0(rootdir,"/ModuleAssignmentsForGOElite.csv"),row.names=FALSE)
NETcolors= SE2modColor    #module color assignments, vector of length equal to number of rows in cleanDat; should have all colors for modules from 1:minimumSizeRank as printed by WGCNA::labels2colors(1:nModules)
nModules <- 9

#inputFile <- "ENDO_MG_TWO_WAY_LIST_NTS_v02b_forGOelite.csv"                                            #Sample File 1 - has full human background
inputFile <- "ModuleAssignmentsForGOElite.csv"  #Sample File 2 - WGCNA kME table for (Dai, et al, 2019)
#INPUT CSV FILE - in the filePath folder.
#Can be formatted as Kme table from WGCNA pipeline, or
#can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
#in this case, the longest list is used as background or the "universe" for the FET contingencies
#  For simple columnwise list input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.

# filePath <- "Z:/Evan/Emory 41 Proteomics/Emory41/GOElite"   #gsub("//","/",outputfigs)
# #Folder that (may) contain the input file specified above, and which will contain the outFilename project Folder.

outFilename <- "SE2-GO"
#SUBFOLDER WITH THIS NAME WILL BE CREATED, and .PDF + .csv file using the same name will be created within this folder.

outputGOeliteInputs=FALSE
#If TRUE, GO Elite background file and module or list-specific input files will be created in the outFilename subfolder.

maxBarsPerOntology=5
#Ontologies per ontology type, used for generating the PDF report; does not limit tabled output

#GMTdatabaseFile="Z:/Evan/ROSMAP_Synaptosome_Proteome_2024/GOElite/Human_GO_AllPathways_noPFOCR_with_GO_iea_September_16_2024_symbol.gmt"   # e.g. "Human_GO_AllPathways_with_GO_iea_June_01_2022_symbol.gmt"
# Current month release will be downloaded if file does not exist.
# **Specify a nonexistent file to always download the current database to this folder.**
# Database .GMT file will be saved to the specified folder with its original date-specific name.
#path/to/filename of ontology database for the appropriate species (no conversion is performed)
#BaderLab website links to their current monthly update of ontologies for Human, Mouse, and Rat, minimally
#http://download.baderlab.org/EM_Genesets/current_release/
#For more information, see documentation:  http://baderlab.org/GeneSets

panelDimensions=c(3,2)    #dimensions of the individual parblots within a page of the main barplot PDF output
pageDimensions=c(8.5,11)  #main barplot PDF output page dimensions, in inches

color=c("darkseagreen3","lightsteelblue1","lightpink4","goldenrod","darkorange","gold")

#scale_fill_manual(values = c("#3FB8AF", "#7FC7AF", "#DAD8A7","#FF9E9D", "#FF3D7F"))

#color=c("#009392FF","#F1EAC8FF","#D0587EFF")
#color=c("#009392FF","#F1EAC8FF","#E5B9ADFF","#D98994FF","#D0587EFF")


#colors respectively for ontology Types:
#"Biological Process","Molecular Function","Cellular Component","Reactome","WikiPathways","MSig.C2"
#must be valid R colors

modulesInMemory=TRUE
#uses cleanDat, net, and kMEdat from WGCNA systems biology pipeline already run, and these variables must be in memory
#inputFile will be ignored
ANOVAgroups=FALSE
#if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
#inputFile will be ignored

############ MUST HAVE AT LEAST 2 THREADS ENABLED TO RUN ############################################################################

parallelThreads=20

removeRedundantGOterms=TRUE
#if true, the 3 GO ontology types are collapased into a minimal set of less redundant terms using the below OBO file
GO.OBOfile<-"C:/Users/ehobby/Documents/SE2_EH/GOelite/go.obo"
#only used and needed if above flag to remove redundant GO terms is TRUE.
#Download from http://current.geneontology.org/ontology/go.obo will commence into the specified folder if the specified filename does not exist.
#Does not appear to be species specific, stores all GO term relations and is periodically updated.

cocluster=TRUE
#If TRUE, output PDF of signed Zscore coclustering on GO:cellular component terms (useful for WGCNA modules)

######################## END OF PARAMETER VARIABLES ###################################################################################


library(piano)
source("GOparallel-FET.R")
GOparallel()  # parameters are set in global environment as above; if not set, the function falls back to defaults and looks for all inputs available.
# priority is given to modulesInMemory































