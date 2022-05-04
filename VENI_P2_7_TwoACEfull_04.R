# ----------------------------------------------------------------------------------------------------------------------
# Program: oneACEvc.R  
#  Author: Hermine Maes (2018-10-22)
# Adapted: Antonie Knigge  
#    Date: 2020-01-16 
#  Update: 2022-03-01
# Twin bivariate ACE model to estimate causes of variation with moderation
# Matrix style model - Raw data - Continuous data
# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|

# Load Libraries & Options
rm(list=ls())
library(psych) 
library(polycor)
library(foreign)
library(Hmisc)
library(readstata13)
library(OpenMx)
library(plyr)
library(dplyr)
library(plotrix)
source("miFunctions.R")
source("colorado-Rfunctions.r")

# Create Output 
filename    <- "VENI_P2_7_TwoACEfull_04"
sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

# Load Data
#setwd("C:/Users/knigg101/surfdrive/Onderzoek/VENI/P2/WorkInProgress/Analyses")
#setwd("/Users/Antonie/surfdrive/Onderzoek/VENI/P2/WorkInProgress/Analyses")

NtrRawData <- read.dta13("NTR_P2_1_06.dta", convert.dates=T, convert.underscore=F, convert.factors=F)
twinData   <- reshape(NtrRawData, idvar="twin_id", timevar="mult_ext2", 
                      v.names=c("Extension", "mult_ext", "lfts14", "lfts16", "lfts18", "age", "c_age", "PID", "male",  
                                "c_cito", "z_cito", "cito_final", "edu12", "edu14", "edu16", "edu18", "edu","edu_alt", "delay",
                                "achie4", "achie4_1", "achie4_2", "achie4_3", "achie4_4",
                                "achie5", "achie5_1", "achie5_2", "achie5_3", "achie5_4", "achie5_5", 
                                "achie6", "achie6_1", "achie6_2", "achie6_3", "achie6_4", "achie6_5", "achie6_6",
                                "achie8", "achie8_1", "achie8_2", "achie8_3", "achie8_4", "achie8_5", "achie8_6", "achie8_7", "achie8_8",
                                "twinsize1b", "twinsize2b", "defm", "cito_nm"), direction="wide", sep = "_")

dim(twinData)
names(twinData)
summary(twinData)
describe(twinData[,1:29], skew=F)

# Give pseudo-missing values for definition variables where only one twin was observed.
# (phenotypes already have missing values for those twin pairs; 
# pseudo-missings were already assigned in stata to twin pairs where both were observed but one has missing info)
twinData$male_2[is.na(twinData$male_2)] <- -999
twinData$defm_2[is.na(twinData$defm_2)] <- -999

# Select Variables for Analysis
vars            <- c('c_cito','edu')         # list of variables names
nv              <- length(vars)	             # number of variables
ntv             <- nv*2                      # number of total variables
selVars         <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="_")
covVars1        <- c('male_1', 'c_yob', 'male_2')
nc1             <- 2                         # number of covariates
#twinData$defm_1 <- twinData$c_cito_1         # create a copy of moderator as a definition variable
#twinData$defm_2 <- twinData$c_cito_2
defVars         <- c('defm_1', 'defm_2')
useVars	        <- c(selVars, covVars1, defVars, 'mz', 'delay_1', 'delay_2', 'cito_final_1', 'cito_final_2')

# Select Cases for Analysis
twinDataNM <- subset(twinData, (male_1!="-999" | male_2!="-999") & c_yob!="NA" & (c_cito_1!="NA" | c_cito_2!="NA" | edu_1!="NA" | edu_2!="NA"), c(useVars, 'FamilyNumber'))
dim(twinDataNM)
summary(twinDataNM)

# Select Data for Analysis
mzData	<- subset(twinDataNM, mz==1, c(selVars, defVars, covVars1))
dzData	<- subset(twinDataNM, mz==0, c(selVars, defVars, covVars1))
dim(mzData)
dim(dzData)

# Create groups based on tracking status (immediate, delayed, missing info)
mz00Data <- subset(       twinDataNM, mz==1 & (delay_1==0     & delay_2==0)     , c(selVars, defVars, covVars1))
mz11Data <- subset(       twinDataNM, mz==1 & (delay_1==1     & delay_2==1)     , c(selVars, defVars, covVars1))
mz01Data <- subset(       twinDataNM, mz==1 & (delay_1==0     & delay_2==1)     , c(selVars, defVars, covVars1))
mz0mData <- select(filter(twinDataNM, mz==1 & (delay_1==0     & is.na(delay_2))), c(selVars, defVars, covVars1))
mz1mData <- select(filter(twinDataNM, mz==1 & (delay_1==1     & is.na(delay_2))), c(selVars, defVars, covVars1))
mzmmData <- select(filter(twinDataNM, mz==1 & (is.na(delay_1) & is.na(delay_2))), c(selVars, defVars, covVars1))
dz00Data <- subset(       twinDataNM, mz==0 & (delay_1==0     & delay_2==0)     , c(selVars, defVars, covVars1))
dz11Data <- subset(       twinDataNM, mz==0 & (delay_1==1     & delay_2==1)     , c(selVars, defVars, covVars1))
dz01Data <- subset(       twinDataNM, mz==0 & (delay_1==0     & delay_2==1)     , c(selVars, defVars, covVars1))
dz0mData <- select(filter(twinDataNM, mz==0 & (delay_1==0     & is.na(delay_2))), c(selVars, defVars, covVars1))
dz1mData <- select(filter(twinDataNM, mz==0 & (delay_1==1     & is.na(delay_2))), c(selVars, defVars, covVars1))
dzmmData <- select(filter(twinDataNM, mz==0 & (is.na(delay_1) & is.na(delay_2))), c(selVars, defVars, covVars1))
dim(mz00Data)
dim(mz11Data)
dim(mz01Data)
dim(mz0mData)
dim(mz1mData)
dim(mzmmData)
dim(dz00Data)
dim(dz11Data)
dim(dz01Data)
dim(dz0mData)
dim(dz1mData)
dim(dzmmData)

# Generate Descriptive Statistics
colMeans(mz00Data,na.rm=TRUE)
colMeans(mz11Data,na.rm=TRUE)
colMeans(mz01Data,na.rm=TRUE)
colMeans(mz0mData,na.rm=TRUE)
colMeans(mz1mData,na.rm=TRUE)
colMeans(mzmmData,na.rm=TRUE)
colMeans(dz00Data,na.rm=TRUE)
colMeans(dz11Data,na.rm=TRUE)
colMeans(dz01Data,na.rm=TRUE)
colMeans(dz0mData,na.rm=TRUE)
colMeans(dz1mData,na.rm=TRUE)
colMeans(dzmmData,na.rm=TRUE)
cov(mz00Data,use="complete")
cov(mz11Data,use="complete")
cov(mz01Data,use="complete")
cov(mz0mData,use="complete")
cov(mz1mData,use="complete")
cov(mzmmData,use="complete")
cov(dz00Data,use="complete")
cov(dz11Data,use="complete")
cov(dz01Data,use="complete")
cov(dz0mData,use="complete")
cov(dz1mData,use="complete")
cov(dzmmData,use="complete")


### Correlations 

## Cito (Performance scores)
# Overall
rMZc <- cor(mzData$c_cito_1, mzData$c_cito_2, use = "pairwise.complete.obs")
rDZc <- cor(dzData$c_cito_1, dzData$c_cito_2, use = "pairwise.complete.obs")

# By tracking type
rMZ00c <- cor(mz00Data$c_cito_1, mz00Data$c_cito_2, use = "pairwise.complete.obs")
rDZ00c <- cor(dz00Data$c_cito_1, dz00Data$c_cito_2, use = "pairwise.complete.obs")
rMZ11c <- cor(mz11Data$c_cito_1, mz11Data$c_cito_2, use = "pairwise.complete.obs")
rDZ11c <- cor(dz11Data$c_cito_1, dz11Data$c_cito_2, use = "pairwise.complete.obs")
rMZ01c <- cor(mz01Data$c_cito_1, mz01Data$c_cito_2, use = "pairwise.complete.obs")
rDZ01c <- cor(dz01Data$c_cito_1, dz01Data$c_cito_2, use = "pairwise.complete.obs")
rMZ0mc <- cor(mz0mData$c_cito_1, mz0mData$c_cito_2, use = "pairwise.complete.obs")
rDZ0mc <- cor(dz0mData$c_cito_1, dz0mData$c_cito_2, use = "pairwise.complete.obs")
rMZ1mc <- cor(mz1mData$c_cito_1, mz1mData$c_cito_2, use = "pairwise.complete.obs")
rDZ1mc <- cor(dz1mData$c_cito_1, dz1mData$c_cito_2, use = "pairwise.complete.obs")
rMZmmc <- cor(mzmmData$c_cito_1, mzmmData$c_cito_2, use = "pairwise.complete.obs")
rDZmmc <- cor(dzmmData$c_cito_1, dzmmData$c_cito_2, use = "pairwise.complete.obs")
cor_cito <- rbind(cbind(rMZc  , rDZc  ),
                  cbind(rMZ00c, rDZ00c), 
                  cbind(rMZ11c, rDZ11c), 
                  cbind(rMZ01c, rDZ01c),
                  cbind(rMZ0mc, rDZ0mc),
                  cbind(rMZ1mc, rDZ1mc),
                  cbind(rMZmmc, rDZmmc))

## Educational attainment
# Overall
rMZe   <- cor(mzData$edu_1, mzData$edu_2, use = "pairwise.complete.obs")
rDZe   <- cor(dzData$edu_1, dzData$edu_2, use = "pairwise.complete.obs")
# By tracking type
rMZ00e <- cor(mz00Data$edu_1, mz00Data$edu_2, use = "pairwise.complete.obs")
rDZ00e <- cor(dz00Data$edu_1, dz00Data$edu_2, use = "pairwise.complete.obs")
rMZ11e <- cor(mz11Data$edu_1, mz11Data$edu_2, use = "pairwise.complete.obs")
rDZ11e <- cor(dz11Data$edu_1, dz11Data$edu_2, use = "pairwise.complete.obs")
rMZ01e <- cor(mz01Data$edu_1, mz01Data$edu_2, use = "pairwise.complete.obs")
rDZ01e <- cor(dz01Data$edu_1, dz01Data$edu_2, use = "pairwise.complete.obs")
rMZ0me <- cor(mz0mData$edu_1, mz0mData$edu_2, use = "pairwise.complete.obs")
rDZ0me <- cor(dz0mData$edu_1, dz0mData$edu_2, use = "pairwise.complete.obs")
rMZ1me <- cor(mz1mData$edu_1, mz1mData$edu_2, use = "pairwise.complete.obs")
rDZ1me <- cor(dz1mData$edu_1, dz1mData$edu_2, use = "pairwise.complete.obs")
rMZmme <- cor(mzmmData$edu_1, mzmmData$edu_2, use = "pairwise.complete.obs")
rDZmme <- cor(dzmmData$edu_1, dzmmData$edu_2, use = "pairwise.complete.obs")
cor_edu  <- rbind(cbind(rMZe  , rDZe  ), 
                  cbind(rMZ00e, rDZ00e), 
                  cbind(rMZ11e, rDZ11e), 
                  cbind(rMZ01e, rDZ01e),
                  cbind(rMZ0me, rDZ0me),
                  cbind(rMZ1me, rDZ1me),
                  cbind(rMZmme, rDZmme))
cor_cito
cor_edu


# ----------------------------------------------------------------------------------------------------------------------
# ACE
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# LABELING & STARTING VALUES

## Labeling

# Main paths
aLabs 	   <- c("aM","aC","aU")
aLabs0	   <- c("aM0","aC0","aU0")
aLabs1	   <- c("aM1","aC1","aU1")
aLabsm     <- c("aMm","aCm","aUm")

cLabs	     <- c("cM","cC","cU")
cLabs0	   <- c("cM0","cC0","cU0")
cLabs1	   <- c("cM1","cC1","cU1")
cLabsm	   <- c("cMm","cCm","cUm")

eLabs 	   <- c("eM","eC","eU")
eLabs0	   <- c("eM0","eC0","eU0")
eLabs1	   <- c("eM1","eC1","eU1")
eLabsm	   <- c("eMm","eCm","eUm")

# Moderating paths
aModLabs   <- c("aMod11","aMod21","aMod22")
aModLabs0  <- c("aMod11_0","aMod21_0","aMod22_0")
aModLabs1  <- c("aMod11_1","aMod21_1","aMod22_1")
aModLabsm  <- c("aMod11_m","aMod21_m","aMod22_m")

cModLabs   <- c("cMod11","cMod21","cMod22")
cModLabs0  <- c("cMod11_0","cMod21_0","cMod22_0")
cModLabs1  <- c("cMod11_1","cMod21_1","cMod22_1")
cModLabsm  <- c("cMod11_m","cMod21_m","cMod22_m")

eModLabs   <- c("eMod11","eMod21","eMod22")
eModLabs0  <- c("eMod11_0","eMod21_0","eMod22_0")
eModLabs1  <- c("eMod11_1","eMod21_1","eMod22_1")
eModLabsm  <- c("eMod11_m","eMod21_m","eMod22_m")

# Means
meanLabs 	 <- c("meanCITO","meanEDU")
meanLabs0  <- c("meanCITO0","meanEDU0")
meanLabs1  <- c("meanCITO1","meanEDU1")
meanLabsm  <- c("meanCITOm","meanEDUm")

## Set Starting Values
svM        <- c(0,2.7)                        # start value for means
svM0       <- c(0,2.6)                        # start value for means
svM1       <- c(0,2.9)                        # start value for means
svMm       <- c(0,2.6)                        # start value for means
svPa       <- c(7.0, 0.6, 0.5)                # start value for path coefficient for a
svPa0      <- c(7.0, 0.6, 0.5)                # start value for path coefficient for a
svPa1      <- c(5.4, 0.6, 0.6)                # start value for path coefficient for a
svPam      <- c(7.1, 0.6, 0.5)                # start value for path coefficient for a
svPc       <- c(2.6, 0.5, 0.0)                # start value for path coefficient for c
svPc0      <- c(3.5, 0.5, 0.1)                # start value for path coefficient for c
svPc1      <- c(1.1, 0.0, 0.0)                # start value for path coefficient for c
svPcm      <- c(2.6, 0.6, 0.1)                # start value for path coefficient for c
svPe       <- c(3.6, 0.2, 0.4)                # start value for path coefficient for e
svPe0      <- c(3.2, 0.2, 0.4)                # start value for path coefficient for e
svPe1      <- c(3.2, 0.2, 0.4)                # start value for path coefficient for e
svPem      <- c(3.9, 0.2, 0.5)                # start value for path coefficient for e
#svPaD     <- vech(diag(svPa,nv,nv))          # start values for diagonal of covariance matrix
#svPcD     <- vech(diag(svPc,nv,nv))          # start values for diagonal of covariance matrix
#svPeD     <- vech(diag(svPe,nv,nv))          # start values for diagonal of covariance matrix
lbPa      <- NA                               # start value for lower bounds
lbPaD     <- diag(lbPa,nv,nv)                 # lower bounds for diagonal of covariance matrix
lbPaD[lower.tri(lbPaD)] <- NA                 # lower bounds for below diagonal elements
lbPaD[upper.tri(lbPaD)] <- NA                 # lower bounds for above diagonal elements
pathModVal = c(0,0,0)                         # start values for moderation coefficients
B_SexVal = 0.5                                # start value for covariate sex
B_YobVal = 0.5                                # start value for covariate


# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create matrices to store observed values for each covariate
Male <- mxMatrix(type="Full", nrow=1, ncol=2, free=FALSE,
                 labels=c("data.male_1","data.male_2"), name="Male") 
Yob  <- mxMatrix(type="Full", nrow=1, ncol=2, free=FALSE,
                 labels=c("data.c_yob","data.c_yob"), name="Yob")

# Create matrices to store coefficients
bMale <- mxMatrix(type="Full", nrow=1, ncol=nv, free=TRUE,
                  values= .01, label=paste("betaMale",1:nv,sep=""), name="bMale")
bYob  <- mxMatrix(type="Full", nrow=1, ncol=nv, free=TRUE,
                  values= .01, label=paste("betaYob",1:nv,sep=""), name="bYob")

# Create Algebra for expected Mean Matrices
meanG0     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=svM0, labels=meanLabs0, name="meanG0" )
meanG1     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=svM1, labels=meanLabs1, name="meanG1" )
meanGm     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=svMm, labels=meanLabsm, name="meanGm" )

#  Add Covariates to the intercept
#  Note the use of the Kronecker-product operator, %x% :
meanG00 <- mxAlgebra(cbind(meanG0, meanG0) + (Male%x%bMale) + (Yob%x%bYob) , name="meanG00")       
meanG11 <- mxAlgebra(cbind(meanG1, meanG1) + (Male%x%bMale) + (Yob%x%bYob) , name="meanG11")       
meanG01 <- mxAlgebra(cbind(meanG0, meanG1) + (Male%x%bMale) + (Yob%x%bYob) , name="meanG01") 
meanG0m <- mxAlgebra(cbind(meanG0, meanGm) + (Male%x%bMale) + (Yob%x%bYob) , name="meanG0m")       
meanG1m <- mxAlgebra(cbind(meanG1, meanGm) + (Male%x%bMale) + (Yob%x%bYob) , name="meanG1m") 
meanGmm <- mxAlgebra(cbind(meanGm, meanGm) + (Male%x%bMale) + (Yob%x%bYob) , name="meanGmm") 

# Matrices a, c, and e to store a, c, and e Path Coefficients
pathA0    <- mxMatrix(name = "a0", type = "Lower", nrow = nv, ncol = nv, free=T, labels = aLabs0, values=svPa0, lbound=lbPaD)
pathA1    <- mxMatrix(name = "a1", type = "Lower", nrow = nv, ncol = nv, free=T, labels = aLabs1, values=svPa1, lbound=lbPaD)
pathAm    <- mxMatrix(name = "am", type = "Lower", nrow = nv, ncol = nv, free=T, labels = aLabsm, values=svPam, lbound=lbPaD)

pathC0    <- mxMatrix(name = "c0", type = "Lower", nrow = nv, ncol = nv, free=T, labels = cLabs0, values=svPc0, lbound=lbPaD)
pathC1    <- mxMatrix(name = "c1", type = "Lower", nrow = nv, ncol = nv, free=T, labels = cLabs1, values=svPc1, lbound=lbPaD)
pathCm    <- mxMatrix(name = "cm", type = "Lower", nrow = nv, ncol = nv, free=T, labels = cLabsm, values=svPcm, lbound=lbPaD)

pathE0    <- mxMatrix(name = "e0", type = "Lower", nrow = nv, ncol = nv, free=T, labels = eLabs0, values=svPe0, lbound=lbPaD)
pathE1    <- mxMatrix(name = "e1", type = "Lower", nrow = nv, ncol = nv, free=T, labels = eLabs1, values=svPe1, lbound=lbPaD)
pathEm    <- mxMatrix(name = "em", type = "Lower", nrow = nv, ncol = nv, free=T, labels = eLabsm, values=svPem, lbound=lbPaD)

# Matrices to store moderator path coefficients
modPathA0 <- mxMatrix(name="aMod0", "Lower", nrow=nv, ncol=nv, free=c(F,F,F), labels=aModLabs0, values=pathModVal )
modPathA1 <- mxMatrix(name="aMod1", "Lower", nrow=nv, ncol=nv, free=c(F,F,F), labels=aModLabs1, values=pathModVal )
modPathAm <- mxMatrix(name="aModm", "Lower", nrow=nv, ncol=nv, free=c(F,F,F), labels=aModLabsm, values=pathModVal )
modPathC0 <- mxMatrix(name="cMod0", "Lower", nrow=nv, ncol=nv, free=c(F,F,F), labels=cModLabs0, values=pathModVal )
modPathC1 <- mxMatrix(name="cMod1", "Lower", nrow=nv, ncol=nv, free=c(F,F,F), labels=cModLabs1, values=pathModVal )
modPathCm <- mxMatrix(name="cModm", "Lower", nrow=nv, ncol=nv, free=c(F,F,F), labels=cModLabsm, values=pathModVal )
modPathE0 <- mxMatrix(name="eMod0", "Lower", nrow=nv, ncol=nv, free=c(F,F,F), labels=eModLabs0, values=pathModVal )
modPathE1 <- mxMatrix(name="eMod1", "Lower", nrow=nv, ncol=nv, free=c(F,F,F), labels=eModLabs1, values=pathModVal )
modPathEm <- mxMatrix(name="eModm", "Lower", nrow=nv, ncol=nv, free=c(F,F,F), labels=eModLabsm, values=pathModVal )

# Matrices to store moderator itself (separate for each twin, because they can have different values on moderator)
mod_tw1  <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels='data.defm_1', name="Mod1")
mod_tw2  <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels='data.defm_2', name="Mod2")

# Matrices generated to hold A, C, and E computed Variance Components, for moderated models
varA1_00    <- mxAlgebra(name = "A1_00",  expression = (a0 + Mod1%x%aMod0) %*% t(a0+ Mod1%x%aMod0))
varC1_00    <- mxAlgebra(name = "C1_00",  expression = (c0 + Mod1%x%cMod0) %*% t(c0+ Mod1%x%cMod0))
varE1_00    <- mxAlgebra(name = "E1_00",  expression = (e0 + Mod1%x%eMod0) %*% t(e0+ Mod1%x%eMod0))

varA1_11    <- mxAlgebra(name = "A1_11",  expression = (a1 + Mod1%x%aMod1) %*% t(a1+ Mod1%x%aMod1))
varC1_11    <- mxAlgebra(name = "C1_11",  expression = (c1 + Mod1%x%cMod1) %*% t(c1+ Mod1%x%cMod1))
varE1_11    <- mxAlgebra(name = "E1_11",  expression = (e1 + Mod1%x%eMod1) %*% t(e1+ Mod1%x%eMod1))

varA1_01    <- mxAlgebra(name = "A1_01",  expression = (a0 + Mod1%x%aMod0) %*% t(a1+ Mod1%x%aMod1))
varC1_01    <- mxAlgebra(name = "C1_01",  expression = (c0 + Mod1%x%cMod0) %*% t(c1+ Mod1%x%cMod1))
varE1_01    <- mxAlgebra(name = "E1_01",  expression = (e0 + Mod1%x%eMod0) %*% t(e1+ Mod1%x%eMod1))

varA1_0m    <- mxAlgebra(name = "A1_0m",  expression = (a0 + Mod1%x%aMod0) %*% t(am+ Mod1%x%aModm))
varC1_0m    <- mxAlgebra(name = "C1_0m",  expression = (c0 + Mod1%x%cMod0) %*% t(cm+ Mod1%x%cModm))
varE1_0m    <- mxAlgebra(name = "E1_0m",  expression = (e0 + Mod1%x%eMod0) %*% t(em+ Mod1%x%eModm))

varA1_1m    <- mxAlgebra(name = "A1_1m",  expression = (a1 + Mod1%x%aMod1) %*% t(am+ Mod1%x%aModm))
varC1_1m    <- mxAlgebra(name = "C1_1m",  expression = (c1 + Mod1%x%cMod1) %*% t(cm+ Mod1%x%cModm))
varE1_1m    <- mxAlgebra(name = "E1_1m",  expression = (e1 + Mod1%x%eMod1) %*% t(em+ Mod1%x%eModm))

varA1_mm    <- mxAlgebra(name = "A1_mm",  expression = (am + Mod1%x%aModm) %*% t(am+ Mod1%x%aModm))
varC1_mm    <- mxAlgebra(name = "C1_mm",  expression = (cm + Mod1%x%cModm) %*% t(cm+ Mod1%x%cModm))
varE1_mm    <- mxAlgebra(name = "E1_mm",  expression = (em + Mod1%x%eModm) %*% t(em+ Mod1%x%eModm))


covA12_00   <- mxAlgebra(name = "A12_00", expression = (a0 + Mod1%x%aMod0) %*% t(a0+ Mod2%x%aMod0))
covC12_00   <- mxAlgebra(name = "C12_00", expression = (c0 + Mod1%x%cMod0) %*% t(c0+ Mod2%x%cMod0))
covE12_00   <- mxAlgebra(name = "E12_00", expression = (e0 + Mod1%x%eMod0) %*% t(e0+ Mod2%x%eMod0))

covA12_11   <- mxAlgebra(name = "A12_11", expression = (a1 + Mod1%x%aMod1) %*% t(a1+ Mod2%x%aMod1))
covC12_11   <- mxAlgebra(name = "C12_11", expression = (c1 + Mod1%x%cMod1) %*% t(c1+ Mod2%x%cMod1))
covE12_11   <- mxAlgebra(name = "E12_11", expression = (e1 + Mod1%x%eMod1) %*% t(e1+ Mod2%x%eMod1))

covA12_01   <- mxAlgebra(name = "A12_01", expression = (a0 + Mod1%x%aMod0) %*% t(a1+ Mod2%x%aMod1))
covC12_01   <- mxAlgebra(name = "C12_01", expression = (c0 + Mod1%x%cMod0) %*% t(c1+ Mod2%x%cMod1))
covE12_01   <- mxAlgebra(name = "E12_01", expression = (e0 + Mod1%x%eMod0) %*% t(e1+ Mod2%x%eMod1))

covA12_0m   <- mxAlgebra(name = "A12_0m", expression = (a0 + Mod1%x%aModm) %*% t(a0+ Mod2%x%aModm))
covC12_0m   <- mxAlgebra(name = "C12_0m", expression = (c0 + Mod1%x%cModm) %*% t(c0+ Mod2%x%cModm))
covE12_0m   <- mxAlgebra(name = "E12_0m", expression = (e0 + Mod1%x%eModm) %*% t(e0+ Mod2%x%eModm))

covA12_1m   <- mxAlgebra(name = "A12_1m", expression = (a1 + Mod1%x%aMod1) %*% t(am+ Mod2%x%aModm))
covC12_1m   <- mxAlgebra(name = "C12_1m", expression = (c1 + Mod1%x%cMod1) %*% t(cm+ Mod2%x%cModm))
covE12_1m   <- mxAlgebra(name = "E12_1m", expression = (e1 + Mod1%x%eMod1) %*% t(em+ Mod2%x%eModm))

covA12_mm   <- mxAlgebra(name = "A12_mm", expression = (am + Mod1%x%aModm) %*% t(am+ Mod2%x%aModm))
covC12_mm   <- mxAlgebra(name = "C12_mm", expression = (cm + Mod1%x%cModm) %*% t(cm+ Mod2%x%cModm))
covE12_mm   <- mxAlgebra(name = "E12_mm", expression = (em + Mod1%x%eModm) %*% t(em+ Mod2%x%eModm))


covA21_00   <- mxAlgebra(name = "A21_00", expression = (a0 + Mod2%x%aMod0) %*% t(a0+ Mod1%x%aMod0))
covC21_00   <- mxAlgebra(name = "C21_00", expression = (c0 + Mod2%x%cMod0) %*% t(c0+ Mod1%x%cMod0))
covE21_00   <- mxAlgebra(name = "E21_00", expression = (e0 + Mod2%x%eMod0) %*% t(e0+ Mod1%x%eMod0))

covA21_11   <- mxAlgebra(name = "A21_11", expression = (a1 + Mod2%x%aMod1) %*% t(a1+ Mod1%x%aMod1))
covC21_11   <- mxAlgebra(name = "C21_11", expression = (c1 + Mod2%x%cMod1) %*% t(c1+ Mod1%x%cMod1))
covE21_11   <- mxAlgebra(name = "E21_11", expression = (e1 + Mod2%x%eMod1) %*% t(e1+ Mod1%x%eMod1))

covA21_01   <- mxAlgebra(name = "A21_01", expression = (a0 + Mod2%x%aMod0) %*% t(a1+ Mod1%x%aMod1))
covC21_01   <- mxAlgebra(name = "C21_01", expression = (c0 + Mod2%x%cMod0) %*% t(c1+ Mod1%x%cMod1))
covE21_01   <- mxAlgebra(name = "E21_01", expression = (e0 + Mod2%x%eMod0) %*% t(e1+ Mod1%x%eMod1))

covA21_0m   <- mxAlgebra(name = "A21_0m", expression = (a0 + Mod2%x%aModm) %*% t(a0+ Mod1%x%aModm))
covC21_0m   <- mxAlgebra(name = "C21_0m", expression = (c0 + Mod2%x%cModm) %*% t(c0+ Mod1%x%cModm))
covE21_0m   <- mxAlgebra(name = "E21_0m", expression = (e0 + Mod2%x%eModm) %*% t(e0+ Mod1%x%eModm))

covA21_1m   <- mxAlgebra(name = "A21_1m", expression = (a1 + Mod2%x%aMod1) %*% t(am+ Mod1%x%aModm))
covC21_1m   <- mxAlgebra(name = "C21_1m", expression = (c1 + Mod2%x%cMod1) %*% t(cm+ Mod1%x%cModm))
covE21_1m   <- mxAlgebra(name = "E21_1m", expression = (e1 + Mod2%x%eMod1) %*% t(em+ Mod1%x%eModm))

covA21_mm   <- mxAlgebra(name = "A21_mm", expression = (am + Mod2%x%aModm) %*% t(am+ Mod1%x%aModm))
covC21_mm   <- mxAlgebra(name = "C21_mm", expression = (cm + Mod2%x%cModm) %*% t(cm+ Mod1%x%cModm))
covE21_mm   <- mxAlgebra(name = "E21_mm", expression = (em + Mod2%x%eModm) %*% t(em+ Mod1%x%eModm))


varA2_00    <- mxAlgebra(name = "A2_00",  expression = (a0 + Mod2%x%aMod0) %*% t(a0+ Mod2%x%aMod0))
varC2_00    <- mxAlgebra(name = "C2_00",  expression = (c0 + Mod2%x%cMod0) %*% t(c0+ Mod2%x%cMod0))
varE2_00    <- mxAlgebra(name = "E2_00",  expression = (e0 + Mod2%x%eMod0) %*% t(e0+ Mod2%x%eMod0))

varA2_11    <- mxAlgebra(name = "A2_11",  expression = (a1 + Mod2%x%aMod1) %*% t(a1+ Mod2%x%aMod1))
varC2_11    <- mxAlgebra(name = "C2_11",  expression = (c1 + Mod2%x%cMod1) %*% t(c1+ Mod2%x%cMod1))
varE2_11    <- mxAlgebra(name = "E2_11",  expression = (e1 + Mod2%x%eMod1) %*% t(e1+ Mod2%x%eMod1))

varA2_01    <- mxAlgebra(name = "A2_01",  expression = (a0 + Mod2%x%aMod0) %*% t(a1+ Mod2%x%aMod1))
varC2_01    <- mxAlgebra(name = "C2_01",  expression = (c0 + Mod2%x%cMod0) %*% t(c1+ Mod2%x%cMod1))
varE2_01    <- mxAlgebra(name = "E2_01",  expression = (e0 + Mod2%x%eMod0) %*% t(e1+ Mod2%x%eMod1))

varA2_0m    <- mxAlgebra(name = "A2_0m",  expression = (a0 + Mod2%x%aMod0) %*% t(am+ Mod2%x%aModm))
varC2_0m    <- mxAlgebra(name = "C2_0m",  expression = (c0 + Mod2%x%cMod0) %*% t(cm+ Mod2%x%cModm))
varE2_0m    <- mxAlgebra(name = "E2_0m",  expression = (e0 + Mod2%x%eMod0) %*% t(em+ Mod2%x%eModm))

varA2_1m    <- mxAlgebra(name = "A2_1m",  expression = (a1 + Mod2%x%aMod1) %*% t(am+ Mod2%x%aModm))
varC2_1m    <- mxAlgebra(name = "C2_1m",  expression = (c1 + Mod2%x%cMod1) %*% t(cm+ Mod2%x%cModm))
varE2_1m    <- mxAlgebra(name = "E2_1m",  expression = (e1 + Mod2%x%eMod1) %*% t(em+ Mod2%x%eModm))

varA2_mm    <- mxAlgebra(name = "A2_mm",  expression = (am + Mod2%x%aModm) %*% t(am+ Mod2%x%aModm))
varC2_mm    <- mxAlgebra(name = "C2_mm",  expression = (cm + Mod2%x%cModm) %*% t(cm+ Mod2%x%cModm))
varE2_mm    <- mxAlgebra(name = "E2_mm",  expression = (em + Mod2%x%eModm) %*% t(em+ Mod2%x%eModm))

# Create a dummy algebra containing A, C, E for selected values of the moderator, such that CIs can be computed for them
myVec         <- mxMatrix("Full",5,1,values=c(-30, -20, -10, 0, 10),name="myVec")
myUnit        <- mxMatrix("Unit",5,1,name="myUnit")
myVarAValuesI  <- mxAlgebra(expression=(((myUnit %x% a0[2,1] + myVec %x% aMod0[2,1]) * (myUnit %x% a0[2,1] + myVec %x% aMod0[2,1])) +
                                        ((myUnit %x% a0[2,2] + myVec %x% aMod0[2,2]) * (myUnit %x% a0[2,2] + myVec %x% aMod0[2,2]))), name="myVarAValuesI")
myVarAValuesD  <- mxAlgebra(expression=(((myUnit %x% a1[2,1] + myVec %x% aMod1[2,1]) * (myUnit %x% a1[2,1] + myVec %x% aMod1[2,1])) +
                                        ((myUnit %x% a1[2,2] + myVec %x% aMod1[2,2]) * (myUnit %x% a1[2,2] + myVec %x% aMod1[2,2]))), name="myVarAValuesD")

# Matrices generated to hold A, C, and E computed Variance Components, for unmoderated models
varA0_00   <- mxAlgebra(name = "A0_00",  expression = a0 %*% t(a0))
varC0_00   <- mxAlgebra(name = "C0_00",  expression = c0 %*% t(c0))
varE0_00   <- mxAlgebra(name = "E0_00",  expression = e0 %*% t(e0))

varA0_11   <- mxAlgebra(name = "A0_11",  expression = a1 %*% t(a1))
varC0_11   <- mxAlgebra(name = "C0_11",  expression = c1 %*% t(c1))
varE0_11   <- mxAlgebra(name = "E0_11",  expression = e1 %*% t(e1))

varA0_01   <- mxAlgebra(name = "A0_01",  expression = a0 %*% t(a1))
varC0_01   <- mxAlgebra(name = "C0_01",  expression = c0 %*% t(c1))
varE0_01   <- mxAlgebra(name = "E0_01",  expression = e0 %*% t(e1))

varA0_0m   <- mxAlgebra(name = "A0_0m",  expression = a0 %*% t(am))
varC0_0m   <- mxAlgebra(name = "C0_0m",  expression = c0 %*% t(cm))
varE0_0m   <- mxAlgebra(name = "E0_0m",  expression = e0 %*% t(em))

varA0_1m   <- mxAlgebra(name = "A0_1m",  expression = a1 %*% t(am))
varC0_1m   <- mxAlgebra(name = "C0_1m",  expression = c1 %*% t(cm))
varE0_1m   <- mxAlgebra(name = "E0_1m",  expression = e1 %*% t(em))

varA0_mm   <- mxAlgebra(name = "A0_mm",  expression = am %*% t(am))
varC0_mm   <- mxAlgebra(name = "C0_mm",  expression = cm %*% t(cm))
varE0_mm   <- mxAlgebra(name = "E0_mm",  expression = em %*% t(em))

# Algebra to compute total variances and standard deviations (diagonal only) per twin
var1_00     <- mxAlgebra( A1_00+C1_00+E1_00, name="V1_00" )
var1_11     <- mxAlgebra( A1_11+C1_11+E1_11, name="V1_11" )
var1_01     <- mxAlgebra( A1_01+C1_01+E1_01, name="V1_01" )
var1_0m     <- mxAlgebra( A1_0m+C1_0m+E1_0m, name="V1_0m" )
var1_1m     <- mxAlgebra( A1_1m+C1_1m+E1_1m, name="V1_1m" )
var1_mm     <- mxAlgebra( A1_mm+C1_mm+E1_mm, name="V1_mm" )

var2_00     <- mxAlgebra( A2_00+C2_00+E2_00, name="V2_00" )
var2_11     <- mxAlgebra( A2_11+C2_11+E2_11, name="V2_11" )
var2_01     <- mxAlgebra( A2_01+C2_01+E2_01, name="V2_01" )
var2_0m     <- mxAlgebra( A2_0m+C2_0m+E2_0m, name="V2_0m" )
var2_1m     <- mxAlgebra( A2_1m+C2_1m+E2_1m, name="V2_1m" )
var2_mm     <- mxAlgebra( A2_mm+C2_mm+E2_mm, name="V2_mm" )

var0_00    <- mxAlgebra( A0_00+C0_00+E0_00, name="V0_00" )
var0_11    <- mxAlgebra( A0_11+C0_11+E0_11, name="V0_11" )
var0_01    <- mxAlgebra( A0_01+C0_01+E0_01, name="V0_01" )
var0_0m    <- mxAlgebra( A0_0m+C0_0m+E0_0m, name="V0_0m" )
var0_1m    <- mxAlgebra( A0_1m+C0_1m+E0_1m, name="V0_1m" )
var0_mm    <- mxAlgebra( A0_mm+C0_mm+E0_mm, name="V0_mm" )

# Standardized variance components for unmoderated model
VarA_imm    <- mxAlgebra( A0_00/V0_00      , name = "A_imm")
VarA_del    <- mxAlgebra( A0_11/V0_11      , name = "A_del")
VarC_imm    <- mxAlgebra( C0_00/V0_00      , name = "C_imm")
VarC_del    <- mxAlgebra( C0_11/V0_11      , name = "C_del")
VarE_imm    <- mxAlgebra( E0_00/V0_00      , name = "E_imm")
VarE_del    <- mxAlgebra( E0_11/V0_11      , name = "E_del")

# Algebra for expected variance/covariance matrices
expCovMZ00 <- mxAlgebra(name = "expCovMZ00", expression = rbind (cbind(A1_00+C1_00+E1_00  , A12_00+C12_00),
                                                                 cbind(A21_00+C21_00      , A2_00+C2_00+E2_00)))

expCovDZ00 <- mxAlgebra(name = "expCovDZ00", expression = rbind (cbind(A1_00+C1_00+E1_00  , 0.5%x%A12_00+C12_00),
                                                                 cbind(0.5%x%A21_00+C21_00, A2_00+C2_00+E2_00))) 


expCovMZ11 <- mxAlgebra(name = "expCovMZ11", expression = rbind (cbind(A1_11+C1_11+E1_11  , A12_11+C12_11),
                                                                 cbind(A21_11+C21_11      , A2_11+C2_11+E2_11)))

expCovDZ11 <- mxAlgebra(name = "expCovDZ11", expression = rbind (cbind(A1_11+C1_11+E1_11  , 0.5%x%A12_11+C12_11),
                                                                 cbind(0.5%x%A21_11+C21_11, A2_11+C2_11+E2_11))) 


expCovMZ01 <- mxAlgebra(name = "expCovMZ01", expression = rbind (cbind(A1_01+C1_01+E1_01  , A12_01+C12_01),
                                                                 cbind(A21_01+C21_01      , A2_01+C2_01+E2_01)))

expCovDZ01 <- mxAlgebra(name = "expCovDZ01", expression = rbind (cbind(A1_01+C1_01+E1_01  , 0.5%x%A12_01+C12_01),
                                                                 cbind(0.5%x%A21_01+C21_01, A2_01+C2_01+E2_01))) 


expCovMZ0m <- mxAlgebra(name = "expCovMZ0m", expression = rbind (cbind(A1_0m+C1_0m+E1_0m  , A12_0m+C12_0m),
                                                                 cbind(A21_0m+C21_0m      , A2_0m+C2_0m+E2_0m)))

expCovDZ0m <- mxAlgebra(name = "expCovDZ0m", expression = rbind (cbind(A1_0m+C1_0m+E1_0m  , 0.5%x%A12_0m+C12_0m),
                                                                 cbind(0.5%x%A21_0m+C21_0m, A2_0m+C2_0m+E2_0m))) 


expCovMZ1m <- mxAlgebra(name = "expCovMZ1m", expression = rbind (cbind(A1_1m+C1_1m+E1_1m  , A12_1m+C12_1m),
                                                                 cbind(A21_1m+C21_1m      , A2_1m+C2_1m+E2_1m)))

expCovDZ1m <- mxAlgebra(name = "expCovDZ1m", expression = rbind (cbind(A1_1m+C1_1m+E1_1m  , 0.5%x%A12_1m+C12_1m),
                                                                 cbind(0.5%x%A21_1m+C21_1m, A2_1m+C2_1m+E2_1m))) 


expCovMZmm <- mxAlgebra(name = "expCovMZmm", expression = rbind (cbind(A1_mm+C1_mm+E1_mm  , A12_mm+C12_mm),
                                                                 cbind(A21_mm+C21_mm      , A2_mm+C2_mm+E2_mm)))

expCovDZmm <- mxAlgebra(name = "expCovDZmm", expression = rbind (cbind(A1_mm+C1_mm+E1_mm  , 0.5%x%A12_mm+C12_mm),
                                                                 cbind(0.5%x%A21_mm+C21_mm, A2_mm+C2_mm+E2_mm))) 

# Create Data Objects for Multiple Groups
dataMZ00    <- mxData( observed=mz00Data, type="raw" )
dataDZ00    <- mxData( observed=dz00Data, type="raw" )
dataMZ11    <- mxData( observed=mz11Data, type="raw" )
dataDZ11    <- mxData( observed=dz11Data, type="raw" )
dataMZ01    <- mxData( observed=mz01Data, type="raw" )
dataDZ01    <- mxData( observed=dz01Data, type="raw" )
dataMZ0m    <- mxData( observed=mz0mData, type="raw" )
dataDZ0m    <- mxData( observed=dz0mData, type="raw" )
dataMZ1m    <- mxData( observed=mz1mData, type="raw" )
dataDZ1m    <- mxData( observed=dz1mData, type="raw" )
dataMZmm    <- mxData( observed=mzmmData, type="raw" )
dataDZmm    <- mxData( observed=dzmmData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ00     <- mxExpectationNormal( covariance="expCovMZ00", means="meanG00", dimnames=selVars )
expDZ00     <- mxExpectationNormal( covariance="expCovDZ00", means="meanG00", dimnames=selVars )
expMZ11     <- mxExpectationNormal( covariance="expCovMZ11", means="meanG11", dimnames=selVars )
expDZ11     <- mxExpectationNormal( covariance="expCovDZ11", means="meanG11", dimnames=selVars )
expMZ01     <- mxExpectationNormal( covariance="expCovMZ01", means="meanG01", dimnames=selVars )
expDZ01     <- mxExpectationNormal( covariance="expCovDZ01", means="meanG01", dimnames=selVars )
expMZ0m     <- mxExpectationNormal( covariance="expCovMZ0m", means="meanG0m", dimnames=selVars )
expDZ0m     <- mxExpectationNormal( covariance="expCovDZ0m", means="meanG0m", dimnames=selVars )
expMZ1m     <- mxExpectationNormal( covariance="expCovMZ1m", means="meanG1m", dimnames=selVars )
expDZ1m     <- mxExpectationNormal( covariance="expCovDZ1m", means="meanG1m", dimnames=selVars )
expMZmm     <- mxExpectationNormal( covariance="expCovMZmm", means="meanGmm", dimnames=selVars )
expDZmm     <- mxExpectationNormal( covariance="expCovDZmm", means="meanGmm", dimnames=selVars )
funML       <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
paths0      <- list( meanG0, pathA0, pathC0, pathE0, modPathA0, modPathC0, modPathE0, myVec, myUnit)
paths1      <- list( meanG1, pathA1, pathC1, pathE1, modPathA1, modPathC1, modPathE1, myVec, myUnit)
pathsm      <- list( meanGm, pathAm, pathCm, pathEm, modPathAm, modPathCm, modPathEm)

pars00      <- list( bMale, bYob,  paths0,                             
                     varA1_00, varC1_00, varE1_00, covA12_00, covC12_00, covE12_00, 
                     varA2_00, varC2_00, varE2_00, covA21_00, covC21_00, covE21_00,
                     var1_00, var2_00, varA0_00, varC0_00, varE0_00, var0_00, VarA_imm, VarC_imm, VarE_imm, myVarAValuesI )

pars11      <- list( bMale, bYob, paths1,                             
                     varA1_11, varC1_11, varE1_11, covA12_11, covC12_11, covE12_11, 
                     varA2_11, varC2_11, varE2_11, covA21_11, covC21_11, covE21_11,
                     var1_11, var2_11, varA0_11, varC0_11, varE0_11, var0_11, VarA_del, VarC_del, VarE_del, myVarAValuesD )

pars01      <- list( bMale, bYob, paths0, paths1,                     
                     varA1_01, varC1_01, varE1_01, covA12_01, covC12_01, covE12_01, 
                     varA2_01, varC2_01, varE2_01, covA21_01, covC21_01, covE21_01,
                     var1_01, var2_01, varA0_01, varC0_01, varE0_01, var0_01  )

pars0m      <- list( bMale, bYob, paths0, pathsm,                     
                     varA1_0m, varC1_0m, varE1_0m, covA12_0m, covC12_0m, covE12_0m, 
                     varA2_0m, varC2_0m, varE2_0m, covA21_0m, covC21_0m, covE21_0m,
                     var1_0m, var2_0m, varA0_0m, varC0_0m, varE0_0m, var0_0m   )

pars1m      <- list( bMale, bYob, paths1, pathsm,                     
                     varA1_1m, varC1_1m, varE1_1m, covA12_1m, covC12_1m, covE12_1m, 
                     varA2_1m, varC2_1m, varE2_1m, covA21_1m, covC21_1m, covE21_1m,
                     var1_1m, var2_1m, varA0_1m, varC0_1m, varE0_1m, var0_1m   )

parsmm      <- list( bMale, bYob, pathsm,                             
                     varA1_mm, varC1_mm, varE1_mm, covA12_mm, covC12_mm, covE12_mm, 
                     varA2_mm, varC2_mm, varE2_mm, covA21_mm, covC21_mm, covE21_mm,
                     var1_mm, var2_mm, varA0_mm, varC0_mm, varE0_mm, var0_mm   )

defs        <- list( Male, Yob, mod_tw1, mod_tw2 )

modelMZ00   <- mxModel( pars00, meanG00, defs, expCovMZ00, dataMZ00, expMZ00, funML, name="MZ00" )
modelMZ11   <- mxModel( pars11, meanG11, defs, expCovMZ11, dataMZ11, expMZ11, funML, name="MZ11" )
modelMZ01   <- mxModel( pars01, meanG01, defs, expCovMZ01, dataMZ01, expMZ01, funML, name="MZ01" )
modelMZ0m   <- mxModel( pars0m, meanG0m, defs, expCovMZ0m, dataMZ0m, expMZ0m, funML, name="MZ0m" )
modelMZ1m   <- mxModel( pars1m, meanG1m, defs, expCovMZ1m, dataMZ1m, expMZ1m, funML, name="MZ1m" )
modelMZmm   <- mxModel( parsmm, meanGmm, defs, expCovMZmm, dataMZmm, expMZmm, funML, name="MZmm" )

modelDZ00   <- mxModel( pars00, meanG00, defs, expCovDZ00, dataDZ00, expDZ00, funML, name="DZ00" )
modelDZ11   <- mxModel( pars11, meanG11, defs, expCovDZ11, dataDZ11, expDZ11, funML, name="DZ11" )
modelDZ01   <- mxModel( pars01, meanG01, defs, expCovDZ01, dataDZ01, expDZ01, funML, name="DZ01" )
modelDZ0m   <- mxModel( pars0m, meanG0m, defs, expCovDZ0m, dataDZ0m, expDZ0m, funML, name="DZ0m" )
modelDZ1m   <- mxModel( pars1m, meanG1m, defs, expCovDZ1m, dataDZ1m, expDZ1m, funML, name="DZ1m" )
modelDZmm   <- mxModel( parsmm, meanGmm, defs, expCovDZmm, dataDZmm, expDZmm, funML, name="DZmm" )

multi       <- mxFitFunctionMultigroup( c("MZ00","DZ00", "MZ11", "DZ11", "MZ01", "DZ01", "MZ0m", "DZ0m", "MZ1m", "DZ1m", "MZmm", "DZmm") )

# Create Confidence Interval Objects
ciEST        <- list(varA0_00, varC0_00, varE0_00, var0_00, VarA_imm, VarC_imm, VarE_imm, pathA0, pathC0, pathE0,
                     varA0_11, varC0_11, varE0_11, var0_11, VarA_del, VarC_del, VarE_del, pathA1, pathC1, pathE1 )
ciACE        <- mxCI( c("A_imm[2,2]", "A_del[2,2]", "C_imm[2,2]", "C_del[2,2]", "E_imm[2,2]", "E_del[2,2]") )

# Build Model with Confidence Intervals
modelFull  <- mxModel( "Bivariate Multigroup",
                       modelMZ00, modelDZ00, modelMZ11, modelDZ11, modelMZ01, modelDZ01,
                       modelMZ0m, modelDZ0m, modelMZ1m, modelDZ1m, modelMZmm, modelDZmm, multi, ciEST, ciACE )


# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run bivariate ACE Model: paths allowed to differ between immediate and delayed tracking; no moderation by cito
# = Table 1, Model 1b
fitFull    <- mxTryHard( modelFull, intervals=T )
sumFull    <- summary( fitFull )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitFull)
fitEstCis(fitFull)

# 
# # ----------------------------------------------------------------------------------------------------------------------
# # RUN SUBMODELS
# 

## Constrain immediate and delayed tracking paths to be equal
# = Table 1, Model 1a
ciEST2     <- list(varA0_00, varC0_00, varE0_00, var0_00, VarA_imm, VarC_imm, VarE_imm, pathA0, pathC0, pathE0 )
ciACE2     <- mxCI( c("A_imm[1,1]", "A_imm[2,2]", "C_imm[1,1]", "C_imm[2,2]", "E_imm[1,1]", "E_imm[2,2]") )
modelBiv   <- mxModel( "Bivariate",
                       modelMZ00, modelDZ00, modelMZ11, modelDZ11, modelMZ01, modelDZ01,
                       modelMZ0m, modelDZ0m, modelMZ1m, modelDZ1m, modelMZmm, modelDZmm, multi, ciEST2, ciACE2 )
modelBiv   <- omxSetParameters( modelBiv, labels = c(meanLabs0, meanLabs1, meanLabsm), free=TRUE, values=svM, newlabels=meanLabs )
modelBiv   <- omxSetParameters( modelBiv, labels = c(aLabs0, aLabs1, aLabsm), free=TRUE, values=svPa, newlabels=aLabs )
modelBiv   <- omxSetParameters( modelBiv, labels = c(cLabs0, cLabs1, cLabsm), free=TRUE, values=svPc, newlabels=cLabs )
modelBiv   <- omxSetParameters( modelBiv, labels = c(eLabs0, eLabs1, eLabsm), free=TRUE, values=svPe, newlabels=eLabs )
fitBiv     <- mxTryHard( modelBiv, intervals=T )
sumBiv     <- summary(fitBiv)
fitGofs(fitBiv); fitEstCis(fitBiv)

## Constrain immediate tracking and missing info on tracking paths to be equal
# Extra check, not in paper
model1_M   <- mxModel( fitFull, name="Biv Multigr Miss = Imm" )
model1_M   <- omxSetParameters( model1_M, labels = c(meanLabs0, meanLabsm), free=TRUE, values=svM, newlabels=meanLabs )
model1_M   <- omxSetParameters( model1_M, labels = c(aLabs0, aLabsm), free=TRUE, values=svPa, newlabels=aLabs )
model1_M   <- omxSetParameters( model1_M, labels = c(cLabs0, cLabsm), free=TRUE, values=svPc, newlabels=cLabs )
model1_M   <- omxSetParameters( model1_M, labels = c(eLabs0, eLabsm), free=TRUE, values=svPe, newlabels=eLabs )
fit1_M     <- mxTryHard( model1_M, intervals=F )
fitGofs(fit1_M); fitEsts(fit1_M, 3)

## Allow moderators to be estimated in full model
# = Table 2, Model 2
modelFullMod<- mxModel( fitFull, name="Bivariate Multigroup Moderated" )
modelFullMod<- omxSetParameters( modelFullMod, labels=c(aModLabs0, aModLabs1, aModLabsm, cModLabs0, cModLabs1, cModLabsm,
                                                        eModLabs0, eModLabs1, eModLabsm), free=c(F,T,T), values=0 )
fitFullMod  <- mxTryHard( modelFullMod, intervals=F )
sumFullMod  <- summary( fitFullMod )
fitGofs(fitFullMod); fitEsts(fitFullMod, 3)

## Allow moderators to be estimated in full model, except constrain c_21 and c_22 in delayed to be non-negative 0
# See note to supplementary Figure 1
twinDataNM2   <- subset(twinData, male_1!="NA" & male_2!="NA" & c_yob!="NA" & c_cito_1!="NA" & c_cito_2!="NA", useVars)
min_cito      <- min(twinDataNM2$c_cito_2, twinDataNM2$c_cito_1)
max_cito      <- max(twinDataNM2$c_cito_2, twinDataNM2$c_cito_1)
c1_constr1min <- mxConstraint(c1[2,1]>-min_cito%x%cMod1[2,1], name = "c1_constr1min")
c1_constr1max <- mxConstraint(c1[2,1]>-max_cito%x%cMod1[2,1], name = "c1_constr1max")
c1_constr2min <- mxConstraint(c1[2,2]>-min_cito%x%cMod1[2,2], name = "c1_constr2min")
c1_constr2max <- mxConstraint(c1[2,2]>-max_cito%x%cMod1[2,2], name = "c1_constr2max")

pars11c       <- list( pars11, c1_constr1min, c1_constr1max, c1_constr2min, c1_constr2max)
pars01c       <- list( pars01, c1_constr1min, c1_constr1max, c1_constr2min, c1_constr2max  )
pars1mc       <- list( pars1m, c1_constr1min, c1_constr1max, c1_constr2min, c1_constr2max)

modelMZ11c    <- mxModel( pars11c, meanG11, defs, expCovMZ11, dataMZ11, expMZ11, funML, name="MZ11c" )
modelMZ01c    <- mxModel( pars01c, meanG01, defs, expCovMZ01, dataMZ01, expMZ01, funML, name="MZ01c" )
modelDZ11c    <- mxModel( pars11c, meanG11, defs, expCovDZ11, dataDZ11, expDZ11, funML, name="DZ11c" )
modelDZ01c    <- mxModel( pars01c, meanG01, defs, expCovDZ01, dataDZ01, expDZ01, funML, name="DZ01c" )

multic        <- mxFitFunctionMultigroup( c("MZ00","DZ00", "MZ11c", "DZ11c", "MZ01c", "DZ01c", "MZ0m", "DZ0m",                 "MZmm", "DZmm") )

modelFullModc <- mxModel( "Biv Multi Mod Constrained",
                        modelMZ00, modelDZ00, modelMZ11c, modelDZ11c, modelMZ01c, modelDZ01c,
                        modelMZ0m, modelDZ0m, modelMZ1m, modelDZ1m, modelMZmm, modelDZmm, multic)
modelFullModc <- omxSetParameters( modelFullModc, labels=c(aModLabs0, aModLabs1, aModLabsm, cModLabs0, cModLabs1, cModLabsm,
                                                           eModLabs0, eModLabs1, eModLabsm), free=c(F,T,T), values=0 )
fitFullModc   <- mxTryHard( modelFullModc, intervals=F )
fitEsts(fitFullModc, 3)

## Constrain immediate and delayed tracking moderating paths to be equal
# Model not shown but referred to at end of section "The impact of delayed tracking moderated by performance"
modelIntMod <- mxModel( fitFullMod, name="Intermediate Moderated" )
modelIntMod <- omxSetParameters( modelIntMod, labels = c(aModLabs0, aModLabs1), free=c(F,T,T), values=0, newlabels=aModLabs )
modelIntMod <- omxSetParameters( modelIntMod, labels = c(cModLabs0, cModLabs1), free=c(F,T,T), values=0, newlabels=cModLabs )
modelIntMod <- omxSetParameters( modelIntMod, labels = c(eModLabs0, eModLabs1), free=c(F,T,T), values=0, newlabels=eModLabs )
fitIntMod   <- mxTryHard( modelIntMod, intervals=F )
fitGofs(fitIntMod); fitEsts(fitIntMod, 3)

## Constrain also moderating paths for those with tracking info missing to be equal
# Extra check, not in paper
modelIntMod2 <- mxModel( fitIntMod, name="Intermediate Moderated 2" )
modelIntMod2 <- omxSetParameters( modelIntMod2, labels = c(aModLabsm, aModLabs), free=c(F,T,T), values=0, newlabels=aModLabs )
modelIntMod2 <- omxSetParameters( modelIntMod2, labels = c(cModLabsm, cModLabs), free=c(F,T,T), values=0, newlabels=cModLabs )
modelIntMod2 <- omxSetParameters( modelIntMod2, labels = c(eModLabsm, eModLabs), free=c(F,T,T), values=0, newlabels=eModLabs )
fitIntMod2   <- mxTryHard( modelIntMod2, intervals=F )
fitGofs(fitIntMod2); fitEsts(fitIntMod2, 3)

## Constrain immediate and delayed tracking main & moderating paths to be equal
# Extra check, not in paper
modelBivMod <- mxModel( fitIntMod2, name="Bivariate Moderated" )
modelBivMod <- omxSetParameters( modelBivMod, labels = c(meanLabs0, meanLabs1, meanLabsm), free=TRUE, values=svM, newlabels=meanLabs )
modelBivMod <- omxSetParameters( modelBivMod, labels = c(aLabs0, aLabs1, aLabsm), free=TRUE, values=svPa, newlabels=aLabs )
modelBivMod <- omxSetParameters( modelBivMod, labels = c(cLabs0, cLabs1, cLabsm), free=TRUE, values=svPc, newlabels=cLabs )
modelBivMod <- omxSetParameters( modelBivMod, labels = c(eLabs0, eLabs1, eLabsm), free=TRUE, values=svPe, newlabels=eLabs )
fitBivMod   <- mxTryHard( modelBivMod, intervals=F )
fitGofs(fitBivMod); fitEsts(fitBivMod, 3)


## Compare fit between models

#  Model with different paths for immediate and delayed tracking VS same paths for both
mxCompare(fitFull, fitBiv)

#  Model with different paths for immediate and missing tracking VS same paths for both
mxCompare(fitFull, fit1_M)

# Model with different moderated paths for delayed and immediate tracking VS unmoderated paths (while differentating between imm. vs del.)
mxCompare(fitFullMod, fitFull)

# Model with full moderation VS moderator of common c path constrained to 0
mxCompare(fitFullMod, fitFullModc)

# Model with different moderated paths for delayed and immediate tracking VS same moderated paths for both
mxCompare(fitFullMod, fitIntMod)
mxCompare(fitFullMod, fitIntMod2) # missing moderating path also same
mxCompare(fitFullModc, fitIntMod)


# # ----------------------------------------------------------------------------------------------------------------------
# # PRINT COEFFICIENTS
# 
#

# Create output table 1 - non-moderated
x1 <-       rbind(sumBiv$parameters[14,'name'],
                  sumBiv$parameters[15,'name'],
                  sumBiv$parameters[ 1,'name'],
                  sumBiv$parameters[ 4,'name'],
                  sumBiv$parameters[ 7,'name'],
                  "Vax",
                  "Vcx",
                  "Vex",
                  "Vtx",
                  sumBiv$parameters[ 2,'name'],
                  sumBiv$parameters[ 5,'name'],
                  sumBiv$parameters[ 8,'name'],
                  sumBiv$parameters[ 3,'name'],
                  sumBiv$parameters[ 6,'name'],
                  sumBiv$parameters[ 9,'name'],
                  "Vay",
                  "Vcy",
                  "Vey",
                  "Vty")
x2 <- round(rbind(sumBiv$parameters[14,'Estimate'],
                  sumBiv$parameters[15,'Estimate'],
                  sumBiv$parameters[ 1,'Estimate'],
                  sumBiv$parameters[ 4,'Estimate'],
                  sumBiv$parameters[ 7,'Estimate'],
                  mxEval(A0_00[1,1],fitBiv),
                  mxEval(C0_00[1,1],fitBiv),
                  mxEval(E0_00[1,1],fitBiv),
                  mxEval(V0_00[1,1],fitBiv),
                  sumBiv$parameters[ 2,'Estimate'],
                  sumBiv$parameters[ 5,'Estimate'],
                  sumBiv$parameters[ 8,'Estimate'],
                  sumBiv$parameters[ 3,'Estimate'],
                  sumBiv$parameters[ 6,'Estimate'],
                  sumBiv$parameters[ 9,'Estimate'],
                  mxEval(A0_00[2,2],fitBiv),
                  mxEval(C0_00[2,2],fitBiv),
                  mxEval(E0_00[2,2],fitBiv),
                  mxEval(V0_00[2,2],fitBiv)), 2)
x3 <- round(rbind(sumBiv$parameters[14,'Std.Error'],
                  sumBiv$parameters[15,'Std.Error'],
                  sumBiv$parameters[ 1,'Std.Error'],
                  sumBiv$parameters[ 4,'Std.Error'],
                  sumBiv$parameters[ 7,'Std.Error'],
                  0,
                  0,
                  0,
                  0,
                  sumBiv$parameters[ 2,'Std.Error'],
                  sumBiv$parameters[ 5,'Std.Error'],
                  sumBiv$parameters[ 8,'Std.Error'],
                  sumBiv$parameters[ 3,'Std.Error'],
                  sumBiv$parameters[ 6,'Std.Error'],
                  sumBiv$parameters[ 9,'Std.Error'],
                  0,
                  0,
                  0,
                  0), 2)
x4 <-       rbind(sumFull$parameters[23,'name'],
                  sumFull$parameters[24,'name'],
                  sumFull$parameters[ 1,'name'],
                  sumFull$parameters[ 4,'name'],
                  sumFull$parameters[ 7,'name'],
                  "Vax",
                  "Vcx",
                  "Vex",
                  "Vtx",
                  sumFull$parameters[ 2,'name'],
                  sumFull$parameters[ 5,'name'],
                  sumFull$parameters[ 8,'name'],
                  sumFull$parameters[ 3,'name'],
                  sumFull$parameters[ 6,'name'],
                  sumFull$parameters[ 9,'name'],
                  "Vay",
                  "Vcy",
                  "Vey",
                  "Vty")
x5 <- round(rbind(sumFull$parameters[23,'Estimate'],
                  sumFull$parameters[24,'Estimate'],
                  sumFull$parameters[ 1,'Estimate'],
                  sumFull$parameters[ 4,'Estimate'],
                  sumFull$parameters[ 7,'Estimate'],
                  mxEval(A0_00[1,1],fitFull),
                  mxEval(C0_00[1,1],fitFull),
                  mxEval(E0_00[1,1],fitFull),
                  mxEval(V0_00[1,1],fitFull),
                  sumFull$parameters[ 2,'Estimate'],
                  sumFull$parameters[ 5,'Estimate'],
                  sumFull$parameters[ 8,'Estimate'],
                  sumFull$parameters[ 3,'Estimate'],
                  sumFull$parameters[ 6,'Estimate'],
                  sumFull$parameters[ 9,'Estimate'],
                  mxEval(A0_00[2,2],fitFull),
                  mxEval(C0_00[2,2],fitFull),
                  mxEval(E0_00[2,2],fitFull),
                  mxEval(V0_00[2,2],fitFull)), 2)
x6 <- round(rbind(sumFull$parameters[23,'Std.Error'],
                  sumFull$parameters[24,'Std.Error'],
                  sumFull$parameters[ 1,'Std.Error'],
                  sumFull$parameters[ 4,'Std.Error'],
                  sumFull$parameters[ 7,'Std.Error'],
                  0,
                  0,
                  0,
                  0,
                  sumFull$parameters[ 2,'Std.Error'],
                  sumFull$parameters[ 5,'Std.Error'],
                  sumFull$parameters[ 8,'Std.Error'],
                  sumFull$parameters[ 3,'Std.Error'],
                  sumFull$parameters[ 6,'Std.Error'],
                  sumFull$parameters[ 9,'Std.Error'],
                  0,
                  0,
                  0,
                  0), 2)
x7 <-       rbind(sumFull$parameters[25,'name'],
                  sumFull$parameters[26,'name'],
                  sumFull$parameters[10,'name'],
                  sumFull$parameters[13,'name'],
                  sumFull$parameters[16,'name'],
                  "Vax",
                  "Vcx",
                  "Vex",
                  "Vtx",
                  sumFull$parameters[11,'name'],
                  sumFull$parameters[14,'name'],
                  sumFull$parameters[17,'name'],
                  sumFull$parameters[12,'name'],
                  sumFull$parameters[15,'name'],
                  sumFull$parameters[18,'name'],
                  "Vay",
                  "Vcy",
                  "Vey",
                  "Vty")
x8 <- round(rbind(sumFull$parameters[25,'Estimate'],
                  sumFull$parameters[26,'Estimate'],
                  sumFull$parameters[10,'Estimate'],
                  sumFull$parameters[13,'Estimate'],
                  sumFull$parameters[16,'Estimate'],
                  mxEval(A0_11[1,1],fitFull),
                  mxEval(C0_11[1,1],fitFull),
                  mxEval(E0_11[1,1],fitFull),
                  mxEval(V0_11[1,1],fitFull),
                  sumFull$parameters[11,'Estimate'],
                  sumFull$parameters[14,'Estimate'],
                  sumFull$parameters[17,'Estimate'],
                  sumFull$parameters[12,'Estimate'],
                  sumFull$parameters[15,'Estimate'],
                  sumFull$parameters[18,'Estimate'],
                  mxEval(A0_11[2,2],fitFull),
                  mxEval(C0_11[2,2],fitFull),
                  mxEval(E0_11[2,2],fitFull),
                  mxEval(V0_11[2,2],fitFull)), 2)
x9 <- round(rbind(sumFull$parameters[25,'Std.Error'],
                  sumFull$parameters[26,'Std.Error'],
                  sumFull$parameters[10,'Std.Error'],
                  sumFull$parameters[13,'Std.Error'],
                  sumFull$parameters[16,'Std.Error'],
                  0,
                  0,
                  0,
                  0,
                  sumFull$parameters[11,'Std.Error'],
                  sumFull$parameters[14,'Std.Error'],
                  sumFull$parameters[17,'Std.Error'],
                  sumFull$parameters[12,'Std.Error'],
                  sumFull$parameters[15,'Std.Error'],
                  sumFull$parameters[18,'Std.Error'],
                  0,
                  0,
                  0,
                  0), 2)
Table2a <- data.frame(par=x1, est=x2, se=x3, par=x4, est=x5, se=x6, par=x7, est=x8, se=x9)
Table2a
write.table(Table2a, file = "TablesFigures/Table2_04.csv", sep = ",", quote = FALSE, row.names=F)


# Create output table 2 - Moderation
x10 <-      rbind(sumFullMod$parameters[23,'name'],
                  sumFullMod$parameters[24,'name'],
                  sumFullMod$parameters[ 1,'name'],
                  sumFullMod$parameters[ 4,'name'],
                  sumFullMod$parameters[ 7,'name'],
                  sumFullMod$parameters[ 2,'name'],
                  sumFullMod$parameters[25,'name'],
                  sumFullMod$parameters[ 5,'name'],
                  sumFullMod$parameters[27,'name'],
                  sumFullMod$parameters[ 8,'name'],
                  sumFullMod$parameters[29,'name'],
                  sumFullMod$parameters[ 3,'name'],
                  sumFullMod$parameters[26,'name'],
                  sumFullMod$parameters[ 6,'name'],
                  sumFullMod$parameters[28,'name'],
                  sumFullMod$parameters[ 9,'name'],
                  sumFullMod$parameters[30,'name'])
x11 <-round(rbind(sumFullMod$parameters[23,'Estimate'],
                  sumFullMod$parameters[24,'Estimate'],
                  sumFullMod$parameters[ 1,'Estimate'],
                  sumFullMod$parameters[ 4,'Estimate'],
                  sumFullMod$parameters[ 7,'Estimate'],
                  sumFullMod$parameters[ 2,'Estimate'],
                  sumFullMod$parameters[25,'Estimate'],
                  sumFullMod$parameters[ 5,'Estimate'],
                  sumFullMod$parameters[27,'Estimate'],
                  sumFullMod$parameters[ 8,'Estimate'],
                  sumFullMod$parameters[29,'Estimate'],
                  sumFullMod$parameters[ 3,'Estimate'],
                  sumFullMod$parameters[26,'Estimate'],
                  sumFullMod$parameters[ 6,'Estimate'],
                  sumFullMod$parameters[28,'Estimate'],
                  sumFullMod$parameters[ 9,'Estimate'],
                  sumFullMod$parameters[30,'Estimate']), 2)
x12 <-round(rbind(sumFullMod$parameters[23,'Std.Error'],
                  sumFullMod$parameters[24,'Std.Error'],
                  sumFullMod$parameters[ 1,'Std.Error'],
                  sumFullMod$parameters[ 4,'Std.Error'],
                  sumFullMod$parameters[ 7,'Std.Error'],
                  sumFullMod$parameters[ 2,'Std.Error'],
                  sumFullMod$parameters[25,'Std.Error'],
                  sumFullMod$parameters[ 5,'Std.Error'],
                  sumFullMod$parameters[27,'Std.Error'],
                  sumFullMod$parameters[ 8,'Std.Error'],
                  sumFullMod$parameters[29,'Std.Error'],
                  sumFullMod$parameters[ 3,'Std.Error'],
                  sumFullMod$parameters[26,'Std.Error'],
                  sumFullMod$parameters[ 6,'Std.Error'],
                  sumFullMod$parameters[28,'Std.Error'],
                  sumFullMod$parameters[ 9,'Std.Error'],
                  sumFullMod$parameters[30,'Std.Error']), 2)
x13 <-      rbind(sumFullMod$parameters[31,'name'],
                  sumFullMod$parameters[32,'name'],
                  sumFullMod$parameters[10,'name'],
                  sumFullMod$parameters[13,'name'],
                  sumFullMod$parameters[16,'name'],
                  sumFullMod$parameters[11,'name'],
                  sumFullMod$parameters[33,'name'],
                  sumFullMod$parameters[14,'name'],
                  sumFullMod$parameters[35,'name'],
                  sumFullMod$parameters[17,'name'],
                  sumFullMod$parameters[37,'name'],
                  sumFullMod$parameters[12,'name'],
                  sumFullMod$parameters[34,'name'],
                  sumFullMod$parameters[15,'name'],
                  sumFullMod$parameters[36,'name'],
                  sumFullMod$parameters[18,'name'],
                  sumFullMod$parameters[38,'name'])
x14 <-round(rbind(sumFullMod$parameters[31,'Estimate'],
                  sumFullMod$parameters[32,'Estimate'],
                  sumFullMod$parameters[10,'Estimate'],
                  sumFullMod$parameters[13,'Estimate'],
                  sumFullMod$parameters[16,'Estimate'],
                  sumFullMod$parameters[11,'Estimate'],
                  sumFullMod$parameters[33,'Estimate'],
                  sumFullMod$parameters[14,'Estimate'],
                  sumFullMod$parameters[35,'Estimate'],
                  sumFullMod$parameters[17,'Estimate'],
                  sumFullMod$parameters[37,'Estimate'],
                  sumFullMod$parameters[12,'Estimate'],
                  sumFullMod$parameters[34,'Estimate'],
                  sumFullMod$parameters[15,'Estimate'],
                  sumFullMod$parameters[36,'Estimate'],
                  sumFullMod$parameters[18,'Estimate'],
                  sumFullMod$parameters[38,'Estimate']), 2)
x15 <-round(rbind(sumFullMod$parameters[31,'Std.Error'],
                  sumFullMod$parameters[31,'Std.Error'],
                  sumFullMod$parameters[10,'Std.Error'],
                  sumFullMod$parameters[13,'Std.Error'],
                  sumFullMod$parameters[16,'Std.Error'],
                  sumFullMod$parameters[11,'Std.Error'],
                  sumFullMod$parameters[33,'Std.Error'],
                  sumFullMod$parameters[14,'Std.Error'],
                  sumFullMod$parameters[35,'Std.Error'],
                  sumFullMod$parameters[17,'Std.Error'],
                  sumFullMod$parameters[37,'Std.Error'],
                  sumFullMod$parameters[12,'Std.Error'],
                  sumFullMod$parameters[34,'Std.Error'],
                  sumFullMod$parameters[15,'Std.Error'],
                  sumFullMod$parameters[36,'Std.Error'],
                  sumFullMod$parameters[18,'Std.Error'],
                  sumFullMod$parameters[38,'Std.Error']), 2)
Table3 <- data.frame(par=x10, est=x11, se=x12, par=x13, est=x14, se=x15)
Table3
write.table(Table3, file = "TablesFigures/Table3_04.csv", sep = ",", quote = FALSE, row.names=F)


## No moderation

# All
round(mxEval(a0,fitBiv$MZ00),2)
round(mxEval(a0^2,fitBiv$MZ00),2)
round(mxEval(A1_00,fitBiv$MZ00),2)
round(mxEval(c0,fitBiv$MZ00),2)
round(mxEval(c0^2,fitBiv$MZ00),2)
round(mxEval(C1_00,fitBiv$MZ00),2)
round(mxEval(e0,fitBiv$MZ00),2)
round(mxEval(e0^2,fitBiv$MZ00),2)
round(mxEval(E1_00,fitBiv$MZ00),2)
round(mxEval(V1_00,fitBiv$MZ00),2)
round(mxEval(A1_00/V1_00,fitBiv$MZ00),2)
round(mxEval(C1_00/V1_00,fitBiv$MZ00),2)
round(mxEval(E1_00/V1_00,fitBiv$MZ00),2)

round(mxEval(a0[2,1]^2/(a0[2,1]^2+a0[2,2]^2),fitBiv$MZ00),2)
round(mxEval(a0[2,2]^2/(a0[2,1]^2+a0[2,2]^2),fitBiv$MZ00),2)

# Immediate
round(mxEval(a0,fitFull$MZ00),2)
round(mxEval(a0^2,fitFull$MZ00),2)
round(mxEval(A1_00,fitFull$MZ00),2)
round(mxEval(c0,fitFull$MZ00),2)
round(mxEval(c0^2,fitFull$MZ00),2)
round(mxEval(C1_00,fitFull$MZ00),2)
round(mxEval(e0,fitFull$MZ00),2)
round(mxEval(e0^2,fitFull$MZ00),2)
round(mxEval(E1_00,fitFull$MZ00),2)
round(mxEval(V1_00,fitFull$MZ00),2)
round(mxEval(A1_00/V1_00,fitFull$MZ00),2)
round(mxEval(C1_00/V1_00,fitFull$MZ00),2)
round(mxEval(E1_00/V1_00,fitFull$MZ00),2)

round(mxEval(a0[2,1]^2/(a0[2,1]^2+a0[2,2]^2),fitFull$MZ00),2)
round(mxEval(a0[2,2]^2/(a0[2,1]^2+a0[2,2]^2),fitFull$MZ00),2)

round(mxEval(a0[2,1]^2/(a0[2,1]^2+c0[2,1]^2+e0[2,1]^2),fitFull$MZ00),2)
round(mxEval(c0[2,1]^2/(a0[2,1]^2+c0[2,1]^2+e0[2,1]^2),fitFull$MZ00),2)
round(mxEval(e0[2,1]^2/(a0[2,1]^2+c0[2,1]^2+e0[2,1]^2),fitFull$MZ00),2)
round(mxEval(a0[2,2]^2/(a0[2,2]^2+c0[2,2]^2+e0[2,2]^2),fitFull$MZ00),2)
round(mxEval(c0[2,2]^2/(a0[2,2]^2+c0[2,2]^2+e0[2,2]^2),fitFull$MZ00),2)
round(mxEval(e0[2,2]^2/(a0[2,2]^2+c0[2,2]^2+e0[2,2]^2),fitFull$MZ00),2)

# Delayed
round(mxEval(a1,fitFull$MZ11),2)
round(mxEval(a1^2,fitFull$MZ11),2)
round(mxEval(A1_11,fitFull$MZ11),2)
round(mxEval(c1,fitFull$MZ11),2)
round(mxEval(c1^2,fitFull$MZ11),2)
round(mxEval(C1_11,fitFull$MZ11),2)
round(mxEval(e1,fitFull$MZ11),2)
round(mxEval(e1^2,fitFull$MZ11),2)
round(mxEval(E1_11,fitFull$MZ11),2)
round(mxEval(V1_11,fitFull$MZ11),2)
round(mxEval(A1_11/V1_11,fitFull$MZ11),2)
round(mxEval(C1_11/V1_11,fitFull$MZ11),2)
round(mxEval(E1_11/V1_11,fitFull$MZ11),2)

round(mxEval(a1[2,1]^2/(a1[2,1]^2+a1[2,2]^2),fitFull$MZ11),2)
round(mxEval(a1[2,2]^2/(a1[2,1]^2+a1[2,2]^2),fitFull$MZ11),2)

round(mxEval(a1[2,1]^2/(a1[2,1]^2+c1[2,1]^2+e1[2,1]^2),fitFull$MZ11),2)
round(mxEval(c1[2,1]^2/(a1[2,1]^2+c1[2,1]^2+e1[2,1]^2),fitFull$MZ11),2)
round(mxEval(e1[2,1]^2/(a1[2,1]^2+c1[2,1]^2+e1[2,1]^2),fitFull$MZ11),2)
round(mxEval(a1[2,2]^2/(a1[2,2]^2+c1[2,2]^2+e1[2,2]^2),fitFull$MZ11),2)
round(mxEval(c1[2,2]^2/(a1[2,2]^2+c1[2,2]^2+e1[2,2]^2),fitFull$MZ11),2)
round(mxEval(e1[2,2]^2/(a1[2,2]^2+c1[2,2]^2+e1[2,2]^2),fitFull$MZ11),2)

# No information on tracking
round(mxEval(A1_mm/V1_mm,fitFull$MZmm),2)
round(mxEval(C1_mm/V1_mm,fitFull$MZmm),2)
round(mxEval(E1_mm/V1_mm,fitFull$MZmm),2)


#=======================================================================#
#   PLOT MODERATION    FULL MODEL                                       #
#=======================================================================#

## Bar plot of immediate vs delayed (Figure 1)

#  Create dataframe of results
ACE_val <- c(sumFull$CI$estimate[1], sumFull$CI$estimate[3], sumFull$CI$estimate[5], 
             sumFull$CI$estimate[2], sumFull$CI$estimate[4], sumFull$CI$estimate[6])
ACE_lb  <- c(sumFull$CI$lbound[1]  , sumFull$CI$lbound[3]  , sumFull$CI$lbound[5], 
             sumFull$CI$lbound[2]  , sumFull$CI$lbound[4]  , sumFull$CI$lbound[6])
ACE_ub  <- c(sumFull$CI$ubound[1]  , sumFull$CI$ubound[3]  , sumFull$CI$ubound[5], 
             sumFull$CI$ubound[2]  , sumFull$CI$ubound[4]  , sumFull$CI$ubound[6])
x_val   <- c("A","C","E", "A", "C", "E")
gr_val  <- factor(c("Immediate", "Immediate", "Immediate", "Delayed", "Delayed", "Delayed"), levels = c("Immediate", "Delayed"))
df_res  <- data.frame(gr_val, x_val, ACE_val, ACE_lb, ACE_ub)
head(df_res)

# Create theme settings for boxplot
My_Theme = theme_minimal() +  
  theme(legend.position="bottom",
        legend.title = element_text(size = 14, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        legend.text  = element_text(size=12),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.x  = element_text(size=12),
        axis.text.y  = element_text(size=10))

# Create plot
p <- ggplot(data=df_res, aes(x=x_val, y=ACE_val, fill=factor(gr_val))) +
     geom_col(position = position_dodge(0.8), width = 0.8)   +
     geom_errorbar(aes(ymin=ACE_lb, ymax=ACE_ub), width=.2,
                   position=position_dodge(0.8)) 

# Add layout settings to the plot
p + My_Theme + 
  labs(fill ="Tracking:", y = "Proportion of Variance") + 
  scale_fill_brewer(palette="Paired") + 
  coord_cartesian(ylim = c(0, 1)) 
ggsave("TablesFigures/Bar_ImmediateVsDelayed_ACE_04.pdf", width = 4, height = 4)


## Graphs of Delayed vs Immediate moderated by Performance 

# Create Arrays of results
selDefs = names(mzData)[1:3]
mzdef1 = as.vector(mzData[, selDefs[1]])
mzdef2 = as.vector(mzData[, selDefs[3]])
dzdef1 = as.vector(dzData[, selDefs[1]])
dzdef2 = as.vector(dzData[, selDefs[3]])
allValuesOfDefVar= c(mzdef1, mzdef2, dzdef1, dzdef2)
cito_sc = sort(unique(allValuesOfDefVar))
cito_length=length(cito_sc)

# In case you want uncentered version of cito on the x-axis
selDefs2 = names(twinDataNM)[13:14]
def1 = as.vector(twinDataNM[, selDefs2[1]])
def2 = as.vector(twinDataNM[, selDefs2[2]])
allValuesOfDefVar2= c(def1, def2)
cito_score = sort(unique(allValuesOfDefVar2))

# Purcell's representation (total variances)
# Immediate = Delayed
Va=array(NA, dim=c(nv,nv,cito_length))
Vc=array(NA, dim=c(nv,nv,cito_length))
Ve=array(NA, dim=c(nv,nv,cito_length))
Vt=array(NA, dim=c(nv,nv,cito_length))
SVa=array(NA, dim=c(nv,nv,cito_length))
SVc=array(NA, dim=c(nv,nv,cito_length))
SVe=array(NA, dim=c(nv,nv,cito_length))
corA = array(NA, dim=c(nv,nv,cito_length))
corC = array(NA, dim=c(nv,nv,cito_length))
corE = array(NA, dim=c(nv,nv,cito_length))
corP = array(NA, dim=c(nv,nv,cito_length))
amod = array(NA, dim=c(nv,nv,cito_length))
cmod = array(NA, dim=c(nv,nv,cito_length))
emod = array(NA, dim=c(nv,nv,cito_length))

# Immediate
Va_imm=array(NA, dim=c(nv,nv,cito_length))
Vc_imm=array(NA, dim=c(nv,nv,cito_length))
Ve_imm=array(NA, dim=c(nv,nv,cito_length))
Vt_imm=array(NA, dim=c(nv,nv,cito_length))
SVa_imm=array(NA, dim=c(nv,nv,cito_length))
SVc_imm=array(NA, dim=c(nv,nv,cito_length))
SVe_imm=array(NA, dim=c(nv,nv,cito_length))
corA_imm = array(NA, dim=c(nv,nv,cito_length))
corC_imm = array(NA, dim=c(nv,nv,cito_length))
corE_imm = array(NA, dim=c(nv,nv,cito_length))
corP_imm = array(NA, dim=c(nv,nv,cito_length))
amod_imm = array(NA, dim=c(nv,nv,cito_length))
cmod_imm = array(NA, dim=c(nv,nv,cito_length))
emod_imm = array(NA, dim=c(nv,nv,cito_length))

# Delayed
Va_del=array(NA, dim=c(nv,nv,cito_length))
Vc_del=array(NA, dim=c(nv,nv,cito_length))
Ve_del=array(NA, dim=c(nv,nv,cito_length))
Vt_del=array(NA, dim=c(nv,nv,cito_length))
SVa_del=array(NA, dim=c(nv,nv,cito_length))
SVc_del=array(NA, dim=c(nv,nv,cito_length))
SVe_del=array(NA, dim=c(nv,nv,cito_length))
corA_del = array(NA, dim=c(nv,nv,cito_length))
corC_del = array(NA, dim=c(nv,nv,cito_length))
corE_del = array(NA, dim=c(nv,nv,cito_length))
corP_del = array(NA, dim=c(nv,nv,cito_length))
amod_del = array(NA, dim=c(nv,nv,cito_length))
cmod_del = array(NA, dim=c(nv,nv,cito_length))
emod_del = array(NA, dim=c(nv,nv,cito_length))

for (i in 1:cito_length) {
  Va[,,i]   <- mxEval((a0 + aMod0%x%cito_sc[i]) %*% t(a0 + aMod0%x%cito_sc[i]), fitBivMod$MZ00)
  Vc[,,i]   <- mxEval((c0 + cMod0%x%cito_sc[i]) %*% t(c0 + cMod0%x%cito_sc[i]), fitBivMod$MZ00)
  Ve[,,i]   <- mxEval((e0 + eMod0%x%cito_sc[i]) %*% t(e0 + eMod0%x%cito_sc[i]), fitBivMod$MZ00)
  Vt[,,i]   <- Va[,,i] + Vc[,,i] + Ve[,,i]
  SVa[,,i]  <- Va[,,i]/Vt[,,i]
  SVc[,,i]  <- Vc[,,i]/Vt[,,i]
  SVe[,,i]  <- Ve[,,i]/Vt[,,i]
  corA[,,i] <- solve(sqrt(diag(2)*Va[,,i]))%*%Va[,,i]%*%solve(sqrt(diag(2)*Va[,,i]))
  corC[,,i] <- solve(sqrt(diag(2)*Vc[,,i]))%*%Vc[,,i]%*%solve(sqrt(diag(2)*Vc[,,i]))
  corE[,,i] <- solve(sqrt(diag(2)*Ve[,,i]))%*%Ve[,,i]%*%solve(sqrt(diag(2)*Ve[,,i]))
  corP[,,i] <- solve(sqrt(diag(2)*Vt[,,i]))%*%Vt[,,i]%*%solve(sqrt(diag(2)*Vt[,,i]))
  amod[,,i] <- mxEval((a0 + aMod0%x%cito_sc[i]), fitBivMod$MZ00)
  cmod[,,i] <- mxEval((c0 + cMod0%x%cito_sc[i]), fitBivMod$MZ00)
  emod[,,i] <- mxEval((e0 + eMod0%x%cito_sc[i]), fitBivMod$MZ00)
  
  Va_imm[,,i]   <- mxEval((a0 + aMod0%x%cito_sc[i]) %*% t(a0 + aMod0%x%cito_sc[i]), fitFullMod$MZ00)
  Vc_imm[,,i]   <- mxEval((c0 + cMod0%x%cito_sc[i]) %*% t(c0 + cMod0%x%cito_sc[i]), fitFullMod$MZ00)
  Ve_imm[,,i]   <- mxEval((e0 + eMod0%x%cito_sc[i]) %*% t(e0 + eMod0%x%cito_sc[i]), fitFullMod$MZ00)
  Vt_imm[,,i]   <- Va_imm[,,i] + Vc_imm[,,i] + Ve_imm[,,i]
  SVa_imm[,,i]  <- Va_imm[,,i]/Vt_imm[,,i]
  SVc_imm[,,i]  <- Vc_imm[,,i]/Vt_imm[,,i]
  SVe_imm[,,i]  <- Ve_imm[,,i]/Vt_imm[,,i]
  corA_imm[,,i] <- solve(sqrt(diag(2)*Va_imm[,,i]))%*%Va_imm[,,i]%*%solve(sqrt(diag(2)*Va_imm[,,i]))
  corC_imm[,,i] <- solve(sqrt(diag(2)*Vc_imm[,,i]))%*%Vc_imm[,,i]%*%solve(sqrt(diag(2)*Vc_imm[,,i]))
  corE_imm[,,i] <- solve(sqrt(diag(2)*Ve_imm[,,i]))%*%Ve_imm[,,i]%*%solve(sqrt(diag(2)*Ve_imm[,,i]))
  corP_imm[,,i] <- solve(sqrt(diag(2)*Vt_imm[,,i]))%*%Vt_imm[,,i]%*%solve(sqrt(diag(2)*Vt_imm[,,i]))
  amod_imm[,,i] <- mxEval((a0 + aMod0%x%cito_sc[i]), fitFullMod$MZ00)
  cmod_imm[,,i] <- mxEval((c0 + cMod0%x%cito_sc[i]), fitFullMod$MZ00)
  emod_imm[,,i] <- mxEval((e0 + eMod0%x%cito_sc[i]), fitFullMod$MZ00)
  
  Va_del[,,i]   <- mxEval((a1 + aMod1%x%cito_sc[i]) %*% t(a1 + aMod1%x%cito_sc[i]), fitFullMod$MZ11)
  Vc_del[,,i]   <- mxEval((c1 + cMod1%x%cito_sc[i]) %*% t(c1 + cMod1%x%cito_sc[i]), fitFullMod$MZ11)
  Ve_del[,,i]   <- mxEval((e1 + eMod1%x%cito_sc[i]) %*% t(e1 + eMod1%x%cito_sc[i]), fitFullMod$MZ11)
  Vt_del[,,i]   <- Va_del[,,i] + Vc_del[,,i] + Ve_del[,,i]
  SVa_del[,,i]  <- Va_del[,,i]/Vt_del[,,i]
  SVc_del[,,i]  <- Vc_del[,,i]/Vt_del[,,i]
  SVe_del[,,i]  <- Ve_del[,,i]/Vt_del[,,i]
  corA_del[,,i] <- solve(sqrt(diag(2)*Va_del[,,i]))%*%Va_del[,,i]%*%solve(sqrt(diag(2)*Va_del[,,i]))
  corC_del[,,i] <- solve(sqrt(diag(2)*Vc_del[,,i]))%*%Vc_del[,,i]%*%solve(sqrt(diag(2)*Vc_del[,,i]))
  corE_del[,,i] <- solve(sqrt(diag(2)*Ve_del[,,i]))%*%Ve_del[,,i]%*%solve(sqrt(diag(2)*Ve_del[,,i]))
  corP_del[,,i] <- solve(sqrt(diag(2)*Vt_del[,,i]))%*%Vt_del[,,i]%*%solve(sqrt(diag(2)*Vt_del[,,i]))
  amod_del[,,i] <- mxEval((a1 + aMod1%x%cito_sc[i]), fitFullMod$MZ11)
  cmod_del[,,i] <- mxEval((c1 + cMod1%x%cito_sc[i]), fitFullMod$MZ11)
  emod_del[,,i] <- mxEval((e1 + eMod1%x%cito_sc[i]), fitFullMod$MZ11)
}

out <- as.matrix(cbind(Va[2,2,],Vc[2,2,],Ve[2,2,],Vt[2,2,],SVa[2,2,],SVc[2,2,],SVe[2,2,], corA[2,1,],corC[2,1,],corE[2,1,], corP[2,1,], 
                       amod[2,1,], amod[2,2,], cmod[2,1,], cmod[2,2,], emod[2,1,], emod[2,2,]))
names(out) = c('Va', 'Vc', 'Ve', 'Vt', 'SVa', 'SVc', 'SVe', 'corA', 'corC', 'corE', 'corP', 'amodC', 'amodU', 'cmodC', 'cmodU', 'emodC', 'emodU')
head(out)

out_imm <- as.matrix(cbind(Va_imm[2,2,],Vc_imm[2,2,],Ve_imm[2,2,],Vt_imm[2,2,],SVa_imm[2,2,],SVc_imm[2,2,],SVe_imm[2,2,], corA_imm[2,1,],corC_imm[2,1,],corE_imm[2,1,], corP_imm[2,1,], 
                           amod_imm[2,1,], amod_imm[2,2,], cmod_imm[2,1,], cmod_imm[2,2,], emod_imm[2,1,], emod_imm[2,2,]))
names(out_imm) = c('Va_imm', 'Vc_imm', 'Ve_imm', 'Vt_imm', 'SVa_imm', 'SVc_imm', 'SVe_imm', 'corA_imm', 'corC_imm', 'corE_imm', 'corP_imm', 'amodC_imm', 'amodU_imm', 'cmodC_imm', 'cmodU_imm', 'emodC_imm', 'emodU_imm')
head(out_imm)

out_del <- as.matrix(cbind(Va_del[2,2,],Vc_del[2,2,],Ve_del[2,2,],Vt_del[2,2,],SVa_del[2,2,],SVc_del[2,2,],SVe_del[2,2,], corA_del[2,1,],corC_del[2,1,],corE_del[2,1,], corP_del[2,1,],
                           amod_del[2,1,], amod_del[2,2,], cmod_del[2,1,], cmod_del[2,2,], emod_del[2,1,], emod_del[2,2,]))
names(out_del) = c('Va_del', 'Vc_del', 'Ve_del', 'Vt_del', 'SVa_del', 'SVc_del', 'SVe_del', 'corA_del', 'corC_del', 'corE_del', 'corP_del', 'amodC_del', 'amodU_del', 'cmodC_del', 'cmodU_del', 'emodC_del', 'emodU_del')
head(out_del)

out_A <- cbind(out_imm[,1], out_del[,1], out_imm[,5], out_del[,5])
out_C <- cbind(out_imm[,2], out_del[,2], out_imm[,6], out_del[,6])
out_E <- cbind(out_imm[,3], out_del[,3], out_imm[,7], out_del[,7])


## Same but now for common and unique variants components separately instead of summed

## Common variance 

# Immediate
Va_immC=array(NA, dim=c(nv,nv,cito_length))
Vc_immC=array(NA, dim=c(nv,nv,cito_length))
Ve_immC=array(NA, dim=c(nv,nv,cito_length))
Vt_immC=array(NA, dim=c(nv,nv,cito_length))
SVa_immC=array(NA, dim=c(nv,nv,cito_length))
SVc_immC=array(NA, dim=c(nv,nv,cito_length))
SVe_immC=array(NA, dim=c(nv,nv,cito_length))

# Delayed
Va_delC=array(NA, dim=c(nv,nv,cito_length))
Vc_delC=array(NA, dim=c(nv,nv,cito_length))
Ve_delC=array(NA, dim=c(nv,nv,cito_length))
Vt_delC=array(NA, dim=c(nv,nv,cito_length))
SVa_delC=array(NA, dim=c(nv,nv,cito_length))
SVc_delC=array(NA, dim=c(nv,nv,cito_length))
SVe_delC=array(NA, dim=c(nv,nv,cito_length))

for (i in 1:cito_length) {
  Va_immC[,,i]   <- mxEval((a0[2,1] + aMod0[2,1]%x%cito_sc[i]) %*% t(a0[2,1] + aMod0[2,1]%x%cito_sc[i]), fitFullMod$MZ00)
  Vc_immC[,,i]   <- mxEval((c0[2,1] + cMod0[2,1]%x%cito_sc[i]) %*% t(c0[2,1] + cMod0[2,1]%x%cito_sc[i]), fitFullMod$MZ00)
  Ve_immC[,,i]   <- mxEval((e0[2,1] + eMod0[2,1]%x%cito_sc[i]) %*% t(e0[2,1] + eMod0[2,1]%x%cito_sc[i]), fitFullMod$MZ00)
  Vt_immC[,,i]   <- Va_immC[,,i] + Vc_immC[,,i] + Ve_immC[,,i]
  SVa_immC[,,i]  <- Va_immC[,,i]/Vt_immC[,,i]
  SVc_immC[,,i]  <- Vc_immC[,,i]/Vt_immC[,,i]
  SVe_immC[,,i]  <- Ve_immC[,,i]/Vt_immC[,,i]
  
  Va_delC[,,i]   <- mxEval((a1[2,1] + aMod1[2,1]%x%cito_sc[i]) %*% t(a1[2,1] + aMod1[2,1]%x%cito_sc[i]), fitFullMod$MZ11)
  Vc_delC[,,i]   <- mxEval((c1[2,1] + cMod1[2,1]%x%cito_sc[i]) %*% t(c1[2,1] + cMod1[2,1]%x%cito_sc[i]), fitFullMod$MZ11)
  Ve_delC[,,i]   <- mxEval((e1[2,1] + eMod1[2,1]%x%cito_sc[i]) %*% t(e1[2,1] + eMod1[2,1]%x%cito_sc[i]), fitFullMod$MZ11)
  Vt_delC[,,i]   <- Va_delC[,,i] + Vc_delC[,,i] + Ve_delC[,,i]
  SVa_delC[,,i]  <- Va_delC[,,i]/Vt_delC[,,i]
  SVc_delC[,,i]  <- Vc_delC[,,i]/Vt_delC[,,i]
  SVe_delC[,,i]  <- Ve_delC[,,i]/Vt_delC[,,i]
}

out_immC <- as.matrix(cbind(Va_immC[2,2,],Vc_immC[2,2,],Ve_immC[2,2,],Vt_immC[2,2,],SVa_immC[2,2,],SVc_immC[2,2,],SVe_immC[2,2,]))
names(out_immC) = c('Va_immC', 'Vc_immC', 'Ve_immC', 'Vt_immC', 'SVa_immC', 'SVc_immC', 'SVe_immC')
head(out_immC)

out_delC <- as.matrix(cbind(Va_delC[2,2,],Vc_delC[2,2,],Ve_delC[2,2,],Vt_delC[2,2,],SVa_delC[2,2,],SVc_delC[2,2,],SVe_delC[2,2,]))
names(out_delC) = c('Va_delC', 'Vc_delC', 'Ve_delC', 'Vt_delC', 'SVa_delC', 'SVc_delC', 'SVe_delC')
head(out_delC)

out_AC <- cbind(out_immC[,1], out_delC[,1], out_immC[,5], out_delC[,5])
out_CC <- cbind(out_immC[,2], out_delC[,2], out_immC[,6], out_delC[,6])
out_EC <- cbind(out_immC[,3], out_delC[,3], out_immC[,7], out_delC[,7])


## Unique variance
# Immediate
Va_immU=array(NA, dim=c(nv,nv,cito_length))
Vc_immU=array(NA, dim=c(nv,nv,cito_length))
Ve_immU=array(NA, dim=c(nv,nv,cito_length))
Vt_immU=array(NA, dim=c(nv,nv,cito_length))
SVa_immU=array(NA, dim=c(nv,nv,cito_length))
SVc_immU=array(NA, dim=c(nv,nv,cito_length))
SVe_immU=array(NA, dim=c(nv,nv,cito_length))

# Delayed
Va_delU=array(NA, dim=c(nv,nv,cito_length))
Vc_delU=array(NA, dim=c(nv,nv,cito_length))
Ve_delU=array(NA, dim=c(nv,nv,cito_length))
Vt_delU=array(NA, dim=c(nv,nv,cito_length))
SVa_delU=array(NA, dim=c(nv,nv,cito_length))
SVc_delU=array(NA, dim=c(nv,nv,cito_length))
SVe_delU=array(NA, dim=c(nv,nv,cito_length))

for (i in 1:cito_length) {
  Va_immU[,,i]   <- mxEval((a0[2,2] + aMod0[2,2]%x%cito_sc[i]) %*% t(a0[2,2] + aMod0[2,2]%x%cito_sc[i]), fitFullMod$MZ00)
  Vc_immU[,,i]   <- mxEval((c0[2,2] + cMod0[2,2]%x%cito_sc[i]) %*% t(c0[2,2] + cMod0[2,2]%x%cito_sc[i]), fitFullMod$MZ00)
  Ve_immU[,,i]   <- mxEval((e0[2,2] + eMod0[2,2]%x%cito_sc[i]) %*% t(e0[2,2] + eMod0[2,2]%x%cito_sc[i]), fitFullMod$MZ00)
  Vt_immU[,,i]   <- Va_immU[,,i] + Vc_immU[,,i] + Ve_immU[,,i]
  SVa_immU[,,i]  <- Va_immU[,,i]/Vt_immU[,,i]
  SVc_immU[,,i]  <- Vc_immU[,,i]/Vt_immU[,,i]
  SVe_immU[,,i]  <- Ve_immU[,,i]/Vt_immU[,,i]
  
  Va_delU[,,i]   <- mxEval((a1[2,2] + aMod1[2,2]%x%cito_sc[i]) %*% t(a1[2,2] + aMod1[2,2]%x%cito_sc[i]), fitFullMod$MZ11)
  Vc_delU[,,i]   <- mxEval((c1[2,2] + cMod1[2,2]%x%cito_sc[i]) %*% t(c1[2,2] + cMod1[2,2]%x%cito_sc[i]), fitFullMod$MZ11)
  Ve_delU[,,i]   <- mxEval((e1[2,2] + eMod1[2,2]%x%cito_sc[i]) %*% t(e1[2,2] + eMod1[2,2]%x%cito_sc[i]), fitFullMod$MZ11)
  Vt_delU[,,i]   <- Va_delU[,,i] + Vc_delU[,,i] + Ve_delU[,,i]
  SVa_delU[,,i]  <- Va_delU[,,i]/Vt_delU[,,i]
  SVc_delU[,,i]  <- Vc_delU[,,i]/Vt_delU[,,i]
  SVe_delU[,,i]  <- Ve_delU[,,i]/Vt_delU[,,i]
}

out_immU <- as.matrix(cbind(Va_immU[2,2,],Vc_immU[2,2,],Ve_immU[2,2,],Vt_immU[2,2,],SVa_immU[2,2,],SVc_immU[2,2,],SVe_immU[2,2,]))
names(out_immU) = c('Va_immU', 'Vc_immU', 'Ve_immU', 'Vt_immU', 'SVa_immU', 'SVc_immU', 'SVe_immU')
head(out_immU)

out_delU <- as.matrix(cbind(Va_delU[2,2,],Vc_delU[2,2,],Ve_delU[2,2,],Vt_delU[2,2,],SVa_delU[2,2,],SVc_delU[2,2,],SVe_delU[2,2,]))
names(out_delU) = c('Va_delU', 'Vc_delU', 'Ve_delU', 'Vt_delU', 'SVa_delU', 'SVc_delU', 'SVe_delU')
head(out_delU)

out_AU <- cbind(out_immU[,1], out_delU[,1], out_immU[,5], out_delU[,5])
out_CU <- cbind(out_immU[,2], out_delU[,2], out_immU[,6], out_delU[,6])
out_EU <- cbind(out_immU[,3], out_delU[,3], out_immU[,7], out_delU[,7])


## Same but now for moderation model with common c path constrained to 0

# Immediate
Va_imm2=array(NA, dim=c(nv,nv,cito_length))
Vc_imm2=array(NA, dim=c(nv,nv,cito_length))
Ve_imm2=array(NA, dim=c(nv,nv,cito_length))
Vt_imm2=array(NA, dim=c(nv,nv,cito_length))
SVa_imm2=array(NA, dim=c(nv,nv,cito_length))
SVc_imm2=array(NA, dim=c(nv,nv,cito_length))
SVe_imm2=array(NA, dim=c(nv,nv,cito_length))
corA_imm2 = array(NA, dim=c(nv,nv,cito_length))
corC_imm2 = array(NA, dim=c(nv,nv,cito_length))
corE_imm2 = array(NA, dim=c(nv,nv,cito_length))
corP_imm2 = array(NA, dim=c(nv,nv,cito_length))
amod_imm2 = array(NA, dim=c(nv,nv,cito_length))
cmod_imm2 = array(NA, dim=c(nv,nv,cito_length))
emod_imm2 = array(NA, dim=c(nv,nv,cito_length))

# Delayed
Va_del2=array(NA, dim=c(nv,nv,cito_length))
Vc_del2=array(NA, dim=c(nv,nv,cito_length))
Ve_del2=array(NA, dim=c(nv,nv,cito_length))
Vt_del2=array(NA, dim=c(nv,nv,cito_length))
SVa_del2=array(NA, dim=c(nv,nv,cito_length))
SVc_del2=array(NA, dim=c(nv,nv,cito_length))
SVe_del2=array(NA, dim=c(nv,nv,cito_length))
corA_del2 = array(NA, dim=c(nv,nv,cito_length))
corC_del2 = array(NA, dim=c(nv,nv,cito_length))
corE_del2 = array(NA, dim=c(nv,nv,cito_length))
corP_del2 = array(NA, dim=c(nv,nv,cito_length))
amod_del2 = array(NA, dim=c(nv,nv,cito_length))
cmod_del2 = array(NA, dim=c(nv,nv,cito_length))
emod_del2 = array(NA, dim=c(nv,nv,cito_length))

for (i in 2:cito_length) {
  Va_imm2[,,i]   <- mxEval((a0 + aMod0%x%cito_sc[i]) %*% t(a0 + aMod0%x%cito_sc[i]), fitFullModc$MZ00)
  Vc_imm2[,,i]   <- mxEval((c0 + cMod0%x%cito_sc[i]) %*% t(c0 + cMod0%x%cito_sc[i]), fitFullModc$MZ00)
  Ve_imm2[,,i]   <- mxEval((e0 + eMod0%x%cito_sc[i]) %*% t(e0 + eMod0%x%cito_sc[i]), fitFullModc$MZ00)
  Vt_imm2[,,i]   <- Va_imm2[,,i] + Vc_imm2[,,i] + Ve_imm2[,,i]
  SVa_imm2[,,i]  <- Va_imm2[,,i]/Vt_imm2[,,i]
  SVc_imm2[,,i]  <- Vc_imm2[,,i]/Vt_imm2[,,i]
  SVe_imm2[,,i]  <- Ve_imm2[,,i]/Vt_imm2[,,i]
  corA_imm2[,,i] <- solve(sqrt(diag(2)*Va_imm2[,,i]))%*%Va_imm2[,,i]%*%solve(sqrt(diag(2)*Va_imm2[,,i]))
  corC_imm2[,,i] <- solve(sqrt(diag(2)*Vc_imm2[,,i]))%*%Vc_imm2[,,i]%*%solve(sqrt(diag(2)*Vc_imm2[,,i]))
  corE_imm2[,,i] <- solve(sqrt(diag(2)*Ve_imm2[,,i]))%*%Ve_imm2[,,i]%*%solve(sqrt(diag(2)*Ve_imm2[,,i]))
  corP_imm2[,,i] <- solve(sqrt(diag(2)*Vt_imm2[,,i]))%*%Vt_imm2[,,i]%*%solve(sqrt(diag(2)*Vt_imm2[,,i]))
  amod_imm2[,,i] <- mxEval((a0 + aMod0%x%cito_sc[i]), fitFullModc$MZ00)
  cmod_imm2[,,i] <- mxEval((c0 + cMod0%x%cito_sc[i]), fitFullModc$MZ00)
  emod_imm2[,,i] <- mxEval((e0 + eMod0%x%cito_sc[i]), fitFullModc$MZ00)
  
  Va_del2[,,i]   <- mxEval((a1 + aMod1%x%cito_sc[i]) %*% t(a1 + aMod1%x%cito_sc[i]), fitFullModc$MZ11c)
  Vc_del2[,,i]   <- mxEval((c1 + cMod1%x%cito_sc[i]) %*% t(c1 + cMod1%x%cito_sc[i]), fitFullModc$MZ11c)
  Ve_del2[,,i]   <- mxEval((e1 + eMod1%x%cito_sc[i]) %*% t(e1 + eMod1%x%cito_sc[i]), fitFullModc$MZ11c)
  Vt_del2[,,i]   <- Va_del2[,,i] + Vc_del2[,,i] + Ve_del2[,,i]
  SVa_del2[,,i]  <- Va_del2[,,i]/Vt_del2[,,i]
  SVc_del2[,,i]  <- Vc_del2[,,i]/Vt_del2[,,i]
  SVe_del2[,,i]  <- Ve_del2[,,i]/Vt_del2[,,i]
  corA_del2[,,i] <- solve(sqrt(diag(2)*Va_del2[,,i]))%*%Va_del2[,,i]%*%solve(sqrt(diag(2)*Va_del2[,,i]))
  corC_del2[,,i] <- solve(sqrt(diag(2)*Vc_del2[,,i]))%*%Vc_del2[,,i]%*%solve(sqrt(diag(2)*Vc_del2[,,i]))
  corE_del2[,,i] <- solve(sqrt(diag(2)*Ve_del2[,,i]))%*%Ve_del2[,,i]%*%solve(sqrt(diag(2)*Ve_del2[,,i]))
  corP_del2[,,i] <- solve(sqrt(diag(2)*Vt_del2[,,i]))%*%Vt_del2[,,i]%*%solve(sqrt(diag(2)*Vt_del2[,,i]))
  amod_del2[,,i] <- mxEval((a1 + aMod1%x%cito_sc[i]), fitFullModc$MZ11c)
  cmod_del2[,,i] <- mxEval((c1 + cMod1%x%cito_sc[i]), fitFullModc$MZ11c)
  emod_del2[,,i] <- mxEval((e1 + eMod1%x%cito_sc[i]), fitFullModc$MZ11c)
}

outc_imm <- as.matrix(cbind(Va_imm2[2,2,],Vc_imm2[2,2,],Ve_imm2[2,2,],Vt_imm2[2,2,],SVa_imm2[2,2,],SVc_imm2[2,2,],SVe_imm2[2,2,], corA_imm2[2,1,],corC_imm2[2,1,],corE_imm2[2,1,], corP_imm2[2,1,], 
                            amod_imm2[2,1,], amod_imm2[2,2,], cmod_imm2[2,1,], cmod_imm2[2,2,], emod_imm2[2,1,], emod_imm2[2,2,]))
names(outc_imm) = c('Va_imm', 'Vc_imm', 'Ve_imm', 'Vt_imm', 'SVa_imm', 'SVc_imm', 'SVe_imm', 'corA_imm', 'corC_imm', 'corE_imm', 'corP_imm', 'amodC_imm', 'amodU_imm', 'cmodC_imm', 'cmodU_imm', 'emodC_imm', 'emodU_imm')
head(outc_imm)

outc_del <- as.matrix(cbind(Va_del2[2,2,],Vc_del2[2,2,],Ve_del2[2,2,],Vt_del2[2,2,],SVa_del2[2,2,],SVc_del2[2,2,],SVe_del2[2,2,], corA_del2[2,1,],corC_del2[2,1,],corE_del2[2,1,], corP_del2[2,1,],
                            amod_del2[2,1,], amod_del2[2,2,], cmod_del2[2,1,], cmod_del2[2,2,], emod_del2[2,1,], emod_del2[2,2,]))
names(outc_del) = c('Va_del', 'Vc_del', 'Ve_del', 'Vt_del', 'SVa_del', 'SVc_del', 'SVe_del', 'corA_del', 'corC_del', 'corE_del', 'corP_del', 'amodC_del', 'amodU_del', 'cmodC_del', 'cmodU_del', 'emodC_del', 'emodU_del')
head(outc_del)

outc_A <- cbind(outc_imm[,1], outc_del[,1], outc_imm[,5], outc_del[,5])
outc_C <- cbind(outc_imm[,2], outc_del[,2], outc_imm[,6], outc_del[,6])
outc_E <- cbind(outc_imm[,3], outc_del[,3], outc_imm[,7], outc_del[,7])


### Actual graphs

## Immediate vs Delayed: A, C, E separately (Supplementary Figure 1)
pdf("TablesFigures/ImmediateVsDelayed_A_C_E_04.pdf", width = 11.06, height = 7.76)
par(mar=c(6.1,5.1,4.1,3.1))
par(mfrow=c(2,3))

# A unstandardized
matplot(cito_sc, out_A[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Genes, Unstandardized")

# C unstandardized
matplot(cito_sc, out_C[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Shared environment, Unstandardized")

# E unstandardized
matplot(cito_sc, out_E[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Non-shared environment, Unstandardized")

par(mar=c(9.5,5.1,4.1,3.1))

# A standardized
matplot(cito_sc, out_A[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Genes, Standardized")
legend("bottom",  inset = c(0,-0.63), c("A immediate","A delayed"), lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# C standardized
matplot(cito_sc, out_C[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Shared environment, Standardized")
legend("bottom",  inset = c(0,-0.63), c("C immediate","C delayed"), lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# E standardized
matplot(cito_sc, out_E[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Non-shared environment, Standardized")
legend("bottom",  inset = c(0,-0.63), c("E immediate","E delayed"), lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)
dev.off()

## Immediate vs Delayed common: A, C, E separately (not shown in paper)
pdf("TablesFigures/ImmediateVsDelayedC_A_C_E_04.pdf", width = 11.06, height = 7.76)
par(mar=c(6.1,5.1,4.1,3.1))
par(mfrow=c(2,3))

# A unstandardized
matplot(cito_sc, out_AC[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Genes, Unstandardized")

# C unstandardized
matplot(cito_sc, out_CC[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Shared environment, Unstandardized")

# E unstandardized
matplot(cito_sc, out_EC[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Non-shared environment, Unstandardized")

par(mar=c(9.5,5.1,4.1,3.1))

# A standardized
matplot(cito_sc, out_AC[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Genes, Standardized")
legend("bottom",  inset = c(0,-0.63), c("A immediate","A delayed"), lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# C standardized
matplot(cito_sc, out_CC[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Shared environment, Standardized")
legend("bottom",  inset = c(0,-0.63), c("C immediate","C delayed"), lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# E standardized
matplot(cito_sc, out_EC[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Non-shared environment, Standardized")
legend("bottom",  inset = c(0,-0.63), c("E immediate","E delayed"), lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)
dev.off()


## Immediate vs Delayed unique: A, C, E separately (Figure 2 paper)
pdf("TablesFigures/ImmediateVsDelayedU_A_C_E_04.pdf", width = 11.06, height = 7.76)
par(mar=c(6.1,5.1,4.1,3.1))
par(mfrow=c(2,3))

# A unstandardized
matplot(cito_sc, out_AU[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Genes, Unstandardized")

# C unstandardized
matplot(cito_sc, out_CU[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Shared environment, Unstandardized")

# E unstandardized
matplot(cito_sc, out_EU[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Non-shared environment, Unstandardized")

par(mar=c(9.5,5.1,4.1,3.1))
# A standardized
matplot(cito_sc, out_AU[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Genes, Standardized")
legend("bottom",  inset = c(0,-0.63), c("A immediate","A delayed"), lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# C standardized
matplot(cito_sc, out_CU[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Shared environment, Standardized")
legend("bottom",  inset = c(0,-0.63), c("C immediate","C delayed"), lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# E standardized
matplot(cito_sc, out_EU[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Non-shared environment, Standardized")
legend("bottom",  inset = c(0,-0.63), c("E immediate","E delayed"), lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)
dev.off()


### Checks related to the note to Supplementary Figure 1

##  Check whether moderated paths remain above zero
pdf("TablesFigures/CheckModerationAbove0_04.pdf", width = 11.06, height = 7.76)
par(mar=c(6.1,5.1,4.1,3.1))
par(mfrow=c(2,3))

# a: immediate
matplot(cito_sc, out_imm[,12:13], type="l", lty=c(4,1), lwd=c(2,2), col=c(4,4), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="a paths immediate")

# c: immediate
matplot(cito_sc, out_imm[,14:15], type="l", lty=c(4,1), lwd=c(2,2), col=c(3,3), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="c paths immediate")

# e: immediate
matplot(cito_sc, out_imm[,16:17], type="l", lty=c(4,1), lwd=c(2,2), col=c(2,2), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="e paths immediate")

par(mar=c(9.5,5.1,4.1,3.1))
# a: delayed
matplot(cito_sc, out_del[,12:13], type="l", lty=c(4,1), lwd=c(2,2), col=c(4,4), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="a paths delayed")
legend("bottom",  inset = c(0,-0.63), c("Common","Unique"), lty=c(4,1), lwd=c(2,2), col=c(4,4), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# c: delayed
matplot(cito_sc, out_del[,14:15], type="l", lty=c(4,1), lwd=c(2,2), col=c(3,3), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="c paths delayed")
legend("bottom",  inset = c(0,-0.63), c("Common","Unique"), lty=c(4,1), lwd=c(2,2), col=c(3,3), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# e: delayed
matplot(cito_sc, out_del[,16:17], type="l", lty=c(4,1), lwd=c(2,2), col=c(2,2), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="e paths delayed")
legend("bottom",  inset = c(0,-0.63), c("Common","Unique"), lty=c(4,1), lwd=c(2,2), col=c(2,2), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)
dev.off()

##  Check whether moderated paths remain above zero after constraining c to 0
pdf("TablesFigures/CheckModerationAbove0con_04.pdf", width = 11.06, height = 7.76)
par(mar=c(6.1,5.1,4.1,3.1))
par(mfrow=c(2,3))

# a: immediate
matplot(cito_sc, outc_imm[,12:13], type="l", lty=c(4,1), lwd=c(2,2), col=c(4,4), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="a paths immediate")

# c: immediate
matplot(cito_sc, outc_imm[,14:15], type="l", lty=c(4,1), lwd=c(2,2), col=c(3,3), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="c paths immediate")

# e: immediate
matplot(cito_sc, outc_imm[,16:17], type="l", lty=c(4,1), lwd=c(2,2), col=c(2,2), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="e paths immediate")

par(mar=c(9.5,5.1,4.1,3.1))
# a: delayed
matplot(cito_sc, outc_del[,12:13], type="l", lty=c(4,1), lwd=c(2,2), col=c(4,4), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="a paths delayed")
legend("bottom",  inset = c(0,-0.63), c("Common","Unique"), lty=c(4,1), lwd=c(2,2), col=c(4,4), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# c: delayed
matplot(cito_sc, outc_del[,14:15], type="l", lty=c(4,1), lwd=c(2,2), col=c(3,3), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="c paths delayed")
legend("bottom",  inset = c(0,-0.63), c("Common","Unique"), lty=c(4,1), lwd=c(2,2), col=c(3,3), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# e: delayed
matplot(cito_sc, outc_del[,16:17], type="l", lty=c(4,1), lwd=c(2,2), col=c(2,2), xlab="Educational Performance (CITO score)",
        ylim = c(-0.5,1.5), ylab="Path effect", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="e paths delayed")
legend("bottom",  inset = c(0,-0.63), c("Common","Unique"), lty=c(4,1), lwd=c(2,2), col=c(2,2), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)
dev.off()

## Immediate vs Delayed: A, C, E separately if c constrained to 0
pdf("TablesFigures/ImmediateVsDelayed_A_C_Econ_04.pdf", width = 11.06, height = 7.76)
par(mar=c(6.1,5.1,4.1,3.1))
par(mfrow=c(2,3))

# A unstandardized
matplot(cito_sc, outc_A[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Genes, Unstandardized")

# C unstandardized
matplot(cito_sc, outc_C[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Shared environment, Unstandardized")

# E unstandardized
matplot(cito_sc, outc_E[,1:2], type="l", lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Non-shared environment, Unstandardized")

par(mar=c(9.5,5.1,4.1,3.1))

# A standardized
matplot(cito_sc, outc_A[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Genes, Standardized")
legend("bottom",  inset = c(0,-0.63), c("A immediate","A delayed"), lty=c(4,1), lwd=c(2,2), col=c("#0072B2","#0072B2"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# C standardized
matplot(cito_sc, outc_C[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Shared environment, Standardized")
legend("bottom",  inset = c(0,-0.63), c("C immediate","C delayed"), lty=c(4,1), lwd=c(2,2), col=c("#009E73","#009E73"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# E standardized
matplot(cito_sc, outc_E[,3:4], type="l", lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.4, cex.axis=1.2, cex.main=1.4, main="Non-shared environment, Standardized")
legend("bottom",  inset = c(0,-0.63), c("E immediate","E delayed"), lty=c(4,1), lwd=c(2,2), col=c("#D55E00","#D55E00"), cex = 1.4, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)
dev.off()


### Extra figures not shown in paper

## ACE in one: All tracking together
pdf("TablesFigures/AllTracking_ACE_04.pdf", width = 7.09, height = 9.43)
par(mar=c(9,5.1,4.1,2.1))
par(mfrow=c(1,2))
matplot(cito_sc, out[,1:4], type="l", lty=4:1, col=4:1, xlab="Educational Performance (CITO score)",
        ylim = c(0,2), ylab="Variance", main="Total - Unstandardized\nvariance components")
legend("bottom",  inset = c(0,-0.72), c("Va","Vc","Ve","Vt"), lty=4:1, col=4:1, cex = 1, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 1)

matplot(cito_sc, out[,5:7], type="l", lty=4:2, col=4:2, xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", main="Standardized\nvariance components")
legend("bottom",  inset = c(0,-0.65), c("A","C","E"), lty=4:1, col=4:1, cex = 1, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)
dev.off()

## ACE in one: Immediate vs Delayed 
pdf("TablesFigures/ImmediateVsDelayed_ACE_04.pdf", width = 7.09, height = 9.43)
par(mar=c(11,5.1,4.1,2.1))
par(mfrow=c(2,2))

# Unstandardized
# Immediate tracking 
matplot(cito_sc, out_imm[,1:4], type="l", lty=4:1, lwd=c(2,2,2,2), col=4:1, xlab="Educational Performance \n(CITO score centered)",
        ylim = c(0,2.5), ylab="Variance", cex.lab=1.1, main="Immediate tracking\nUnstandardized")
legend("bottom",  inset = c(0,-0.77), c("Genetic","Shared Environmental","Non-shared Environmental","Total"), lty=4:1, lwd=c(2,2,2,2), col=4:1, cex = 1.1, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 1)

# Delayed tracking
matplot(cito_sc, out_del[,1:4], type="l", lty=4:1, lwd=c(2,2,2,2), col=4:1, xlab="Educational Performance (CITO score centered)",
        ylim = c(0,2.5), ylab="Variance", cex.lab=1.1, main="Delayed tracking\nUnstandardized")
legend("bottom",  inset = c(0,-0.77), c("Genetic","Shared Environmental","Non-shared Environmental","Total"), lty=4:1, lwd=c(2,2,2,2), col=4:1, cex = 1.1, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 1)

par(mar=c(10,5.1,5.1,2.1))

# Standardized
# Immediate tracking
matplot(cito_sc, out_imm[,5:7], type="l", lty=4:2, lwd=c(2,2,2), col=4:2, xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.1, main="Immediate tracking\nStandardized")
legend("bottom",  inset = c(0,-0.67), c("Genetic (A)","Shared Environmental (C)","Non-shared Environmental (E)"), lty=4:1, lwd=c(2,2,2), col=4:1, cex = 1.1, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)

# Delayed tracking
matplot(cito_sc, out_del[,5:7], type="l", lty=4:2, lwd=c(2,2,2), col=4:2, xlab="Educational Performance (CITO score)",
        ylim = c(0,1), ylab="Proportion of Variance", cex.lab=1.1, main="Delayed tracking\nStandardized")
legend("bottom",  inset = c(0,-0.67), c("Genetic (A)","Shared Environmental (C)","Non-shared Environmental (E)"), lty=4:1, lwd=c(2,2,2), col=4:1, cex = 1.1, y.intersp = 1, ncol=1, bty = "n", xpd=TRUE, yjust = 0)
dev.off()


## Correlation between performance and attainment
pdf("TablesFigures/CorrPerfAttain_04.pdf", width = 11.06, height = 7.76)
par(mfrow=c(1,2))
matplot(cito_sc, out_imm[,8:11], type='l', lty=4:1, col=4:1, xlab='Cito score',
        ylab="Correlation Performance and Attainment", main="Immediate")
legend("bottomright", c("rA","rC","rE","rP"), lty=4:1, col=4:1, cex = 0.8, y.intersp = 1)

matplot(cito_sc, out_del[,8:11], type='l', lty=4:1, col=4:1, xlab='Cito score',
        ylab="Correlation Performance and Attainment", main="Delayed")
legend("bottomright", c("rA","rC","rE","rP"), lty=4:1, col=4:1, cex = 0.8, y.intersp = 1)
dev.off()

# If c constrained to 0
pdf("TablesFigures/CorrPerfAttain_con_04.pdf", width = 11.06, height = 7.76)
par(mfrow=c(1,2))
matplot(cito_sc, outc_imm[,8:11], type='l', lty=4:1, col=4:1, xlab='Cito score',
        ylim = c(0,1), ylab="Correlation Performance and Attainment", main="Immediate")
legend("bottomright", c("rA","rC","rE","rP"), lty=4:1, col=4:1, cex = 0.8, y.intersp = 1)

matplot(cito_sc, outc_del[,8:11], type='l', lty=4:1, col=4:1, xlab='Cito score',
        ylim = c(0,1), ylab="Correlation Performance and Attainment", main="Delayed")
legend("bottomright", c("rA","rC","rE","rP"), lty=4:1, col=4:1, cex = 0.8, y.intersp = 1)
dev.off()


# ----------------------------------------------------------------------------------------------------------------------
sink()
save.image(paste(filename,".Rdata",sep=""))
