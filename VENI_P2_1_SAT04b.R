# ----------------------------------------------------------------------------------------------------------------------
# Program: oneSATc.R  
#  Author: Hermine Maes (2018-01-04)
# Adapted: Antonie Knigge  
#    Date: 2018-11-27 
#    vs03: 2020-10-09, Use newer version of Openmx (same as main analyses)
#    vs04: 2021-02-03, Used dataset version 6 with also pairs where one twin has info on cito included
#
# Twin Univariate Saturated model to estimate means and (co)variances across multiple groups
# Matrix style model - Raw data - Continuous data
# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|

# Load Libraries & Options
rm(list=ls())
library(OpenMx)
library(psych); library(polycor)
library(readstata13)
source("miFunctions.R")

# Set the right working directory
#setwd("C:/Users/knigg101/surfdrive/Onderzoek/VENI/P2/WorkInProgress/Analyses")
#setwd("/Users/Antonie/surfdrive/Onderzoek/VENI/P2/WorkInProgress/Analyses")

# Create Output 
filename    <- "VENI_P2_1_SAT04b"
sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

# Load Data

NtrRawData <- read.dta13("NTR_P2_1_06.dta", convert.dates=T, convert.underscore=F, convert.factors=F)
twinData   <- reshape(NtrRawData, idvar="twin_id", timevar="mult_ext", 
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
describe(twinData[,1:46], skew=F)

# Give pseudo-missing values for definition variables where only one twin was observed.
# (phenotypes already have missing values for those twin pairs; 
# pseudo-missings were already assigned in stata to twin pairs where both were observed but one has missing info)
twinData$male_1[is.na(twinData$male_1)] <- -999
twinData$defm_1[is.na(twinData$defm_1)] <- -999
twinData$male_2[is.na(twinData$male_2)] <- -999
twinData$defm_2[is.na(twinData$defm_2)] <- -999

# Select Variables for Analysis
vars            <- 'c_cito'                  # list of variables names
nv              <- length(vars)	             # number of variables
ntv             <- nv*2                      # number of total variables
selVars         <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="_")
covVars1        <- c('male_1', 'c_yob', 'male_2')
nc1             <- 2                         # number of covariates
twinData$defm_1 <- twinData$c_cito_1         # create a copy of moderator as a definition variable
twinData$defm_2 <- twinData$c_cito_2
defVars         <- c('defm_1', 'defm_2')
useVars	        <- c(selVars, covVars1, defVars, 'twzyg', 'delay_1', 'delay_2')

# Select Cases for Analysis
twinDataNM <- subset(twinData, (male_1!="-999" | male_2!="-999") & c_yob!="NA" & (c_cito_1!="NA" | c_cito_2!="NA" | edu_1!="NA" | edu_2!="NA"), c(useVars, 'FamilyNumber'))
dim(twinDataNM)
summary(twinDataNM)

# ----------------------------------------------------------------------------------------------------------------------
# SATURATED, multigroup for all possible zygosity types, so separately for sex
# ----------------------------------------------------------------------------------------------------------------------

# Select Data for Analysis
mzmData	<- subset(twinDataNM, twzyg==1, selVars)
dzmData	<- subset(twinDataNM, twzyg==2, selVars)
mzfData	<- subset(twinDataNM, twzyg==3, selVars)
dzfData	<- subset(twinDataNM, twzyg==4, selVars)
dmfData	<- subset(twinDataNM, twzyg==5, selVars) # DZ male-female
dfmData	<- subset(twinDataNM, twzyg==6, selVars) # DZ female-male

# Generate Descriptive Statistics 
colMeans(mzmData,na.rm=TRUE) 
colMeans(dzmData,na.rm=TRUE) 
colMeans(mzfData,na.rm=TRUE) 
colMeans(dzfData,na.rm=TRUE) 
colMeans(dmfData,na.rm=TRUE) 
colMeans(dfmData,na.rm=TRUE) 
cov(mzmData,use="complete") 
cov(dzmData,use="complete")
cov(mzfData,use="complete") 
cov(dzfData,use="complete")
cov(dmfData,use="complete")
cov(dfmData,use="complete")
cor(mzmData,use="complete") 
cor(dzmData,use="complete")
cor(mzfData,use="complete") 
cor(dzfData,use="complete")
cor(dmfData,use="complete")
cor(dfmData,use="complete")

# Set Starting Values
svMe      <-  0                        # start value for means
svVa      <-  69                       # start value for variance
lbVa      <- .0001                     # lower bound for variance


# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create Algebra for expected Mean Matrices
meanMZm    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZm1" ,"mMZm2" ), name="meanMZm"  )
meanDZm    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZm1" ,"mDZm2" ), name="meanDZm"  )
meanMZf    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZf1" ,"mMZf2" ), name="meanMZf"  )
meanDZf    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZf1" ,"mDZf2" ), name="meanDZf"  )
meanDZmf   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZmf1","mDZmf2"), name="meanDZmf" )
meanDZfm   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZfm1","mDZfm2"), name="meanDZfm" )

# Create Algebra for expected Variance/Covariance Matrices
covMZm     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vMZm1" ,"cMZm21" ,"vMZm2" ), name="covMZm"  )
covDZm     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vDZm1" ,"cDZm21" ,"vDZm2" ), name="covDZm"  )
covMZf     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vMZf1" ,"cMZf21" ,"vMZf2" ), name="covMZf"  )
covDZf     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vDZf1" ,"cDZf21" ,"vDZf2" ), name="covDZf"  )
covDZmf    <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vDZmf1","cDZmf21","vDZmf2"), name="covDZmf" )
covDZfm    <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vDZfm1","cDZfm21","vDZfm2"), name="covDZfm" )

# Create Data Objects for Multiple Groups
dataMZm    <- mxData( observed=mzmData, type="raw" )
dataDZm    <- mxData( observed=dzmData, type="raw" )
dataMZf    <- mxData( observed=mzfData, type="raw" )
dataDZf    <- mxData( observed=dzfData, type="raw" )
dataDZmf   <- mxData( observed=dmfData, type="raw" )
dataDZfm   <- mxData( observed=dfmData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZm     <- mxExpectationNormal( covariance="covMZm" , means="meanMZm" , dimnames=selVars )
expDZm     <- mxExpectationNormal( covariance="covDZm" , means="meanDZm" , dimnames=selVars )
expMZf     <- mxExpectationNormal( covariance="covMZf" , means="meanMZf" , dimnames=selVars )
expDZf     <- mxExpectationNormal( covariance="covDZf" , means="meanDZf" , dimnames=selVars )
expDZmf    <- mxExpectationNormal( covariance="covDZmf", means="meanDZmf", dimnames=selVars )
expDZfm    <- mxExpectationNormal( covariance="covDZfm", means="meanDZfm", dimnames=selVars )
funML      <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
modelMZm   <- mxModel( meanMZm , covMZm , dataMZm , expMZm , funML, name="MZm"  )
modelDZm   <- mxModel( meanDZm , covDZm , dataDZm , expDZm , funML, name="DZm"  )
modelMZf   <- mxModel( meanMZf , covMZf , dataMZf , expMZf , funML, name="MZf"  )
modelDZf   <- mxModel( meanDZf , covDZf , dataDZf , expDZf , funML, name="DZf"  )
modelDZmf  <- mxModel( meanDZmf, covDZmf, dataDZmf, expDZmf, funML, name="DZmf" )
modelDZfm  <- mxModel( meanDZfm, covDZfm, dataDZfm, expDZfm, funML, name="DZfm" )
multi      <- mxFitFunctionMultigroup( c("MZm","DZm", "MZf", "DZf", "DZmf", "DZfm") )

# Create Confidence Interval Objects
ciCov     <- mxCI( c('MZm.covMZm','DZm.covDZm','MZf.covMZf','DZf.covDZf','DZmf.covDZmf','DZfm.covDZfm') )
ciMean    <- mxCI( c('MZm.meanMZm','DZm.meanDZm','MZf.meanMZf','DZf.meanDZf','DZmf.meanDZmf','DZfm.meanDZfm') )

# Build Saturated Model with Confidence Intervals
modelSAT  <- mxModel( "oneSATc", modelMZm, modelDZm, modelMZf, modelDZf, modelDZmf, modelDZfm, multi, ciCov, ciMean )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run Saturated Model
fitSAT    <- mxTryHard( modelSAT, intervals=F )
#summary(fitSAT)

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitSAT)
fitEsts(fitSAT, 2)
#mxGetExpected( fitSAT, c("means","covariance") )

# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# Constrain Variances to be equal across Twin Order
modelEVO <- mxModel( fitSAT, name="oneEVOc" )
modelEVO <- omxSetParameters( modelEVO, label=c("vMZm1" ,"vMZm2" ), free=TRUE, values=svVa, newlabels='vMZm'  )
modelEVO <- omxSetParameters( modelEVO, label=c("vDZm1" ,"vDZm2" ), free=TRUE, values=svVa, newlabels='vDZm'  )
modelEVO <- omxSetParameters( modelEVO, label=c("vMZf1" ,"vMZf2" ), free=TRUE, values=svVa, newlabels='vMZf'  )
modelEVO <- omxSetParameters( modelEVO, label=c("vDZf1" ,"vDZf2" ), free=TRUE, values=svVa, newlabels='vDZf'  )
modelEVO <- omxSetParameters( modelEVO, label=c("vDZmf1","vDZmf2"), free=TRUE, values=svVa, newlabels='vDZmf' )
modelEVO <- omxSetParameters( modelEVO, label=c("vDZfm1","vDZfm2"), free=TRUE, values=svVa, newlabels='vDZfm' )
fitEVO   <- mxRun( modelEVO, intervals=F )
fitGofs(fitEVO); fitEsts(fitEVO, 2)

# Constrain Variances to be equal across Twin Order and Zygosity
modelEVOZ <- mxModel( fitEVO, name="oneEVOZc" )
modelEVOZ <- omxSetParameters( modelEVOZ, label=c("vMZm" ,"vDZm" ), free=TRUE, values=svVa, newlabels='vZm' )
modelEVOZ <- omxSetParameters( modelEVOZ, label=c("vMZf" ,"vDZf" ), free=TRUE, values=svVa, newlabels='vZf' )
modelEVOZ <- omxSetParameters( modelEVOZ, label=c("vDZmf","vDZfm"), free=TRUE, values=svVa, newlabels='vZo' )
fitEVOZ   <- mxRun( modelEVOZ, intervals=F )
fitGofs(fitEVOZ); fitEsts(fitEVOZ, 2)

# Constrain Variances to be equal across Twin Order, Zygosity, and Sex
modelEVOZS <- mxModel( fitEVOZ, name="oneEVOZSc" )
modelEVOZS <- omxSetParameters( modelEVOZS, label=c("vZm","vZf", "vZo"), free=TRUE, values=svVa, newlabels='vZ' )
fitEVOZS   <- mxRun( modelEVOZS, intervals=F )
fitGofs(fitEVOZS); fitEsts(fitEVOZS, 2)

# Print Comparative Fit Statistics
mxCompare( fitSAT    , subs <- list(fitEVO, fitEVOZ, fitEVOZS) )
mxCompare( fitEVO    , fitEVOZ )
mxCompare( fitEVOZ   , fitEVOZS )


## Also constrain means in between

# Constrain Variances and Means to be equal across Twin Order
modelEVMO  <- mxModel( fitEVO, name="oneEVMOc" )
modelEVMO  <- omxSetParameters( modelEVMO, label=c("mMZm1" ,"mMZm2" ), free=TRUE, values=svMe, newlabels='mMZm'  )
modelEVMO  <- omxSetParameters( modelEVMO, label=c("mDZm1" ,"mDZm2" ), free=TRUE, values=svMe, newlabels='mDZm'  )
modelEVMO  <- omxSetParameters( modelEVMO, label=c("mMZf1" ,"mMZf2" ), free=TRUE, values=svMe, newlabels='mMZf'  )
modelEVMO  <- omxSetParameters( modelEVMO, label=c("mDZf1" ,"mDZf2" ), free=TRUE, values=svMe, newlabels='mDZf'  )
modelEVMO  <- omxSetParameters( modelEVMO, label=c("mDZmf1","mDZmf2"), free=TRUE, values=svMe, newlabels='mDZmf' )
modelEVMO  <- omxSetParameters( modelEVMO, label=c("mDZfm1","mDZfm2"), free=TRUE, values=svMe, newlabels='mDZfm' )
fitEVMO    <- mxRun( modelEVMO, intervals=F )
fitGofs(fitEVMO); fitEsts(fitEVMO, 2)

# Constrain Variances and Means to be equal across Twin Order and Variances across Zygosity
modelEVZ <- mxModel( fitEVMO, name="oneEVZc" )
modelEVZ <- omxSetParameters( modelEVZ, label=c("vMZm" ,"vDZm" ), free=TRUE, values=svVa, newlabels='vZm' )
modelEVZ <- omxSetParameters( modelEVZ, label=c("vMZf" ,"vDZf" ), free=TRUE, values=svVa, newlabels='vZf' )
modelEVZ <- omxSetParameters( modelEVZ, label=c("vDZmf","vDZfm"), free=TRUE, values=svVa, newlabels='vZo' )
fitEVZ   <- mxRun( modelEVZ, intervals=F )
fitGofs(fitEVZ); fitEsts(fitEVZ, 2)

# Constrain Variances and Means to be equal across Twin Order and Zygosity
modelEVMZ <- mxModel( fitEVZ, name="oneEVMZc" )
modelEVMZ <- omxSetParameters( modelEVMZ, label=c("mMZm" ,"mDZm" ), free=TRUE, values=svMe, newlabels='mZm' )
modelEVMZ <- omxSetParameters( modelEVMZ, label=c("mMZf" ,"mDZf" ), free=TRUE, values=svMe, newlabels='mZf' )
modelEVMZ <- omxSetParameters( modelEVMZ, label=c("mDZmf","mDZfm"), free=TRUE, values=svMe, newlabels='mZo' )
fitEVMZ   <- mxRun( modelEVMZ, intervals=F )
fitGofs(fitEVMZ); fitEsts(fitEVMZ, 2)

# Constrain Variances and Means to be equal across Twin Order and Variances across Sex
modelEVS <- mxModel( fitEVMO, name="oneEVSc" )
modelEVS <- omxSetParameters( modelEVS, label=c("vMZm","vMZf"), free=TRUE, values=svVa, newlabels='vMZ' )
modelEVS <- omxSetParameters( modelEVS, label=c("vDZm","vDZf","vDZmf","vDZfm"), free=TRUE, values=svVa, newlabels='vDZ' )
fitEVS   <- mxRun( modelEVS, intervals=F )
fitGofs(fitEVS); fitEsts(fitEVS, 2)

# Constrain Variances and Means to be equal across Twin Order and Sex
modelEVMS <- mxModel( fitEVS, name="oneEVMSc" )
modelEVMS <- omxSetParameters( modelEVMS, label=c("mMZm","mMZf"), free=TRUE, values=svMe, newlabels='mMZ' )
modelEVMS <- omxSetParameters( modelEVMS, label=c("mDZm","mDZf","mDZmf","mDZfm"), free=TRUE, values=svMe, newlabels='mDZ' )
fitEVMS   <- mxRun( modelEVMS, intervals=F )
fitGofs(fitEVMS); fitEsts(fitEVMS, 2)

# Constrain Variances and Means to be equal across Twin Order and Zygosity and Variances across Sex
modelEVMZ_VS <- mxModel( fitEVMZ, name="oneEVMZ_VSc" )
modelEVMZ_VS <- omxSetParameters( modelEVMZ_VS, label=c("vZm","vZf", "vZo"), free=TRUE, values=svVa, newlabels='vZ' )
fitEVMZ_VS   <- mxRun( modelEVMZ_VS, intervals=F )
fitGofs(fitEVMZ_VS); fitEsts(fitEVMZ_VS, 2)

# Constrain Variances and Means to be equal across Twin Order, Zygosity, and Sex
modelEVMZ_VMS <- mxModel( fitEVMZ_VS, name="oneEVMZ_VMSc" )
modelEVMZ_VMS <- omxSetParameters( modelEVMZ_VMS, label=c("mZm","mZf", "mZo"), free=TRUE, values=svMe, newlabels='mZ' )
fitEVMZ_VMS   <- mxTryHard( modelEVMZ_VMS, intervals=F )
fitGofs(fitEVMZ_VMS); fitEsts(fitEVMZ_VMS, 2)


# Print Comparative Fit Statistics
mxCompare( fitSAT    , subs <- list(fitEVO, fitEVMO, fitEVZ, fitEVMZ, fitEVS, fitEVMS, fitEVMZ_VS, fitEVMZ_VMS) )
mxCompare( fitEVO    , fitEVMO )
mxCompare( fitEVMO   , fitEVZ )
mxCompare( fitEVZ    , fitEVMZ )
mxCompare( fitEVMZ   , fitEVMZ_VS )
mxCompare( fitEVMZ_VS, fitEVMZ_VMS )

# Extra check
mxCompare( fitEVMO   , fitEVS )
mxCompare( fitEVS    , fitEVMS )


# ----------------------------------------------------------------------------------------------------------------------
# SATURATED, multigroup for MZ and DZ zygosity types with covariate for sex
# ----------------------------------------------------------------------------------------------------------------------

# Select Covariates for Analysis
twinDataSub  <- twinDataNM
dim(twinDataSub)

# Select Data for Analysis
mzData	<- subset(twinDataSub, twzyg==1 | twzyg==3                      , useVars)
dzData	<- subset(twinDataSub, twzyg==2 | twzyg==4 | twzyg==5 | twzyg==6, useVars)

# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")

# Set Starting Values
svBe      <- 0.01                      # start value for regressions
svMe      <- -0.5                      # start value for means
svVa      <- 69                        # start value for variance
lbVa      <- 0.0001                    # lower bound for variance

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create Matrices for Covariates and linear Regression Coefficients
defL      <- mxMatrix( type="Full", nrow=1, ncol=2, free=FALSE, labels=c("data.male_1", "data.male_2"), name="defL" )
pathBl    <- mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=svBe, labels =c("b11","b11"), name="b1" )

# Create Algebra for expected Mean Matrices
meanMZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZ1","mMZ2"), name="meanMZ" )
meanDZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZ1","mDZ2"), name="meanDZ" )
expMeanMZ <- mxAlgebra( expression= meanMZ + b1*defL, name="expMeanMZ" )
expMeanDZ <- mxAlgebra( expression= meanDZ + b1*defL, name="expMeanDZ" )

# Create Algebra for expected Variance/Covariance Matrices
covMZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vMZ1","cMZ21","vMZ2"), name="covMZ" )
covDZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv), labels=c("vDZ1","cDZ21","vDZ2"), name="covDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="covMZ", means="expMeanMZ", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="covDZ", means="expMeanDZ", dimnames=selVars )
funML_c   <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars       <- list( pathBl )
defs       <- list( defL )
modelMZ    <- mxModel( pars, defs, meanMZ, expMeanMZ, covMZ, dataMZ, expMZ, funML_c, name="MZ" )
modelDZ    <- mxModel( pars, defs, meanDZ, expMeanDZ, covDZ, dataDZ, expDZ, funML_c, name="DZ" )
multi_c    <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Confidence Interval Objects
ciCov_c     <- mxCI( c('MZ.covMZ','DZ.covDZ') )
ciMean_c    <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )

# Build Saturated Model with Confidence Intervals
modelSAT_c  <- mxModel( "oneSATca", pars, modelMZ, modelDZ, multi_c, ciCov_c, ciMean_c )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run Saturated Model
fitSAT_c    <- mxTryHard( modelSAT_c, intervals=F )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitSAT_c)
fitEsts(fitSAT_c, 2)
#mxGetExpected( fitSAT_c, c("means","covariance") )

# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# Test Significance of Covariate
modelCOV  <- mxModel( fitSAT_c, name="oneCOVcc" )
modelCOV  <- omxSetParameters( modelCOV, label="b11", free=FALSE, values=0 )
fitCOV    <- mxTryHard( modelCOV )

# Constrain Variances to be equal across Twin Order
modelEVO_c <- mxModel( fitSAT_c, name="oneEVOcc" )
modelEVO_c <- omxSetParameters( modelEVO_c, label=c("vMZ1","vMZ2"), free=TRUE, values=svVa, newlabels='vMZ' )
modelEVO_c <- omxSetParameters( modelEVO_c, label=c("vDZ1","vDZ2"), free=TRUE, values=svVa, newlabels='vDZ' )
fitEVO_c   <- mxTryHard( modelEVO_c, intervals=F )
fitGofs(fitEVO_c); fitEsts(fitEVO_c, 2)

# Constrain Variances and Means to be equal across Twin Order
modelEVMO_c  <- mxModel( fitEVO_c, name="oneEVMOcc" )
modelEVMO_c  <- omxSetParameters( modelEVMO_c, label=c("mMZ1","mMZ2"), free=TRUE, values=svMe, newlabels='mMZ' )
modelEVMO_c  <- omxSetParameters( modelEVMO_c, label=c("mDZ1","mDZ2"), free=TRUE, values=svMe, newlabels='mDZ' )
fitEVMO_c    <- mxTryHard( modelEVMO_c, intervals=F )
fitGofs(fitEVMO_c); fitEsts(fitEVMO_c, 2)

# Constrain Variances and Means to be equal across Twin Order and Variances across Zygosity
modelEVZ_c <- mxModel( fitEVMO_c, name="oneEVZcc" )
modelEVZ_c <- omxSetParameters( modelEVZ_c, label=c("vMZ","vDZ"), free=TRUE, values=svVa, newlabels='vZ' )
fitEVZ_c   <- mxTryHard( modelEVZ_c, intervals=F )
fitGofs(fitEVZ_c); fitEsts(fitEVZ_c, 2)

# Constrain Variances and Means to be equal across Twin Order and Zygosity
modelEVMZ_c <- mxModel( fitEVZ_c, name="oneEVMZcc" )
modelEVMZ_c <- omxSetParameters( modelEVMZ_c, label=c("mMZ","mDZ"), free=TRUE, values=svMe, newlabels='mZ' )
fitEVMZ_c   <- mxTryHard( modelEVMZ_c, intervals=F )
fitGofs(fitEVMZ_c); fitEsts(fitEVMZ_c, 2)


# Print Comparative Fit Statistics
mxCompare( fitSAT_c , subs <- list(fitCOV, fitEVO_c, fitEVMO_c, fitEVZ_c, fitEVMZ_c) )
mxCompare( fitEVO_c , fitEVMO_c )
mxCompare( fitEVMO_c, fitEVZ_c )
mxCompare( fitEVZ_c , fitEVMZ_c )


# ----------------------------------------------------------------------------------------------------------------------
sink()
save.image(paste(filename,".Rdata",sep=""))
