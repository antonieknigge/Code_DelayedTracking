/*	
	In order to run this file, one needs to set the working directory correctly in the Stata command window.
	Typically, when still working on the file, this will be something like:
	cd "C:\Users\knigg101\surfdrive\Onderzoek\VENI"
	cd "/Users/Antonie/surfdrive/Onderzoek/VENI"
*/

capture log close
log using "P2/WorkInProgress/Analyses/VENI_P2_DescriptiveAnalyses05", replace text 

/*
File: 			VENI_P2_DescriptiveAnalyses05.do
Author(s): 		a.knigge@uu.nl 
Date: 			2020-01-14
Purpose: 		Clean dataset and create variables.
Source data:	NTR_P2_1_06.dta
Contains: 		Data from the Netherlands Twin Register prepared for analyses
Remarks:		Comments starting and ending with !! point out that action may be required if you change anything to this do-file
*/

//	STEP 0 ==========================================================================================================
//	Program Setup

version 13.1
clear all
set linesize 140
macro drop _all
set more off


*	We make locals for the path(s) of the folder(s) in which the data that we use and save can be found. 
*	The purpose? If there is a change in the directory structure, the paths only need to be updated once, namely here.
local path1 = "P2/WorkInProgress/Analyses/"

/*	
	The dataset(s) that is/are used may change. Then, they get a different version number. 
	!!Change numbers below to use the newest version of a dataset!! 
*/
local VersionUF1 = "06"		// UF1 = UsedFile1 

use "`path1'NTR_P2_1_`VersionUF1'", clear

/*
	We save in this do-file some graphs  
	If anything changes in the do-file: 
	!!Update below the version number of the figures that are saved!!
*/
local VersionSF1 = "03"		// SF1 = SavedFigure1


//	STEP 1 ==========================================================================================================
//	Descriptives

*	Generate descriptives of the variables for table in paper
*	(cito_nm>0 indicates that at least one of the twin pair should have non-missing on cito)
sum edu cito_final delay yob male age 
sum edu cito_final delay yob male age if cito_nm>0 & c_cito<.			
sum edu cito_final delay yob male age if cito_nm>0 & c_cito<. & mz==1
sum edu cito_final delay yob male age if cito_nm>0 & c_cito<. & mz==0

*	Overview of number of pairs in the analyses
codebook twin_id 
codebook twin_id if cito_nm>0 & c_cito<.
codebook twin_id if cito_nm>0 & c_cito<. & mz==1
codebook twin_id if cito_nm>0 & c_cito<. & mz==0

*	Size of timing of tracking groups
codebook delay if cito_nm>0 & c_cito<.

*	Descriptive information on selectivity for selectivity section paper
tab oplm       if cito_nm>0 & c_cito<. & oplm>-1
tab oplv 	   if cito_nm>0 & c_cito<. & oplv>-1
tab allochtoon if cito_nm>0 & c_cito<. & allochtoon>-1


//	STEP 2 ==========================================================================================================
//	Relate choice for homogeneous vs heterogeneous class to parental education and own achievement 

*	Necessary to take dependence between twins into account
svyset twin_id

*	Relation with parental education
	 tab edup4c 	delay if cito_nm>0 & c_cito<., row chi2
svy: tab edup4c 	delay if cito_nm>0 & c_cito<., row			// Table 4
ttest eduparents 		  if cito_nm>0 & c_cito<., by(delay)

*	Relation with own achievement
	 tab achie5   	delay if cito_nm>0 & c_cito<., row chi2
svy: tab achie5   	delay if cito_nm>0 & c_cito<., row
	 tab achie6 	delay if cito_nm>0 & c_cito<., row chi2 
svy: tab achie6 	delay if cito_nm>0 & c_cito<., row 
	 tab achie8  	delay if cito_nm>0 & c_cito<., row chi2
svy: tab achie8  	delay if cito_nm>0 & c_cito<., row			// Table 5
	 tab cito_final delay if cito_nm>0 & c_cito<., row chi2
svy: tab cito_final delay if cito_nm>0 & c_cito<., row
ttest cito_final          if cito_nm>0 & c_cito<., by(delay)

corr edu cito_final achie8 achie6 achie5 if cito_nm>0 & c_cito<.

*	Relation with parental education given own achievement (not in paper)
bysort achie8: tab edup4c delay if cito_nm>0 & c_cito<., chi2 row 
table achie8 edup4c             if cito_nm>0 & c_cito<., contents(mean delay)


//	STEP 3 ==========================================================================================================
//	Get an idea how cito score relates to educational attainment (not in paper)

lowess edu12 	cito_final if cito_nm>0 & c_cito<. & mult_ext==1,  jitter(3) 
graph export "`path1'TablesFigures/VENI_P2_LowessEdu12_`VersionSF1'.pdf", replace 
lowess edu   	cito_final if cito_nm>0 & c_cito<. & mult_ext==1,  jitter(3)
graph export "`path1'TablesFigures/VENI_P2_LowessEdu_`VersionSF1'.pdf", replace 
lowess edu_alt  cito_final if cito_nm>0 & c_cito<. & mult_ext==1,  jitter(3)
graph export "`path1'TablesFigures/VENI_P2_LowessEdu_alt_`VersionSF1'.pdf", replace 


//	STEP FINALIZE ===================================================================================================
//	Metadata

*	!!If you make a change to the do-file, make a short entry of the change in the log below!!
/*
	Date		Author			Description of change
	2020-09-09	a.knigge@uu.nl	Used the new scoring of edu; used achie8 instead of achie6 in the tables;
								Added input for descriptive table of paper; Added input for selectivity numbers paper;
								Based all results on analytic sample as much as possible;
								Fixed path of log.
	2020-12-21	a.knigge@uu.nl	Fixed Chi2 test for dependence of observations	
	2021-02-02	a.knigge@uu.nl	Used dataset version 6: which includes those pairs where one twin has info on CITO.
	2021-12-09	a.knigge@uu.nl	Added info on size of immediate/delayed groups.
*/

log close
