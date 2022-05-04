/*	
	In order to run this file, one needs to set the working directory correctly in the Stata command window.
	Typically, when still working on the file, this will be something like:
	cd "C:\Users\knigg101\surfdrive\Onderzoek\VENI"
	cd "/Users/Antonie/surfdrive/Onderzoek/VENI"
*/

capture log close
log using "P2/WorkInProgress/DataCleaning/VENI_P2_1_Clean06", replace text 

/*
File: 			VENI_P2_1_Clean06.do
Author(s): 		a.knigge@uu.nl 
Date: 			2018-10-19
Purpose: 		Clean dataset and create variables.
Source data:	NTR_P2_0_01.dta
Contains: 		Source data as received from the Netherlands Twin Register 
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
local VersionUF1 = "01"		// UF1 = UsedFile1 

/*
	We save in this do-file: .dta  
	If anything changes in the do-file: 
	!!Update below the version number of the data that are saved!!
*/
local VersionSF1 = "06"		// SF1 = SavedFile1

use "`path1'NTR_P2_0_`VersionUF1'", clear


//	STEP 1 ==========================================================================================================
//	Create Outcome Variables


//	STEP 1.1 --------------------------------------------------------------------------------------------------------
//	Educational attainment age 12 filled in by the mother
 
** 	Explore which categories there are over the years

*  	Which survey is taken when?
tab invjrm12 if hschm12_v2<.
tab invjrm12 if hschm12_v3<.
tab invjrm12 if hschm12   <.
tab invjrm12 if hschm12_1n<.
tab invjrm12 if hschm12_2n<.
tab invjrm12 if hschm12_3n<.
tab invjrm12 if vervolgm12<.	// Same results when twzyg<. is added

*	What years of birth do each survey cover?
tab yob if hschm12_v2<.
tab yob if hschm12_v3<.
tab yob if hschm12   <.
tab yob if hschm12_1n<.
tab yob if hschm12_2n<.
tab yob if hschm12_3n<.
tab yob if vervolgm12<.			// Same results when twzyg<. is added

*   Which categories are there in each survey?
tab hschm12_v2
tab hschm12_v3
tab hschm12
tab hschm12_1n
tab hschm12_2n
tab hschm12_3n
tab vervolgm12

/*
Variabele	Invuljaar	Geboortejaar

hschm12_v2	1999-2000	1986-1987		1999		Lijst12_versie2_ouders_groenBoekje_1999			SVO/LBO/MBO/	 BRUGKLAS  /MAVO/HAVO/VWO-GYM
hschm12_v3	1999-2004	1986-1991		2000		Lijst12_versie3_ouders_paarsBoekje_2000			SVO/LBO/MBO/VMBO/BRUGKLASx3/MAVO/HAVO/VWO-GYM
hschm12		2003-2015	1991-2002		2003		Lijst12_versie4a_ouders_geelBoekje_2003 		SVO/PRO/VMBOx4/HAVO/VWO (meerdere) 
hschm12_1n	2009-2015	1995-2002		2006		Lijst12_versie5_ouders_geel-grijsBoekje_2006	SVO/PRO/VMBOx4/HAVO/VWO	
hschm12_2n	2009-2014	1997-2001		/2009		/Lijst12_versie6_ouders1Lijst_druk2_2009			" "
hschm12_3n	2009-2014	1997-2001

vervolgm12	2006-2015	1993-2002		2006/2009
*/


** 	Check whether there are strange things going on

*  	There is no child that has information entered on two occassions by the mother in case of v2 / v3
list PID if hschm12_v2<. & ((hschm12_v3>-1 & hschm12_v3<.) | hschm12<. | hschm12_1n<.)
list PID if (hschm12_v3>-1 & hschm12_v3<.) & 				(hschm12<. | hschm12_1n<.)
tab hschm12 hschm12_1n 
/*
	However, all values of hschm12_1n have been transferred to hschm12
	I have checked, and they have the same levels for both variables in each instance except for one case.
 	So it does not matter which variable is used. (except that hschm12 reads "dubbel antwoord"
	when someone also reported something for hschm12_2n).
	Part of hschm12 has also been transferred to hschm12_1n.
	For 1991 until 1994: all cases with value for hschm12 have no value for hschm12_1n; for 1995: 634 cases; for 1996: 29.
	For 1995, 298 cases are the same for hschm12 and hschm12_1n; for 1996: 1,090; for 1997 until 2002: all cases. 

*/
*list PID hschm12 hschm12_1n if hschm12<. & hschm12_1n<.
tab yob if hschm12<.
tab yob if hschm12<. & hschm12_1n==.
tab yob if hschm12<. & hschm12_1n<.
*	(people who have info on hschm_v2 have -1 on hschm12_v3)
*list PID if hschm12_v2<. & ( hschm12_v3<. | hschm12<. | hschm12_1n<.)

/*	
	Children that do not indicate to go to school, do not indicate
	a level either. Except for PID==22130 & 22131. I leave their level as is.
	Moreover, they indicate a level if they indicate to go to school next year: I don't
	think this is problematic.
*/
tab hschm12   	if vervolgm12==1
tab hschm12_1n 	if vervolgm12==1
tab hschm12   	if vervolgm12==-1
tab hschm12_1n 	if vervolgm12==-1
tab hschm12   	if vervolgm12==3

*	If missing on hschm12_1n, then also missing on hschm12_2n
*	If missing on hschm12_1n or hschm12_2n, then also missing on hschm12_3n
tab hschm12_2n if hschm12_1n<1 | hschm12_1n>9
tab hschm12_3n if hschm12_1n<1 | hschm12_1n>9 | hschm12_2n<1 | hschm12_2n>9 


**	Check which combination of levels occur when two options are indicated
tab hschm12_1n hschm12_2n if hschm12_3n<1
/*
-1: -1
 1:	-1				7
 2: -1
 3:	-1	4 	5		7
 4:	-1		
 5:	-1			6
 6:	-1				7
 7:	-1					8
 8:	-1						9
 9:	-1
*/

**	Check which combination of levels occur when three options are indicated
tab hschm12_3n
tab hschm12_1n hschm12_2n if hschm12_3n==5
tab hschm12_1n hschm12_2n if hschm12_3n==5 & twzyg<.
tab hschm12_1n hschm12_2n if hschm12_3n==6
tab hschm12_1n hschm12_2n if hschm12_3n==6 & twzyg<.
tab hschm12_1n hschm12_2n if hschm12_3n==8
tab hschm12_1n hschm12_2n if hschm12_3n==8 & twzyg<.
tab hschm12_1n hschm12_2n if hschm12_3n==9
tab hschm12_1n hschm12_2n if hschm12_3n==9 & twzyg<.

/*
3: 4+[5]
4:		 5+[6]
6:				7+[8]
7:						8+[9]
*/

** 	Create the actual variable 

*	Mentioned one particular level
gen edu12_m =.
replace edu12_m = .c  if hschm12_v2 ==1 | hschm12_v3 ==1  | hschm12 ==1  					// Nog niet van toepassing
replace edu12_m = .c  if vervolgm12==1 & hschm12    <=1  									// Nog niet van toepassing
replace edu12_m = .c  if vervolgm12==1 & hschm12_1n < 1  									// Nog niet van toepassing
replace edu12_m = .b  if hschm12_v2 ==2 | hschm12_v3 ==2  | hschm12 ==2  | hschm12_1n==1	// SVO
replace edu12_m = 0   if hschm12_v2 ==3 | hschm12_v3 ==3 									// LBO
replace edu12_m = 0   if 								  	hschm12 ==3  | hschm12_1n==2	// PRO
replace edu12_m = 1.5 if 				  hschm12_v3 ==5									// VMBO
replace edu12_m = 0   if 									hschm12 ==7  | hschm12_1n==6 	// VMBO basisberoepsgerichte
replace edu12_m = 1   if 									hschm12 ==6  | hschm12_1n==5 	// VMBO kaderberoepsgerichte
replace edu12_m = 2   if 									hschm12 ==5  | hschm12_1n==4 	// VMBO gemengde
replace edu12_m = 2   if 									hschm12 ==4  | hschm12_1n==3 	// VMBO theoretische
replace edu12_m = 2   if hschm12_v2 ==6 | hschm12_v3 ==9  									// MAVO
replace edu12_m = 3   if hschm12_v2 ==7 | hschm12_v3 ==10 | hschm12 ==8  | hschm12_1n==7	// HAVO
replace edu12_m = 4   if 									hschm12 ==9  | hschm12_1n==8	// VWO
replace edu12_m = 4   if hschm12_v2 ==8 | hschm12_v3 ==11									// VWO-Gymnasium
replace edu12_m = 4   if 									hschm12 ==10 | hschm12_1n==9	// Gymnasium
replace edu12_m = 1.5 if hschm12_v2 ==4 | hschm12_v3 ==4 									// MBO

* 	Explicitly mentioned two levels
replace edu12_m =  .a if hschm12_v2 ==5 | hschm12_v3 ==6 									// brugklas (anders dan beroepsonderwijs)
replace edu12_m = 1.5 if 									hschm12 ==12					// VMBO verschillende
replace edu12_m = 2.5 if 				  hschm12_v3 ==7 									// MAVO/HAVO
replace edu12_m = 2.5 if 									hschm12 ==13					// VMBO b+HAVO
replace edu12_m = 2.5 if 									hschm12 ==11 					// VMBO/HAVO
replace edu12_m = 3   if 									hschm12 ==15					// VMBO/HAVO/VWO
replace edu12_m = 3.5 if 				  hschm12_v3 ==8  | hschm12 ==16					// HAVO/VWO
replace edu12_m = 4   if 									hschm12 ==14					// VWO / Gymnasium
replace edu12_m = 0.5 if hschm12_1n ==5 & hschm12_2n ==6									// VMBO b+k
replace edu12_m = 1.5 if hschm12_1n ==4 & hschm12_2n ==5									// VMBO g+k
replace edu12_m = 1.5 if hschm12_1n ==3 & hschm12_2n ==5									// VMBO t+k
replace edu12_m = 2   if hschm12_1n ==3 & hschm12_2n ==4									// VMBO t+g
replace edu12_m = 2.5 if hschm12_1n ==6 & hschm12_2n ==7									// VMBO/HAVO
replace edu12_m = .b  if hschm12_1n ==1 & hschm12_2n ==7									// SVO / HAVO***
replace edu12_m = 2.5 if hschm12_1n ==3 & hschm12_2n ==7									// VMBO t / HAVO
replace edu12_m = 3.5 if hschm12_1n ==7 & hschm12_2n ==8									// HAVO/VWO
replace edu12_m = 4   if hschm12_1n ==8 & hschm12_2n ==9									// VWO / Gymnasium
*** This is the one case that differed between hschm12 and hschm12_1n (fixed this way).

* 	Explicitly mentioned three levels
replace edu12_m = 1   if hschm12_1n ==4 & hschm12_2n ==5 & hschm12_3n ==6					// VMBO b+k+g
replace edu12_m = 1.5 if hschm12_1n ==3 & hschm12_2n ==4 & hschm12_3n ==5					// VMBO k+g+t
replace edu12_m = 3   if hschm12_1n ==6 & hschm12_2n ==7 & hschm12_3n ==8					// VMBO b / HAVO / VWO
replace edu12_m = 3.5 if hschm12_1n ==7 & hschm12_2n ==8 & hschm12_3n ==9					// HAVO / VWO / GYM    

		
//	STEP 1.2 --------------------------------------------------------------------------------------------------------
//	Educational attainment age 12 filled in by the father												
												
** 	Explore which categories there are over the years

*  	Which survey is taken when?
tab invjrv12 if hschv12_v2<.
tab invjrv12 if hschv12_v3<.
tab invjrv12 if hschv12   <.
tab invjrv12 if hschv12_1n<.
tab invjrv12 if vervolgv12<.

*	What years of birth do each survey cover?
tab yob if hschv12_v2<.
tab yob if hschv12_v3<.
tab yob if hschv12   <.
tab yob if hschv12_1n<.
tab yob if vervolgv12<.

*   Which categories are there in each survey?
tab hschv12_v2
tab hschv12_v3
tab hschv12
tab hschv12_1n
tab vervolgv12

* 	Add labels to hschv12
label define hschv12 1 "Nog niet van toepassing" 2 "SVO" 3 "PRO" 4 "VMBO-t" 5 "VMBO-g" 6 "VMBO-k" 7 "VMBO-b" 8 "HAVO" 9 "VWO", add
label values hschv12 hschv12
tab hschv12

/*
Variabele	Invuljaar	Geboortejaar

hschv12_v2	1999-2000	1986-1987		1999		Lijst12_versie2_ouders_groenBoekje_1999			SVO/LBO/MBO/	 BRUGKLAS  /MAVO/HAVO/VWO-GYM
hschv12_v3	1999-2004	1986-1991		2000		Lijst12_versie3_ouders_paarsBoekje_2000			SVO/LBO/MBO/VMBO/BRUGKLASx3/MAVO/HAVO/VWO-GYM
hschv12		2003-2009	1991-1996		2003		Lijst12_versie4a_ouders_geelBoekje_2003 		SVO/PRO/VMBOx4/HAVO/VWO (meerdere)
	"			"			"			/2006		/Lijst12_versie5_ouders_geel-grijsBoekje_2006	SVO/PRO/VMBOx4/HAVO/VWO
hschv12_1n	2014-2017	2001-2002		/2009		/Lijst12_versie6_ouders1Lijst_druk2_2009			" "	
hschv12_?n	?
vervolgm12	2006-2009	1993-1996		2006/2009
		   +2014-2017
*/


** 	Check whether there are strange things going on

*	There is no child that has information entered on two occassions by the father 
list PID if hschv12_v2<. & ((hschv12_v3>-1 & hschv12_v3<.) | hschv12<. | hschv12_1n<.)
list PID if (hschv12_v3>-1 & hschv12_v3<.) & 				(hschv12<. | hschv12_1n<.)
list PID if hschv12<. & 												 hschv12_1n<.
*	(in contrast to mothers, here twins who have info on hschv_v2 do not have -1 on hschv12_v3)
list PID if hschv12_v2<. & ( hschv12_v3<. | hschv12<. | hschv12_1n<.)

/*	
	Children that do not indicate to go to sec. school, do not indicate
	a level either. For one exception, where mother does indicate that child goes to sec. school. 
	I leave this as is. 
*/
tab hschv12 if vervolgv12== 1
tab hschv12_1n if vervolgv12== 1
tab hschv12 if vervolgv12==-1
tab hschv12_1n if vervolgv12==-1


**	Check which combination of levels occur when more than one option is possible
/*
In this version of the data, there is not variable for the father with second option.
Will have to ask about this!!
tab hschv12_1n hschv12b

/*	
	Of the 1,521 children that indicate in their first choice that they do not go to school yet:
	- 1,019 are missing
	- 491 score -1
	- 1 scores 1
	- 10 score a real level
	on the second option (see also below). These may be children that already know which level 
	they go to next year.
	(There is only 1 child that scores 1 on the second option, scores also 1 on the first).
*/

**	Check which combination of levels occur when two options are indicated
tab hschm12_1n hschm12_2n if hschm12_3n<1
	
-2:		-1
-1:		-1		
 1:		-1	1	4				8	9	
 2:		-1
 3:		-1
 4:	 -2	-1			5	6	7	8
 5:		-1				6
 6:		-1					7
 7:		-1						8
 8:		-1							9
 9:		-1							9
*/


** 	Create the actual variable 

*	Mentioned one particular level
gen 	edu12_v = .
replace edu12_v = .c  if hschv12_v2 ==1 | hschv12_v3 ==1  | hschv12 ==1  					// Nog niet van toepassing
replace edu12_v = .c  if vervolgv12==1 & hschv12_1n <1										// Nog niet van toepassing
replace edu12_v = .c  if vervolgv12==1 & hschv12 ==. & hschv12_1n ==.						// Nog niet van toepassing

replace edu12_v = .b  if hschv12_v2 ==2 | hschv12_v3 ==2  | hschv12 ==2  | hschv12_1n ==1	// SVO
replace edu12_v = 0   if hschv12_v2 ==3 | hschv12_v3 ==3 									// LBO
replace edu12_v = 0   if 								    hschv12 ==3  | hschv12_1n ==2	// PRO
replace edu12_v = 0   if 								    hschv12 ==7  | hschv12_1n ==6	// VMBO basisberoepsgerichte
replace edu12_v = 1   if 								    hschv12 ==6  | hschv12_1n ==5	// VMBO kaderberoepsgerichte
replace edu12_v = 2   if 								    hschv12 ==5  | hschv12_1n ==4	// VMBO gemengde
replace edu12_v = 2   if 								    hschv12 ==4  | hschv12_1n ==3	// VMBO theoretische
replace edu12_v = 1.5 if 				  hschv12_v3 ==5									// VMBO
replace edu12_v = 2   if hschv12_v2 ==6 | hschv12_v3 ==9  									// MAVO
replace edu12_v = 3   if hschv12_v2 ==7 | hschv12_v3 ==10 | hschv12 ==8  | hschv12_1n ==7 	// HAVO
replace edu12_v = 4   if 									hschv12 ==9  | hschv12_1n ==8   // VWO
replace edu12_v = 4   if hschv12_v2 ==8 | hschv12_v3 ==11 									// VWO-Gymnasium
replace edu12_v = 4   if 									hschv12 ==10 | hschv12_1n ==9	// Gymnasium
replace edu12_v = 1.5 if hschv12_v2 ==4 | hschv12_v3 ==4 									// MBO

* 	Explicitly mentioned to be in a class with more than one level
replace edu12_v =  .a if hschv12_v2 ==5 | hschv12_v3 ==6 									// brugklas (anders dan beroepsonderwijs)
replace edu12_v = 1.5 if hschv12 ==12														// VMBO verschillende leerwegen
replace edu12_v = 2.5 if hschv12 ==13 														// VMBO b+HAVO
replace edu12_v = 2.5 if hschv12 ==11  														// VMBO (t)+HAVO
replace edu12_v = 2.5 if 				  hschv12_v3 ==7 									// MAVO/HAVO
replace edu12_m = 3   if hschv12 ==15														// VMBO/HAVO/VWO
replace edu12_v = 3.5 if 				  hschv12_v3 ==8 									// HAVO/VWO
replace edu12_v = 3.5 if hschv12 ==16 														//  " "
replace edu12_v = 4   if hschv12 ==14														// VWO/GYM


//	STEP 1.3 --------------------------------------------------------------------------------------------------------
//	Educational attainment age 12, combination of mother's and father's report	

*	The mean and distribution of mother and father's report are pretty similar
codebook edu12_m edu12_v		
tab edu12_m
tab edu12_v 

*	Check how much parents differ in their reports on the individual level
tab edu12_m edu12_v
gen edu_dif_mv = edu12_m-edu12_v
tab edu_dif_mv
gen edu_adif_mv = abs(edu_dif_mv)
tab edu_adif_mv

/*	
	If both parents report a level:
	- 87.7% report the same level 
	- 95.5% report a difference of 1 or less (+7.8)
	- 98.9% report a difference of 2 or less (+3.4)
	- 99.8% report a difference of 3 or less (+0.9)
	Fathers and mothers report a higher level to the same extent (about 6.1% each). 
*/

/*
*	Get a feeling for in what kind of cases parents report different edu levels
list agem12 invjrm12 hschm12 	agev12 invjrv12 hschv12 hschv12_1n 	if edu_dif_mv==-5
list agem12 invjrm12 hschm12_v3 agev12 invjrv12 hschv12_v3 			if edu_dif_mv==-4.5
list agem12 invjrm12 hschm12_v2 agev12 invjrv12 hschv12_v2 			if edu_dif_mv==-3 & hschm12_v2<.
list agem12 invjrm12 hschm12_v3 agev12 invjrv12 hschv12_v3 			if edu_dif_mv==-3 & hschm12_v3<. & hschm12_v3>-1
list agem12 invjrm12 hschm12 	agev12 invjrv12 hschv12 hschv12_1n  if edu_dif_mv==-3 & hschm12<.

list agem12 invjrm12 hschm12_v2 agev12 invjrv12 hschv12_v2 if edu12_m!=edu12_v & (edu12_m<. & edu12_v<.) & (hschm12_v2<. | hschv12_v2<.)
list agem12 invjrm12 hschm12_v3 agev12 invjrv12 hschv12_v3 if edu12_m!=edu12_v & (edu12_m<. & edu12_v<.) & (hschm12_v3<. | hschv12_v3<.)
*/

*	Variable that simply takes average of what both parents say
egen edu12 = rowmean(edu12_m edu12_v)
replace edu12 = .a if (edu12_m==.a & (edu12_v==.a | edu12_v==.b | edu12_v==.c | edu12_v==.)) | (edu12_v==.a & (edu12_m==.b | edu12_m==.c | edu12_m==.))
replace edu12 = .b if (edu12_m==.b & (edu12_v==.b | edu12_v==.c | edu12_v==.))               | (edu12_v==.b & (edu12_m==.c | edu12_m==.))
replace edu12 = .c if (edu12_m==.c & (edu12_v==.c | edu12_v==.))                             | (edu12_v==.c & edu12_m==.)
codebook edu12
tab edu12


//	STEP 1.4 --------------------------------------------------------------------------------------------------------
//	Educational attainment age 14 filled in by student
 
** 	Explore which categories there are over the years

*  	Which survey is taken when?
tab invjaars14 if onderwijss14  <.
tab invjaars14 if onderw1s14    <.
tab invjaars14 if onderw2s14    <.
tab invjaars14 if onderw3s14    <.
tab invjaars14 if onderwijsNs14 <.

*	What years of birth do each survey cover?
tab yob if onderwijss14  <. & twzyg<.
tab yob if onderw1s14    <. & twzyg<.
tab yob if onderwijsNs14 <. & twzyg<.

tab yob if onderwijss14  <. 
tab yob if onderw1s14    <. 
tab yob if onderwijsNs14 <. 

*   Which categories are there in each survey?
tab onderwijss14 	if twzyg<.
tab onderw1s14		if twzyg<.
tab onderw2s14		if twzyg<.
tab onderw3s14		if twzyg<.
tab onderw4s14		if twzyg<.
tab onderw5s14		if twzyg<.
tab onderwijsNs14	if twzyg<.

tab onderwijss14
tab onderw1s14
tab onderw2s14
tab onderw3s14
tab onderw4s14
tab onderw5s14
tab onderwijsNs14


/*
Variabele		Invuljaar	Geboortejaar

onderwijss14	2005-2010	1989-1994 (1979-1998)	2004		DHBQ_14-16-18_versie1_rood-grijsBoekje			VMBO (mavo, lbo)/HAVO/VWO-GYM
													2005		DHBQ_14-16-18_versie2_blauwBoekje_druk1
onderw1s14		2008-2012	1993-1998 (1983-1999)	2009		DHBQ14_versie3_ONLINE_2009						VMBO/VMBO-t/HAVO/VWO/Gymnasium 
onderwXs14			""								2011		DHBQ14_versie4_ONLINE_2011						+ 	MBO1/MBO2-4/HBO/UNI/Anders, nl
													2012		DHBQ14_versie5_ONLINE_2012							(meerdere)
onderwijsNs14	2013		1998-1999 (1989-2001)	2013		DHBQ14_versie6_ONLINE_2013 						VMBO/VMBO-t/HAVO/VWO/Gymnasium 
*/


** 	Check whether there are strange things going on

*  	There is no child that has information entered on two occassions  
list PID if onderwijss14<. & (onderw1s14<. | onderwijsNs14<.)
list PID if onderw1s14<. & onderwijsNs14<.		

/*	
	Should ask NTR for variables that indicate whether someone is still in school,
	and about which diplomas were attained.
*/


**	Check which combination of levels occur when more than one option is possible

* 	First those who give 2 answers, N=522/(651)
tab onderw1s14 onderw2s14 if onderw2s14>-1 & onderw3s14==-1
tab onderw1s14 onderw2s14 if onderw2s14>-1 & onderw3s14==-1 & twzyg<.
* 	(I've checked: when people have -1 on onderw(X)s14, they have -1 on onderw(X+1)s14.

/*
 1 :	2	3	4	   (6)	7			10
 2 : 		3		   (6)	7			10
 3 :			4	5		7  (8) (9)	10
 4 :				5		   (8) (9)	10
 5 :							   (9)	10
(7):						   (8)	   (10)	

Between brackets if siblings are included
*/

* 	Those who give 3 answers, N=42/(75)
tab onderw2s14 onderw3s14 if onderw1s14==1  & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==1  & onderw3s14>-1 & onderw4s14==-1 & twzyg<.
tab onderw2s14 onderw3s14 if onderw1s14==2  & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==2  & onderw3s14>-1 & onderw4s14==-1 & twzyg<.
tab onderw2s14 onderw3s14 if onderw1s14==3  & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==3  & onderw3s14>-1 & onderw4s14==-1 & twzyg<.
tab onderw2s14 onderw3s14 if onderw1s14==4  & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==4  & onderw3s14>-1 & onderw4s14==-1 & twzyg<.
tab onderw2s14 onderw3s14 if onderw1s14==5  & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==5  & onderw3s14>-1 & onderw4s14==-1 & twzyg<.
tab onderw2s14 onderw3s14 if onderw1s14==6  & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==6  & onderw3s14>-1 & onderw4s14==-1 & twzyg<.
tab onderw2s14 onderw3s14 if onderw1s14==7  & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==7  & onderw3s14>-1 & onderw4s14==-1 & twzyg<.
tab onderw2s14 onderw3s14 if onderw1s14==8  & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==8  & onderw3s14>-1 & onderw4s14==-1 & twzyg<.
tab onderw2s14 onderw3s14 if onderw1s14==9  & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==9  & onderw3s14>-1 & onderw4s14==-1 & twzyg<.
tab onderw2s14 onderw3s14 if onderw1s14==10 & onderw3s14>-1 & onderw4s14==-1
tab onderw2s14 onderw3s14 if onderw1s14==10 & onderw3s14>-1 & onderw4s14==-1 & twzyg<.


/*
 1 :	2+[3 4(6)(7)10]  3+[4     (7)(8)10]									6+[(7)]	 7+[(8)]    
 2 : 			         3+[4 5(6)(7)(8)10]											 7+[(8)]
 3 :										  4+[5(8)(9)10]		
 4 :														 5+[(9) 10 ]    8+[(9)]
 5 :														 
 6 :		
 7 :		
 8 :
 9 :
 10:
 
Between brackets if siblings are included
*/

*	Those who give 4 answers, N=1/(4)
list onderw1s14 onderw2s14 onderw3s14 onderw4s14 			if onderw4s14>-1 & onderw5s14==-1 
list onderw1s14 onderw2s14 onderw3s14 onderw4s14 			if onderw4s14>-1 & onderw5s14==-1 & twzyg<.
/*
 2+3+4+[(7) 10] 
(3  +    6+7+8)
*/

list onderw1s14 onderw2s14 onderw3s14 onderw4s14 onderw5s14 if onderw5s14>-1 & onderw5s14<.
list onderw1s14 onderw2s14 onderw3s14 onderw4s14 onderw5s14 if onderw5s14>-1 & onderw5s14<. & twzyg<.

*	Those who give 5 answers, N=1/(1)
* 1+3+4+5+10


** 	Create the actual variable 

*	Mentioned one particular level
gen edu14 =.
replace edu14 = 0.5 if 					 (onderw1s14== 1 & onderw2s14==-1) | onderwijsNs14==1	// VMBO
replace edu14 = 1.5 if onderwijss14==1															// VMBO (mavo, lbo)
replace edu14 = 2   if 					 (onderw1s14== 2 & onderw2s14==-1) | onderwijsNs14==2 	// VMBO theoretische
replace edu14 = 3   if onderwijss14==2 | (onderw1s14== 3 & onderw2s14==-1) | onderwijsNs14==3 	// HAVO
replace edu14 = 4   if 					 (onderw1s14== 4 & onderw2s14==-1) | onderwijsNs14==4	// VWO
replace edu14 = 4   if onderwijss14==3 															// VWO-Gymnasium
replace edu14 = 4   if 					 (onderw1s14== 5 & onderw2s14==-1) | onderwijsNs14==5	// Gymnasium
replace edu14 = 0   if 					 (onderw1s14== 6 & onderw2s14==-1) 						// MBO1
replace edu14 = 1.5 if 					 (onderw1s14== 7 & onderw2s14==-1) 						// MBO2-4
replace edu14 = 3   if 					 (onderw1s14== 8 & onderw2s14==-1) 						// HBO
replace edu14 = 4   if 					 (onderw1s14== 9 & onderw2s14==-1) 						// Universiteit
replace edu14 = .a  if 					 (onderw1s14==10 & onderw2s14==-1) 						// Anders

* 	Mentioned two levels
replace edu14 = 1.5 if onderw1s14==1 & onderw2s14== 2 & onderw3s14==-1							// VMBO   + VMBO-t
replace edu14 = 2.5 if onderw1s14==1 & onderw2s14== 3 & onderw3s14==-1							// VMBO   + HAVO
replace edu14 = 3   if onderw1s14==1 & onderw2s14== 4 & onderw3s14==-1							// VMBO   + VWO
replace edu14 = 0   if onderw1s14==1 & onderw2s14== 6 & onderw3s14==-1							// VMBO   + MBO1
replace edu14 = 1   if onderw1s14==1 & onderw2s14== 7 & onderw3s14==-1							// VMBO   + MBO2-4
replace edu14 = 0.5 if onderw1s14==1 & onderw2s14==10 & onderw3s14==-1							// VMBO   + Anders
replace edu14 = 2.5 if onderw1s14==2 & onderw2s14== 3 & onderw3s14==-1							// VMBO-t + HAVO
replace edu14 = 2   if onderw1s14==2 & onderw2s14== 6 & onderw3s14==-1							// VMBO-t + MBO1
replace edu14 = 2   if onderw1s14==2 & onderw2s14== 7 & onderw3s14==-1							// VMBO-t + MBO2-4
replace edu14 = 2   if onderw1s14==2 & onderw2s14==10 & onderw3s14==-1							// VMBO-t + Anders
replace edu14 = 3.5 if onderw1s14==3 & onderw2s14== 4 & onderw3s14==-1							// HAVO   + VWO
replace edu14 = 3.5 if onderw1s14==3 & onderw2s14== 5 & onderw3s14==-1							// HAVO   + Gymnasium
replace edu14 = 3   if onderw1s14==3 & onderw2s14== 7 & onderw3s14==-1							// HAVO   + MBO2-4
replace edu14 = 3   if onderw1s14==3 & onderw2s14== 8 & onderw3s14==-1							// HAVO   + HBO
replace edu14 = 3   if onderw1s14==3 & onderw2s14== 9 & onderw3s14==-1							// HAVO   + Universiteit
replace edu14 = 3   if onderw1s14==3 & onderw2s14==10 & onderw3s14==-1							// HAVO   + Anders
replace edu14 = 4   if onderw1s14==4 & onderw2s14== 5 & onderw3s14==-1							// VWO    + Gymnasium
replace edu14 = 4   if onderw1s14==4 & onderw2s14== 8 & onderw3s14==-1							// VWO    + HBO
replace edu14 = 4   if onderw1s14==4 & onderw2s14== 9 & onderw3s14==-1							// VWO    + Universiteit
replace edu14 = 4   if onderw1s14==4 & onderw2s14==10 & onderw3s14==-1							// VWO    + Anders
replace edu14 = 4   if onderw1s14==5 & onderw2s14== 9 & onderw3s14==-1							// Gymna. + Universiteit
replace edu14 = 4   if onderw1s14==5 & onderw2s14==10 & onderw3s14==-1							// Gymna. + Anders
replace edu14 = 1.5 if onderw1s14==6 & onderw2s14== 7 & onderw3s14==-1							// MBO1   + MBO2-4
replace edu14 = 2   if onderw1s14==7 & onderw2s14== 8 & onderw3s14==-1							// MBO2-4 + HBO
replace edu14 = 1.5 if onderw1s14==7 & onderw2s14==10 & onderw3s14==-1							// MBO2-4 + Anders

* 	Mentioned three levels
replace edu14 = 2   if onderw1s14==1 & onderw2s14== 2 & onderw3s14== 3							// VMBO   + VMBO-t 	+ HAVO
replace edu14 = 3   if onderw1s14==1 & onderw2s14== 2 & onderw3s14== 4							// VMBO   + VMBO-t 	+ VWO
replace edu14 = 1   if onderw1s14==1 & onderw2s14== 2 & onderw3s14== 6							// VMBO   + VMBO-t 	+ MBO1
replace edu14 = 1.5 if onderw1s14==1 & onderw2s14== 2 & onderw3s14== 7							// VMBO   + VMBO-t 	+ MBO2-4
replace edu14 = 1.5 if onderw1s14==1 & onderw2s14== 2 & onderw3s14==10							// VMBO   + VMBO-t 	+ Anders
replace edu14 = 3   if onderw1s14==1 & onderw2s14== 3 & onderw3s14== 4							// VMBO   + HAVO 	+ VWO
replace edu14 = 2.5 if onderw1s14==1 & onderw2s14== 3 & onderw3s14== 7							// VMBO   + HAVO 	+ MBO2-4
replace edu14 = 3   if onderw1s14==1 & onderw2s14== 3 & onderw3s14== 8							// VMBO   + HAVO 	+ HBO
replace edu14 = 2.5 if onderw1s14==1 & onderw2s14== 3 & onderw3s14==10							// VMBO   + HAVO 	+ Anders
replace edu14 = 2   if onderw1s14==1 & onderw2s14== 7 & onderw3s14== 8							// VMBO   + MBO2-4 	+ HBO
replace edu14 = 2.5 if onderw1s14==1 & onderw2s14== 8 & onderw3s14== 9							// VMBO   + HBO 	+ Uni
replace edu14 = 3   if onderw1s14==2 & onderw2s14== 3 & onderw3s14== 4							// VMBO-t + HAVO	+ VWO
replace edu14 = 3   if onderw1s14==2 & onderw2s14== 3 & onderw3s14== 5							// VMBO-t + HAVO	+ Gym
replace edu14 = 2.5 if onderw1s14==2 & onderw2s14== 3 & onderw3s14== 6							// VMBO-t + HAVO	+ MBO1
replace edu14 = 2.5 if onderw1s14==2 & onderw2s14== 3 & onderw3s14== 7							// VMBO-t + HAVO	+ MBO2-4
replace edu14 = 2.5 if onderw1s14==2 & onderw2s14== 3 & onderw3s14== 8							// VMBO-t + HAVO	+ HBO
replace edu14 = 2.5 if onderw1s14==2 & onderw2s14== 3 & onderw3s14==10							// VMBO-t + HAVO	+ Anders
replace edu14 = 2   if onderw1s14==2 & onderw2s14== 7 & onderw3s14== 8							// VMBO-t + MBO2-4	+ HBO
replace edu14 = 2   if onderw1s14==2 & onderw2s14== 7 & onderw3s14==10							// VMBO-t + MBO2-4	+ Anders
replace edu14 = 2   if onderw1s14==2 & onderw2s14==10 & onderw3s14==10							// VMBO-t + Anders	+ Anders
replace edu14 = 3.5 if onderw1s14==3 & onderw2s14== 4 & onderw3s14== 5							// HAVO	  + VWO		+ Gym
replace edu14 = 3.5 if onderw1s14==3 & onderw2s14== 4 & onderw3s14== 8							// HAVO	  + VWO 	+ HBO
replace edu14 = 3.5 if onderw1s14==3 & onderw2s14== 4 & onderw3s14== 9							// HAVO	  + VWO		+ Uni
replace edu14 = 3.5 if onderw1s14==3 & onderw2s14== 4 & onderw3s14==10							// HAVO	  + VWO		+ Anders
replace edu14 = 4   if onderw1s14==4 & onderw2s14== 5 & onderw3s14== 9							// VWO	  +	Gym		+ Uni
replace edu14 = 4   if onderw1s14==4 & onderw2s14== 5 & onderw3s14==10							// VWO	  +	Gym		+ Anders
replace edu14 = 4   if onderw1s14==4 & onderw2s14== 8 & onderw3s14== 9							// VWO	  +	HBO		+ Uni

* 	Mentioned four levels
replace edu14 = 2.5 if onderw1s14==2 & onderw2s14== 3 & onderw3s14== 4 & onderw4s14== 7			// VMBO-t + HAVO	+ VWO 	 + MBO2-4
replace edu14 = 3   if onderw1s14==2 & onderw2s14== 3 & onderw3s14== 4 & onderw4s14==10			// VMBO-t + HAVO	+ VWO 	 + Anders
replace edu14 = 3   if onderw1s14==3 & onderw2s14== 6 & onderw3s14== 7 & onderw4s14== 8			// HAVO	  + MBO1 	+ MBO2-4 + HBO

* 	Mentioned five levels
replace edu14 = 3   if onderw1s14==1 & onderw2s14== 3 & onderw3s14== 4	& onderw4s14== 5 & onderw5s14== 10	// VMBO + HAVO + VWO + Gym + Anders

codebook edu14
codebook edu14 if edu12<.


//	STEP 1.5 --------------------------------------------------------------------------------------------------------
//	Educational attainment age 16 filled in by student

** 	Explore which categories there are over the years

*  	Which survey is taken when?
tab invjaars16 if onderwijss16  <.
tab invjaars16 if onderw1s16    <.
tab invjaars16 if onderw2s16    <.
tab invjaars16 if onderw3s16    <.
tab invjaars16 if onderw4s16    <.
tab invjaars16 if onderw5s16    <.
tab invjaars16 if onderwijsNs16 <.

*	What years of birth do each survey cover?
tab yob if onderwijss16  <. & twzyg<.
tab yob if onderw1s16    <. & twzyg<.
tab yob if onderwijsNs16 <. & twzyg<.

tab yob if onderwijss16  <. 
tab yob if onderw1s16    <. 
tab yob if onderwijsNs16 <. 

*   Which categories are there in each survey?
tab onderwijss16 	if twzyg<.
tab onderw1s16		if twzyg<.
tab onderw2s16		if twzyg<.
tab onderw3s16		if twzyg<.
tab onderw4s16		if twzyg<.
tab onderw5s16		if twzyg<.
tab onderwijsNs16	if twzyg<.

tab onderwijss16
tab onderw1s16
tab onderw2s16
tab onderw3s16
tab onderw4s16
tab onderw5s16
tab onderwijsNs16


/*
Variabele		Invuljaar	Geboortejaar

onderwijss16	2005-2008	1988-1990 (1981-1995)	2004		DHBQ_14-16-18_versie1_rood-grijsBoekje			VMBO/HAVO/VWO-GYM
													2005		DHBQ_14-16-18_versie2_blauwBoekje_druk1
onderw1s16		2009-2014	1991-1996 (1983-1999)	2009		DHBQ14_versie3_ONLINE_2009						VMBO/VMBO-t/HAVO/VWO/Gymnasium 
onderwXs16			""								2011		DHBQ14_versie4_ONLINE_2011						+ 	MBO1/MBO2-4/HBO/UNI/Anders, nl
													2012		DHBQ14_versie5_ONLINE_2012							(meerdere)
onderwijsNs16	2013		1995-1997 (1994-1999)	2013		DHBQ14_versie6_ONLINE_2013 						VMBO/VMBO-t/HAVO/VWO/Gymnasium 
*/


** 	Check whether there are strange things going on

*  	There is no child that has information entered on two occassions  
list PID if onderwijss16<. & (onderw1s16<. | onderwijsNs16<.)
list PID if onderw1s16<. & onderwijsNs16<.		


**	Check which combination of levels occur when more than one option is possible

* 	First those who give 2 answers, N=431/(510)
tab onderw1s16 onderw2s16 if onderw2s16>-1 & onderw3s16==-1
tab onderw1s16 onderw2s16 if onderw2s16>-1 & onderw3s16==-1 & twzyg<.
* 	(I've checked: when people have -1 on onderw(X)s16, they have -1 on onderw(X+1)s16.

/*
 1 :	2	3		    6	7	8		10
 2 : 		3	4	    6	7	8		10
 3 :			4	5	6	7   8  	
 4 :				5		   (8)  9	10
 5 :							    9	10
 6 :				        7
 7 :						   (8)	    10	

Between brackets if siblings are included
*/

* 	Those who give 3 answers, N=23/(42)
tab onderw2s16 onderw3s16 if onderw1s16==1  & onderw3s16>-1 & onderw4s16==-1
tab onderw2s16 onderw3s16 if onderw1s16==1  & onderw3s16>-1 & onderw4s16==-1 & twzyg<.
tab onderw2s16 onderw3s16 if onderw1s16==2  & onderw3s16>-1 & onderw4s16==-1
tab onderw2s16 onderw3s16 if onderw1s16==2  & onderw3s16>-1 & onderw4s16==-1 & twzyg<.
tab onderw2s16 onderw3s16 if onderw1s16==3  & onderw3s16>-1 & onderw4s16==-1
tab onderw2s16 onderw3s16 if onderw1s16==3  & onderw3s16>-1 & onderw4s16==-1 & twzyg<.
tab onderw2s16 onderw3s16 if onderw1s16==4  & onderw3s16>-1 & onderw4s16==-1
tab onderw2s16 onderw3s16 if onderw1s16==4  & onderw3s16>-1 & onderw4s16==-1 & twzyg<.
tab onderw2s16 onderw3s16 if onderw1s16==5  & onderw3s16>-1 & onderw4s16==-1
tab onderw2s16 onderw3s16 if onderw1s16==5  & onderw3s16>-1 & onderw4s16==-1 & twzyg<.
tab onderw2s16 onderw3s16 if onderw1s16==6  & onderw3s16>-1 & onderw4s16==-1				
tab onderw2s16 onderw3s16 if onderw1s16==6  & onderw3s16>-1 & onderw4s16==-1 & twzyg<. 		
tab onderw2s16 onderw3s16 if onderw1s16==7  & onderw3s16>-1 & onderw4s16==-1
tab onderw2s16 onderw3s16 if onderw1s16==7  & onderw3s16>-1 & onderw4s16==-1 & twzyg<.
tab onderw2s16 onderw3s16 if onderw1s16==8  & onderw3s16>-1 & onderw4s16==-1
tab onderw2s16 onderw3s16 if onderw1s16==8  & onderw3s16>-1 & onderw4s16==-1 & twzyg<.
tab onderw2s16 onderw3s16 if onderw1s16==9  & onderw3s16>-1 & onderw4s16==-1
tab onderw2s16 onderw3s16 if onderw1s16==9  & onderw3s16>-1 & onderw4s16==-1 & twzyg<.
tab onderw2s16 onderw3s16 if onderw1s16==10 & onderw3s16>-1 & onderw4s16==-1
tab onderw2s16 onderw3s16 if onderw1s16==10 & onderw3s16>-1 & onderw4s16==-1 & twzyg<.


/*
 1 :  2+[3   6 7]   3+[     (8)10]    								6+[7]	  (7+[ 8	 ])    
 2 :  			    3+[4 6 7 8 10]									6+[7]	   7+[(8)(10)]		
 3 :								4+[5(8)]						6+[8] 	   7+[    10 ]	(8+[9])
 4 :											5+[(7)(8) 10 ]  							(8+[9])
 5 :																 									(9+[10])
 6 :		
 7 :		
 8 :
 9 :
 10:
 
Between brackets if siblings are included
*/

*	Those who give 4 answers, N=2/(4)
list onderw1s16 onderw2s16 onderw3s16 onderw4s16 			if onderw4s16>-1 & onderw5s16==-1 
list onderw1s16 onderw2s16 onderw3s16 onderw4s16 			if onderw4s16>-1 & onderw5s16==-1 & twzyg<.
/*
 1+2+3   +6
   2+3+4 +6 
*/

list onderw1s16 onderw2s16 onderw3s16 onderw4s16 onderw5s16 if onderw5s16>-1 & onderw5s16<.
*	Those who give 5 answers, N=1/(1)
* 1+2+3+6+7


** 	Create the actual variable 

*	Mentioned one particular level
gen edu16 =.
replace edu16 = 0.5 if 					 (onderw1s16== 1 & onderw2s16==-1) | onderwijsNs16==1	// VMBO (vmbo-t ook optie)
replace edu16 = 1.5 if onderwijss16==1															// VMBO (enige optie)
replace edu16 = 2   if 					 (onderw1s16== 2 & onderw2s16==-1) | onderwijsNs16==2 	// VMBO theoretische
replace edu16 = 3   if onderwijss16==2 | (onderw1s16== 3 & onderw2s16==-1) | onderwijsNs16==3 	// HAVO
replace edu16 = 4   if 					 (onderw1s16== 4 & onderw2s16==-1) | onderwijsNs16==4	// VWO
replace edu16 = 4   if onderwijss16==3 															// VWO-Gymnasium
replace edu16 = 4   if 					 (onderw1s16== 5 & onderw2s16==-1) | onderwijsNs16==5	// Gymnasium
replace edu16 = 0   if 					 (onderw1s16== 6 & onderw2s16==-1) 						// MBO1
replace edu16 = 1.5 if 					 (onderw1s16== 7 & onderw2s16==-1) 						// MBO2-4
replace edu16 = 3   if 					 (onderw1s16== 8 & onderw2s16==-1) 						// HBO
replace edu16 = 4   if 					 (onderw1s16== 9 & onderw2s16==-1) 						// Universiteit
replace edu16 = .a  if 					 (onderw1s16==10 & onderw2s16==-1) 						// Anders

* 	Mentioned two levels
replace edu16 = 1.5 if onderw1s16==1 & onderw2s16== 2 & onderw3s16==-1							// VMBO   + VMBO-t
replace edu16 = 2.5 if onderw1s16==1 & onderw2s16== 3 & onderw3s16==-1							// VMBO   + HAVO
replace edu16 = 3   if onderw1s16==1 & onderw2s16== 4 & onderw3s16==-1							// VMBO   + VWO
replace edu16 = 0   if onderw1s16==1 & onderw2s16== 6 & onderw3s16==-1							// VMBO   + MBO1
replace edu16 = 1   if onderw1s16==1 & onderw2s16== 7 & onderw3s16==-1							// VMBO   + MBO2-4
replace edu16 = 2   if onderw1s16==1 & onderw2s16== 8 & onderw3s16==-1							// VMBO   + HBO
replace edu16 = 0.5 if onderw1s16==1 & onderw2s16==10 & onderw3s16==-1							// VMBO   + Anders
replace edu16 = 2.5 if onderw1s16==2 & onderw2s16== 3 & onderw3s16==-1							// VMBO-t + HAVO
replace edu16 = 3   if onderw1s16==2 & onderw2s16== 4 & onderw3s16==-1							// VMBO-t + VWO
replace edu16 = 2   if onderw1s16==2 & onderw2s16== 6 & onderw3s16==-1							// VMBO-t + MBO1
replace edu16 = 2   if onderw1s16==2 & onderw2s16== 7 & onderw3s16==-1							// VMBO-t + MBO2-4
replace edu16 = 2   if onderw1s16==2 & onderw2s16== 8 & onderw3s16==-1							// VMBO-t + HBO
replace edu16 = 2   if onderw1s16==2 & onderw2s16==10 & onderw3s16==-1							// VMBO-t + Anders
replace edu16 = 3.5 if onderw1s16==3 & onderw2s16== 4 & onderw3s16==-1							// HAVO   + VWO
replace edu16 = 3.5 if onderw1s16==3 & onderw2s16== 5 & onderw3s16==-1							// HAVO   + Gymnasium
replace edu16 = 3   if onderw1s16==3 & onderw2s16== 6 & onderw3s16==-1							// HAVO   + MBO1
replace edu16 = 3   if onderw1s16==3 & onderw2s16== 7 & onderw3s16==-1							// HAVO   + MBO2-4
replace edu16 = 3   if onderw1s16==3 & onderw2s16== 8 & onderw3s16==-1							// HAVO   + HBO
replace edu16 = 3   if onderw1s16==3 & onderw2s16== 9 & onderw3s16==-1							// HAVO   + Universiteit
replace edu16 = 3   if onderw1s16==3 & onderw2s16==10 & onderw3s16==-1							// HAVO   + Anders
replace edu16 = 4   if onderw1s16==4 & onderw2s16== 5 & onderw3s16==-1							// VWO    + Gymnasium
replace edu16 = 4   if onderw1s16==4 & onderw2s16== 8 & onderw3s16==-1							// VWO    + HBO
replace edu16 = 4   if onderw1s16==4 & onderw2s16== 9 & onderw3s16==-1							// VWO    + Universiteit
replace edu16 = 4   if onderw1s16==4 & onderw2s16==10 & onderw3s16==-1							// VWO    + Anders
replace edu16 = 4   if onderw1s16==5 & onderw2s16== 9 & onderw3s16==-1							// Gymna. + Universiteit
replace edu16 = 4   if onderw1s16==5 & onderw2s16==10 & onderw3s16==-1							// Gymna. + Anders
replace edu16 = 1.5 if onderw1s16==6 & onderw2s16== 7 & onderw3s16==-1							// MBO1   + MBO2-4
replace edu16 = 2   if onderw1s16==7 & onderw2s16== 8 & onderw3s16==-1							// MBO2-4 + HBO
replace edu16 = 1.5 if onderw1s16==7 & onderw2s16==10 & onderw3s16==-1							// MBO2-4 + Anders


* 	Mentioned three levels
replace edu16 = 2   if onderw1s16==1 & onderw2s16== 2 & onderw3s16== 3	& onderw4s16==.			// VMBO   + VMBO-t 	+ HAVO
replace edu16 = 3   if onderw1s16==1 & onderw2s16== 2 & onderw3s16== 4	& onderw4s16==.			// VMBO   + VMBO-t 	+ VWO
replace edu16 = 1   if onderw1s16==1 & onderw2s16== 2 & onderw3s16== 6	& onderw4s16==.			// VMBO   + VMBO-t 	+ MBO1
replace edu16 = 1.5 if onderw1s16==1 & onderw2s16== 2 & onderw3s16== 7	& onderw4s16==.			// VMBO   + VMBO-t 	+ MBO2-4
replace edu16 = 1.5 if onderw1s16==1 & onderw2s16== 2 & onderw3s16==10	& onderw4s16==.			// VMBO   + VMBO-t 	+ Anders
replace edu16 = 3   if onderw1s16==1 & onderw2s16== 3 & onderw3s16== 4	& onderw4s16==.			// VMBO   + HAVO 	+ VWO
replace edu16 = 2.5 if onderw1s16==1 & onderw2s16== 3 & onderw3s16== 7	& onderw4s16==.			// VMBO   + HAVO 	+ MBO2-4
replace edu16 = 3   if onderw1s16==1 & onderw2s16== 3 & onderw3s16== 8	& onderw4s16==.			// VMBO   + HAVO 	+ HBO
replace edu16 = 2.5 if onderw1s16==1 & onderw2s16== 3 & onderw3s16==10	& onderw4s16==.			// VMBO   + HAVO 	+ Anders
replace edu16 = 1   if onderw1s16==1 & onderw2s16== 6 & onderw3s16== 7	& onderw4s16==.			// VMBO   + MBO1 	+ MBO2-4
replace edu16 = 2   if onderw1s16==1 & onderw2s16== 7 & onderw3s16== 8	& onderw4s16==.			// VMBO   + MBO2-4 	+ HBO
replace edu16 = 2.5 if onderw1s16==1 & onderw2s16== 8 & onderw3s16== 9	& onderw4s16==.			// VMBO   + HBO 	+ Uni
replace edu16 = 3   if onderw1s16==2 & onderw2s16== 3 & onderw3s16== 4	& onderw4s16==.			// VMBO-t + HAVO	+ VWO
replace edu16 = 3   if onderw1s16==2 & onderw2s16== 3 & onderw3s16== 5	& onderw4s16==.			// VMBO-t + HAVO	+ Gym
replace edu16 = 2.5 if onderw1s16==2 & onderw2s16== 3 & onderw3s16== 6	& onderw4s16==.			// VMBO-t + HAVO	+ MBO1
replace edu16 = 2.5 if onderw1s16==2 & onderw2s16== 3 & onderw3s16== 7	& onderw4s16==.			// VMBO-t + HAVO	+ MBO2-4
replace edu16 = 2.5 if onderw1s16==2 & onderw2s16== 3 & onderw3s16== 8	& onderw4s16==.			// VMBO-t + HAVO	+ HBO
replace edu16 = 2.5 if onderw1s16==2 & onderw2s16== 3 & onderw3s16==10	& onderw4s16==.			// VMBO-t + HAVO	+ Anders
replace edu16 = 2   if onderw1s16==2 & onderw2s16== 6 & onderw3s16== 7	& onderw4s16==.			// VMBO-t + MBO1 	+ MBO2-4
replace edu16 = 2   if onderw1s16==2 & onderw2s16== 7 & onderw3s16== 8	& onderw4s16==.			// VMBO-t + MBO2-4	+ HBO
replace edu16 = 2   if onderw1s16==2 & onderw2s16== 7 & onderw3s16==10	& onderw4s16==.			// VMBO-t + MBO2-4	+ Anders
replace edu16 = 2   if onderw1s16==2 & onderw2s16==10 & onderw3s16==10	& onderw4s16==.			// VMBO-t + Anders	+ Anders
replace edu16 = 3.5 if onderw1s16==3 & onderw2s16== 4 & onderw3s16== 5	& onderw4s16==.			// HAVO	  + VWO		+ Gym
replace edu16 = 3.5 if onderw1s16==3 & onderw2s16== 4 & onderw3s16== 8	& onderw4s16==.			// HAVO	  + VWO 	+ HBO
replace edu16 = 3.5 if onderw1s16==3 & onderw2s16== 4 & onderw3s16== 9	& onderw4s16==.			// HAVO	  + VWO		+ Uni
replace edu16 = 3.5 if onderw1s16==3 & onderw2s16== 4 & onderw3s16==10	& onderw4s16==.			// HAVO	  + VWO		+ Anders
replace edu16 = 3   if onderw1s16==3 & onderw2s16== 6 & onderw3s16== 8	& onderw4s16==.			// HAVO	  + MBO1	+ HBO
replace edu16 = 3   if onderw1s16==3 & onderw2s16== 7 & onderw3s16==10	& onderw4s16==.			// HAVO	  + MBO2-4	+ Anders
replace edu16 = 3   if onderw1s16==3 & onderw2s16== 8 & onderw3s16== 9	& onderw4s16==.			// HAVO	  + HBO		+ Uni
replace edu16 = 4   if onderw1s16==4 & onderw2s16== 5 & onderw3s16== 7	& onderw4s16==.			// VWO	  +	Gym		+ MBO2-4
replace edu16 = 4   if onderw1s16==4 & onderw2s16== 5 & onderw3s16== 8	& onderw4s16==.			// VWO	  +	Gym		+ HBO
replace edu16 = 4   if onderw1s16==4 & onderw2s16== 5 & onderw3s16== 9	& onderw4s16==.			// VWO	  +	Gym		+ Uni
replace edu16 = 4   if onderw1s16==4 & onderw2s16== 5 & onderw3s16==10	& onderw4s16==.			// VWO	  +	Gym		+ Anders
replace edu16 = 4   if onderw1s16==4 & onderw2s16== 8 & onderw3s16== 9	& onderw4s16==.			// VWO	  +	HBO		+ Uni
replace edu16 = 4   if onderw1s16==5 & onderw2s16== 9 & onderw3s16==10	& onderw4s16==.			// Gym	  +	Uni		+ Anders


* 	Mentioned four levels
replace edu16 = 1.5 if onderw1s16==1 & onderw2s16== 2 & onderw3s16== 3 & onderw4s16== 6			// VMBO   + VMBO-t	+ HAVO 	 + MBO1
replace edu16 = 2   if onderw1s16==2 & onderw2s16== 3 & onderw3s16== 4 & onderw4s16== 6			// VMBO-t + HAVO	+ VWO 	 + MBO1
replace edu16 = 2.5 if onderw1s16==2 & onderw2s16== 3 & onderw3s16== 4 & onderw4s16== 7			// VMBO-t + HAVO	+ VWO 	 + MBO2-4
replace edu16 = 3   if onderw1s16==2 & onderw2s16== 3 & onderw3s16== 4 & onderw4s16==10			// VMBO-t + HAVO	+ VWO 	 + Anders
replace edu16 = 3   if onderw1s16==3 & onderw2s16== 6 & onderw3s16== 7 & onderw4s16== 8			// HAVO	  + MBO1 	+ MBO2-4 + HBO

* 	Mentioned five levels
replace edu16 = 2   if onderw1s16==1 & onderw2s16== 2 & onderw3s16== 3	& onderw4s16== 6 & onderw5s16==  7	// VMBO + VMBO-t + HAVO + MBO1 + MBO2-4
replace edu16 = 3   if onderw1s16==1 & onderw2s16== 3 & onderw3s16== 4	& onderw4s16== 5 & onderw5s16== 10	// VMBO + HAVO   + VWO  + Gym  + Anders

codebook edu16
codebook edu16 if edu12<.
codebook edu16 if edu12<. & edu14<.


//	STEP 1.6 --------------------------------------------------------------------------------------------------------
//	Educational attainment age 18 filled in by student

** 	Explore which categories there are over the years

*  	Which survey is taken when?
tab invjaars18 if onderwijss18  <.
*tab invjaars18 if onderw1s18    <.
*tab invjaars18 if onderw2s18    <.
*tab invjaars18 if onderw3s18    <.
*tab invjaars18 if onderw4s18    <.
*tab invjaars18 if onderw5s18    <.
*tab invjaars18 if onderwijsNs18 <.

*	What years of birth do each survey cover?
tab yob if onderwijss18  <. & twzyg<.
*tab yob if onderw1s18    <. & twzyg<.
*tab yob if onderwijsNs18 <. & twzyg<.

tab yob if onderwijss18  <. 
*tab yob if onderw1s18    <. 
*tab yob if onderwijsNs18 <. 

*   Which categories are there in each survey?
tab onderwijss18 	if twzyg<.
*tab onderw1s18		if twzyg<.
*tab onderw2s18		if twzyg<.
*tab onderw3s18		if twzyg<.
*tab onderw4s18		if twzyg<.
*tab onderw5s18		if twzyg<.
*tab onderwijsNs18	if twzyg<.

tab onderwijss18
*tab onderw1s18
*tab onderw2s18
*tab onderw3s18
*tab onderw4s18
*tab onderw5s18
*tab onderwijsNs18


/*
Variabele		Invuljaar	Geboortejaar

onderwijss18	2004-2007	1986-1988 (1981-1994)	2004		DHBQ_14-16-18_versie1_rood-grijsBoekje			VMBO/HAVO/VWO-GYM
													2005		DHBQ_14-16-18_versie2_blauwBoekje_druk1
???
onderw1s18		2009-2014	1991-1996 (1983-1999)	2009		DHBQ14_versie3_ONLINE_2009						VMBO/VMBO-t/HAVO/VWO/Gymnasium 
onderwXs18			""								2011		DHBQ14_versie4_ONLINE_2011						+ 	MBO1/MBO2-4/HBO/UNI/Anders, nl
													2012		DHBQ14_versie5_ONLINE_2012							(meerdere)
onderwijsNs18	2013		1995-1997 (1994-1999)	2003		DHBQ14_versie6_ONLINE_2013 						VMBO/VMBO-t/HAVO/VWO/Gymnasium 
*/


** 	Check whether there are strange things going on

*  	There is no child that has information entered on two occassions  
*list PID if onderwijss18<. & (onderw1s18<. | onderwijsNs18<.)
*list PID if onderw1s18<. & onderwijsNs18<.		

/*
**	Check which combination of levels occur when more than one option is possible

* 	First those who give 2 answers, N=431/(510)
tab onderw1s18 onderw2s18 if onderw2s18>-1 & onderw3s18==-1
tab onderw1s18 onderw2s18 if onderw2s18>-1 & onderw3s18==-1 & twzyg<.
* 	(I've checked: when people have -1 on onderw(X)s18, they have -1 on onderw(X+1)s18.

/*
 1 :	2	3		    6	7	8		10
 2 : 		3	4	    6	7	8		10
 3 :			4	5	6	7   8  	
 4 :				5		   (8)  9	10
 5 :							    9	10
 6 :				        7
 7 :						   (8)	    10	

Between brackets if siblings are included
*/

* 	Those who give 3 answers, N=23/(42)
tab onderw2s18 onderw3s18 if onderw1s18==1  & onderw3s18>-1 & onderw4s18==-1
tab onderw2s18 onderw3s18 if onderw1s18==1  & onderw3s18>-1 & onderw4s18==-1 & twzyg<.
tab onderw2s18 onderw3s18 if onderw1s18==2  & onderw3s18>-1 & onderw4s18==-1
tab onderw2s18 onderw3s18 if onderw1s18==2  & onderw3s18>-1 & onderw4s18==-1 & twzyg<.
tab onderw2s18 onderw3s18 if onderw1s18==3  & onderw3s18>-1 & onderw4s18==-1
tab onderw2s18 onderw3s18 if onderw1s18==3  & onderw3s18>-1 & onderw4s18==-1 & twzyg<.
tab onderw2s18 onderw3s18 if onderw1s18==4  & onderw3s18>-1 & onderw4s18==-1
tab onderw2s18 onderw3s18 if onderw1s18==4  & onderw3s18>-1 & onderw4s18==-1 & twzyg<.
tab onderw2s18 onderw3s18 if onderw1s18==5  & onderw3s18>-1 & onderw4s18==-1
tab onderw2s18 onderw3s18 if onderw1s18==5  & onderw3s18>-1 & onderw4s18==-1 & twzyg<.
tab onderw2s18 onderw3s18 if onderw1s18==6  & onderw3s18>-1 & onderw4s18==-1				
tab onderw2s18 onderw3s18 if onderw1s18==6  & onderw3s18>-1 & onderw4s18==-1 & twzyg<. 		
tab onderw2s18 onderw3s18 if onderw1s18==7  & onderw3s18>-1 & onderw4s18==-1
tab onderw2s18 onderw3s18 if onderw1s18==7  & onderw3s18>-1 & onderw4s18==-1 & twzyg<.
tab onderw2s18 onderw3s18 if onderw1s18==8  & onderw3s18>-1 & onderw4s18==-1
tab onderw2s18 onderw3s18 if onderw1s18==8  & onderw3s18>-1 & onderw4s18==-1 & twzyg<.
tab onderw2s18 onderw3s18 if onderw1s18==9  & onderw3s18>-1 & onderw4s18==-1
tab onderw2s18 onderw3s18 if onderw1s18==9  & onderw3s18>-1 & onderw4s18==-1 & twzyg<.
tab onderw2s18 onderw3s18 if onderw1s18==10 & onderw3s18>-1 & onderw4s18==-1
tab onderw2s18 onderw3s18 if onderw1s18==10 & onderw3s18>-1 & onderw4s18==-1 & twzyg<.


/*
 1 :  2+[3   6 7]   3+[     (8)10]    								6+[7]	  (7+[ 8	 ])    
 2 :  			    3+[4 6 7 8 10]									6+[7]	   7+[(8)(10)]		
 3 :								4+[5(8)]						6+[8] 	   7+[    10 ]	(8+[9])
 4 :											5+[(7)(8) 10 ]  							(8+[9])
 5 :																 									(9+[10])
 6 :		
 7 :		
 8 :
 9 :
 10:
 
Between brackets if siblings are included
*/

*	Those who give 4 answers, N=2/(4)
list onderw1s18 onderw2s18 onderw3s18 onderw4s18 			if onderw4s18>-1 & onderw5s18==-1 
list onderw1s18 onderw2s18 onderw3s18 onderw4s18 			if onderw4s18>-1 & onderw5s18==-1 & twzyg<.
/*
 1+2+3   +6
   2+3+4 +6 
*/

list onderw1s18 onderw2s18 onderw3s18 onderw4s18 onderw5s18 if onderw5s18>-1 & onderw5s18<.
*	Those who give 5 answers, N=1/(1)
* 1+2+3+6+7
*/

** 	Create the actual variable 

*	Mentioned one particular level
gen edu18 =.
replace edu18 = 1.5 if onderwijss18==1															// VMBO (enige optie)
replace edu18 = 3   if onderwijss18==2 	 														// HAVO
replace edu18 = 4   if onderwijss18==3 															// VWO-Gymnasium
** 	More than one level mentioned not available - check whether it's really not available?

codebook edu18
codebook edu18 if edu12<.
codebook edu18 if edu12<. & edu14<.


gen edu_dif_14_12 		= edu14-edu12
gen edu_difabs_14_12	= abs(edu_dif_14_12)
gen edu_dif_16_14 		= edu16-edu14
gen edu_difabs_16_14	= abs(edu_dif_16_14)

tab edu_dif_14_12
tab edu_difabs_14_12
tab edu_dif_16_14
tab edu_difabs_16_14

tab edu_dif_14_12		if twzyg<.
tab edu_difabs_14_12	if twzyg<.
tab edu_dif_16_14		if twzyg<.
tab edu_difabs_16_14	if twzyg<.

egen edu14_16 = rowmean(edu14 edu16)
tab edu14 if edu16<.
tab edu14_16
codebook edu16 if edu12<.
codebook edu14_16 if edu12<.
tab edu16 if edu12<.
tab edu14_16 if edu12<.

tab o9h0_1 
tab oplh0_2 
tab o7hg1 
tab op9h2
tab op9h3
tab o13h4
tab o13h5 
tab o13h6 
tab oplhu7 
tab o9h8 
tab o9h10 


//	STEP 1.7 --------------------------------------------------------------------------------------------------------
//	Educational attainment: 

tab diplvo1s18 diplvo2s18 	// Checked: diplvo2s18>diplvo1s18 if diplvo1s18==0,1,2,3 & diplvo2s18==1,2,3
tab diplvo2s18 diplvo3s18	// Checked: diplvo3s18>diplvo2s18 if diplvo2s18==0,1,2,3 & diplvo3s18==1,2,3

gen dip18 = .
replace dip18 = 1.5 if diplvo1s18==1	
replace dip18 = 3   if diplvo1s18==2	
replace dip18 = 4   if diplvo1s18==3	
replace dip18 = 1.5 if diplvo2s18==1  	 
replace dip18 = 3   if diplvo2s18==2	
replace dip18 = 4   if diplvo2s18==3	
replace dip18 = 1.5 if diplvo3s18==1  	 
replace dip18 = 3   if diplvo3s18==2	
replace dip18 = 4   if diplvo3s18==3	

tab diplvo1s16 diplvo2s16 	// Checked: diplvo2s16>diplvo1s16 if diplvo1s16==0,1,2,3 & diplvo2s16==1,2,3
tab diplvo2s16 diplvo3s16	// Checked: diplvo3s16>diplvo2s16 if diplvo2s16==0,1,2,3 & diplvo3s16==1,2,3

gen dip16 = .
replace dip16 = 1.5 if diplvo1s16==1	
replace dip16 = 3   if diplvo1s16==2	
replace dip16 = 4   if diplvo1s16==3	
replace dip16 = 1.5 if diplvo2s16==1  	 
replace dip16 = 3   if diplvo2s16==2	
replace dip16 = 4   if diplvo2s16==3	
replace dip16 = 1.5 if diplvo3s16==1  	 
replace dip16 = 3   if diplvo3s16==2	
replace dip16 = 4   if diplvo3s16==3	

*	The few cases where dip16>edu16 if edu16<. are not very realistic given the age, so I won't
* 	do anything with them. I will use dip16 only when edu16==. 
list onderwijss14 onderwijss16 lfts16 diplvo1s16 diplvo2s16 diplvo3s16 edu16 dip16 ///
	 if dip16<. & (dip16>edu16 | edu16==.)				

gen edu = . 
replace edu = edu18 if edu18<.
replace edu = dip18 if dip18<. & (dip18>edu18 | edu18==.) 
replace edu = edu16 if edu18==. & dip18==. & edu16<.			// "set to missing" are missings changed to specific missings
replace edu = dip16 if edu18==. & dip18==. & dip16<.  & (edu16==. | edu16 ==.a)
replace edu = edu14 if edu18==. & dip18==. & dip16==. & (edu16==. | edu16 ==.a) & edu14<.
tab edu
codebook edu

* 	Create age at which education is being measured.
gen age = .
replace age = lfts18 if edu18<.
replace age = lfts18 if dip18<. & (dip18>edu18 | edu18==.) 
replace age = lfts16 if edu18==. & dip18==.	& edu16<.			// "set to missing" are missings changed to specific missings
replace age = lfts16 if edu18==. & dip18==. & dip16<.  & (edu16==. | edu16 ==.a)
replace age = lfts14 if edu18==. & dip18==. & dip16==. & (edu16==. | edu16 ==.a) & edu14<.
replace age = . if edu==.
list PID edu14 lfts14 edu16 lfts16 edu18 lfts18 edu age if edu<. & age==. // 4 children have edu16 but lfts16 missing
codebook lfts16
replace age = 16.67 if edu16<. & lfts16==.	// Give 4 children with age missing the average age 
tab age
codebook age

*	Create an alternative measure of educational attainment to get better prediction of unclear "VMBO" category
gen edu_alt = edu									
replace edu_alt = edu12 if edu_alt==1.5 & edu12<.	// Use information at age 12 to make a better prediction


//	STEP 2 ==========================================================================================================
//	Definite vs delayed tracking

* Homogenous vs heterogenous class, based on mother
gen delay_m = .

labelbook hschm12
replace delay_m = 0 if hschm12>= 3 & hschm12<=10									// SVO not included
replace delay_m = 1 if hschm12>=11 & hschm12<=17		

labelbook hschm12_v2		
replace delay_m = 0 if hschm12_v2>= 3 & hschm12_v2<=8								// SVO not included
replace delay_m = 1 if hschm12_v2== 5												// brugklas (anders dan beroepsonderwijs)
	
labelbook hschm12_v3			
replace delay_m = 0 if hschm12_v3>= 3 & hschm12_v3<=11								// SVO not included
replace delay_m = 1 if hschm12_v3>= 6 & hschm12_v3<= 8								// brugklas (anders dan beroepsonderwijs); MAVO/HAVO; HAVO/VWO

labelbook hschm12_1n hschm12_2n 																	
replace delay_m = 0 if hschm12_1n>= 2 & hschm12_1n<=9								// SVO not included
replace delay_m = 1 if hschm12_2n>= 2 & hschm12_2n<=9 & hschm12_1n!=1				// See earlier check that hschm12_2n/_3n are missing if _1n missing

codebook delay_m
tab delay_m

* Homogenous vs heterogenous class, based on father
gen delay_v = .

labelbook hschv12
replace delay_v = 0 if hschv12>= 3 & hschv12<=10									// SVO not included
replace delay_v = 1 if hschv12>=11 & hschv12<=16		

labelbook hschv12_v2		
replace delay_v = 0 if hschv12_v2>= 3 & hschv12_v2<=8								// SVO not included
replace delay_v = 1 if hschv12_v2== 5												// brugklas (anders dan beroepsonderwijs)
	
labelbook hschv12_v3			
replace delay_v = 0 if hschv12_v3>= 3 & hschv12_v3<=11								// SVO not included
replace delay_v = 1 if hschv12_v3>= 6 & hschv12_v3<= 8								// brugklas (anders dan beroepsonderwijs); MAVO/HAVO; HAVO/VWO

labelbook hschv12_1n // hschv12_2n is still missing in the data																
replace delay_v = 0 if hschv12_1n>= 2 & hschv12_1n<=9								// SVO not included
*replace delay_v = 1 if hschv12_2n>= 2 & hschv12_2n<=9 & hschv12_1n>= 2				// See earlier check that hschm12_2n/_3n are missing if _1n missing
  
codebook delay_v
tab delay_v

gen delay = delay_m
replace delay = delay_v if delay_m==.
replace delay = .b if edu12==.b
replace delay = .c if edu12==.c
tab delay
codebook delay


//	STEP 3 ==========================================================================================================
//	Achievement level dummies

* Based on cito-scores (see stcrt-2014-36955.pdf)
recode cito_final (500/523 = 1 "VMBO Basis") 						///
				  (524/528 = 2 "VMBO Kader") 						///
				  (529/536 = 3 "VMBO Theoretisch & Gemengd") 		///
				  (537/544 = 4 "HAVO")								///
				  (545/550 = 5 "VWO"), gen(achie5)
tab achie5, gen(achie5_)
tab achie5_1	
tab achie5_2	
tab achie5_3	
tab achie5_4	
tab achie5_5


* Based on cito-scores extended (see, e.g., stcrt-2016-10096.pdf)
recode cito_final (500/518 = 1 "VMBO Basis") 						///
				  (519/525 = 2 "VMBO Basis / Kader") 				///
				  (526/528 = 3 "VMBO Kader") 						///
				  (529/532 = 4 "VMBO Theoretisch & Gemengd") 		///
				  (533/536 = 5 "VMBO Theoretisch & Gemengd / HAVO") ///
				  (537/539 = 6 "HAVO")								///
				  (540/544 = 7 "HAVO / VWO")						///
				  (545/550 = 8 "VWO"), gen(achie8)
tab achie8, gen(achie8_)
tab achie8_1	
tab achie8_2	
tab achie8_3	
tab achie8_4	
tab achie8_5
tab achie8_6
tab achie8_7
tab achie8_8

* 	Based on cito-scores extended (in-between form)
recode cito_final (500/528 = 1 "VMBO Basis + Kader") 				///
				  (529/532 = 2 "VMBO Theoretisch & Gemengd") 		///
				  (533/536 = 3 "VMBO Theoretisch & Gemengd / HAVO") ///
				  (537/539 = 4 "HAVO")								///
				  (540/544 = 5 "HAVO / VWO")						///
				  (545/550 = 6 "VWO"), gen(achie6)
tab achie6, gen(achie6_)
tab achie6_1	
tab achie6_2	
tab achie6_3	
tab achie6_4	
tab achie6_5
tab achie6_6

* Based on cito-scores collapsed
recode cito_final (500/528 = 1 "VMBO Basis + Kader") 				///
				  (529/536 = 2 "VMBO Theoretisch & Gemengd") 		///
				  (537/544 = 3 "HAVO")								///
				  (545/550 = 4 "VWO"), gen(achie4)
tab achie4, gen(achie4_)
tab achie4_1	
tab achie4_2	
tab achie4_3	
tab achie4_4	


//	STEP 4 ==========================================================================================================
//	Educational attainment parents

**	Age 7 - education mother
tab oplm7 oplmoe7, mis  // missing oplm also have missing on oplmoe
fre oplm7 oplmoe7		// use oplm (more categories and slightly more valid cases

tab invjrm12 if oplm7  <.			// 2002-2018
tab yob 	 if oplm7  <.			// 1991-2006
fre oplm7							// 13 categories

**	Age 10 - education mother
tab oplm10 oplmoe10, mis  // missing oplm also have missing on oplmoe
fre oplm10 oplmoe10		  // use oplm (more categories and slightly more valid cases

tab invjrm12 if oplm10  <.			// 2002-2018
tab yob 	 if oplm10  <.			// 1991-2006
fre oplm10

tab oplm7 oplm10, mis
* 11,523 missings on oplm7  of which 3356 have a valid answer on oplm10
* 12,799 missings on oplm10 of which 4632 have a valid answer on oplm7

**	Age 7 - education father	
tab oplv7 oplvad7, mis  // missing oplv also have missing on oplvad
fre oplv7 oplvad7		// use oplv (more categories and slightly more valid cases

tab invjrm12 if oplv7  <.			// 2002-2018
tab yob 	 if oplv7  <.			// 1991-2006
fre oplv7							// 13 categories

**	Age 10 - education father
tab oplv10 oplvad10, mis  // missing oplv also have missing on oplvad
fre oplv10 oplvad10		  // use oplv (more categories and slightly more valid cases

tab invjrm12 if oplv10  <.			// 2002-2018
tab yob 	 if oplv10  <.			// 1991-2006
fre oplv10

tab oplv7 oplv10, mis
* 12,360 missings on oplv7  of which 3505 have a valid answer on oplv10
* 13,488 missings on oplv10 of which 4633 have a valid answer on oplv7

* 	Education father and mother
tab oplm7  oplv7,  mis				// 11,308 missing for both, 1052 and 215 only one parent
tab oplm10 oplv10, mis				// 12,595 missing for both,  893 and 204 only one parent


sum oplm7 oplv7 oplm10 oplv10 
corr oplm7 oplv7 oplm10 oplv10 
* use age 10, if missing age 7 (correlate highly over time) 

gen oplm = oplm10
	replace oplm = oplm7 if oplm10 ==.
gen oplv = oplv10
	replace oplv = oplv7 if oplv10 ==.

*	Parental education (highest in family)
gen eduparents = oplv
	replace eduparents = oplm if oplv<oplm & oplm!=.
replace eduparents =.a if eduparents == -1

recode eduparents (1/5 = 1 "Primary/lower secondary: basis, mulo, mavo, lts") ///
				  (6/9 = 2 "Upper secondary: havo, vwo, mbo") ///
				  (10/11=3 "Tertiary: hbo, few years hbo/uni") ///
				  (12/13=4 "Tertiary: uni, postdoctoral"), gen(edup4c)
tab edup4c


//	STEP 5 ==========================================================================================================
//	Generate several IDs

sort FamilyNumber Extension

*	Generate a dummy indicating whether a twin pair is MZ (1) or DZ (0)
gen mz = 1 if (twzyg==1 | twzyg==3)
replace mz = 0 if (twzyg==2 | twzyg==4 | twzyg==5 | twzyg==6)
tab mz

* 	Give "sex" a more informative name
tab sex
labelbook sex				// males coded as 1; females as 2
gen male = .				// give missing value to 4 people indicated as gender neutral
replace male = 1 if sex==1
replace male = 0 if sex==2
tab sex male

*	Give within a family two children that are part of the same twin pair the same number.
tab Extension if multiple_type==2				// There are multipe twin pairs in some families
gen     same_twin = 0 if multiple_type==2
replace same_twin = 1 if multiple_type==2 & (Extension==1 | Extension ==2)
replace same_twin = 2 if multiple_type==2 & (Extension==3 | Extension ==4)
replace same_twin = 3 if multiple_type==2 & (Extension==5 | Extension ==6)
tab same_twin

* 	Give each twin pair a unique id.
egen twin_id = group(FamilyNumber same_twin)
codebook twin_id

*	Calculate the number of a twin pair that is observed 
by twin_id, sort: egen twinsize1a = count(PID) if twin_id<. 
tab twinsize1a

*	Give within a family three children that are part of the same triplet the same number.
tab Extension if multiple_type==3				// In families with a triplet, there is only one triplet
gen same_trip  = 1 if multiple_type==3
tab same_trip

* 	Give each triplet a unique id.
egen trip_id = group(FamilyNumber same_trip)
codebook trip_id

*	Calculate the number of a triplet that is observed 
by trip_id, sort: egen tripsize = count(PID) if trip_id<. 
tab tripsize
gen     twin_trip = 0 if multiple_type==3 
replace twin_trip = 1 if multiple_type==3 & tripsize==2
replace twin_trip = 1 if multiple_type==3 & tripsize==3 & mult_ext <= 2	// They lack info on twzyg, so not included this way
tab twin_trip 

* 	Identify each complete twin pair and the first complete pair within a triplet
gen pair_complete = .
replace pair_complete = 1 if twinsize1a==2 | twin_trip==1
tab pair_complete

*	Give temporary id to each twin pair and to the first complete pair within a triplet
gen 	 twtr_id = twin_id 			 if multiple_type==2
replace  twtr_id = trip_id + 1000000 if multiple_type==3
codebook twtr_id

* 	Give id to each complete twin pair and to the first complete pair within a triplet
egen pair_id = group(twtr_id pair_complete) 
codebook pair_id

* 	Fix a wrong year of birth
sort FamilyNumber Extension
list FamilyNumber if FamilyNumber==FamilyNumber[_n+1] & twtr_id==twtr_id[_n+1] & yob!=yob[_n+1] & twtr_id<. // One family: 1995 vs 1996
replace yob=1996 if FamilyNumber==15140	 // Based on ages at filling in other surveys, correct yob seems to be 1996  


//	STEP 6 ==========================================================================================================
//	Select cases to be analyzed

* 	Select only twins (no triplets or higher included nor siblings)
keep if multiple_type==2					
tab twinsize1a
codebook twin_id 

*	Select birth cohorts for which information on cito and edu was collected
keep if yob>1985 & yob<2000
tab twinsize1a
codebook twin_id 

* 	Select only twins that entered high school already at the time of survey
keep if edu12 != .c
tab twinsize1a
codebook twin_id 

* 	Select only twins that are in a regular education (not in special education)
keep if edu12 != .b
tab twinsize1a
codebook twin_id  

*	Show the number of complete / incomplete pairs after these selections
by twin_id, sort: egen twinsize1b = count(PID) 
tab twinsize1b
codebook twin_id if twinsize1b==2

*	Calculate the number of complete pairs with twzyg and male not missing
by twin_id, sort: egen twinsize2a = count(PID) if twzyg<. & male<. 
tab twinsize2a
codebook twin_id if twinsize2a==2

*	Create variable that helps identify analytic sample
*	(=only those pairs where at least one has non-missing for cito and only persons that have non-missing for cito)
by twin_id, sort: egen cito_nm = count(cito_final) 
tab cito_nm

*	Calculate the number of (in)complete pairs of analytical sample 
by twin_id, sort: egen twinsize2b = count(PID) if twzyg<. & male<. & cito_final<. & cito_nm>0
tab twinsize2b
codebook twinsize2b
codebook twin_id if twinsize2b<.
codebook twin_id if twinsize2b==2	// number of complete pairs for twzyg, male, cito_final

*	Calculate the number of a twin pair that is observed with none of the analyzed variables missing
by twin_id, sort: egen twinsize3 = count(PID) if twin_id<. & twzyg<. & male<. & cito_final<. ///
														   & delay<. & edu<.
tab twinsize3
codebook twin_id if twinsize3==2

*	Select only twins for which zygosity is known
keep if twzyg<.
tab yob 
codebook twin_id 

*keep if pair_id < . 						// Only complete twin pairs and the first complete pair within a triplet
tab Extension
tab mult_ext


//	STEP 7 ==========================================================================================================
//	Make openmx ready: give pseudo-missings, center variables, etc.

*	Center variables
center age 		  if male<. & cito_final<. & cito_nm>0  
center cito_final if male<. & cito_final<. & cito_nm>0, gen(c_cito)
center cito_final if male<. & cito_final<. & cito_nm>0, gen(z_cito) s
sum yob , detail
gen c_yob = yob - 1992 // subtract the mean / median integer (subtracting exact mean leads to rounding problems)

/* 	Create definition variable for cito and give all definition variables pseudo-missing codes.
	Openmx is not able to analyze cases with missing info on definition variables. 
	Sometimes one twin has missing info but the co-twin does not, and it would be a waste to kick
	both out. A trick to prevent this in openmx is to give a pseudo-missing value (e.g., -999) to 
	the twin that has missing on the definiton variable but set the dependent variables to actual
	missing for them to prevent that they end up in the analyses with the weird pseudo-missing value.
	Cito is in this case both a dependent variable and a definition variable (when it functions as moderator)
*/
gen     defm   = c_cito
replace defm   = -999 if c_cito==.
replace edu    = .    if c_cito==.	// This makes sure those with missings are not actually analyzed
replace edu    = .    if   male==.	// This makes sure those with missings are not actually analyzed
replace c_cito = .    if   male==.	// This makes sure those with missings are not actually analyzed
replace male   = -999 if   male==.

*	Generate a new "extension" variable based on order what children respond for variable delay
sort twin_id delay mult_ext
bysort twin_id: gen mult_ext2 = _n
codebook mult_ext2

* 	Select variables to be kept
keep PID FamilyNumber pair_id twin_id Extension multiple_type mult_ext mult_ext2 mz twzyg male yob c_yob 		///
	 allochtoon nkids agem12 agev12 lfts14 lfts16 lfts18 age c_age edup4c eduparents oplv oplm 			 		///
	 edu12 edu14 edu16 edu18 edu edu_alt delay c_cito z_cito cito_final achie4 achie5 achie6 achie8      		///
	 achie4_1-achie4_4 achie5_1-achie5_5 achie6_1-achie6_6 achie8_1-achie8_8 twinsize1b twinsize2b defm cito_nm


//	STEP FINALIZE ===================================================================================================
//	Metadata

*	!!Update this info if you make a change!!
label data "NTR_P2_1_`VersionSF1' 2021-02-02"
notes: VENI Project 2 VENI_P2_clean`VersionSF1' \ VENI_P2_1_Clean06.do \ 		///
NTR data cleaned for P2, Version `VersionSF1' \ a.knigge@uu.nl \ 2021-02-02
datasignature set, reset

*	!!If you make a change to the do-file, make a short entry of the change in the log below!!
/*
	Date		Author			Description of change
	2019-07-31	a.knigge@uu.nl	Based on P1_1_Clean01.do; Added edu16 and edu18; fixed mistakes in edu12; 
								changed creation of delay
	2019-09-26	a.knigge@uu.nl	Did not delete pairs who have sex missing as was the case previously.	
	2020-01-16	a.knigge@uu.nl	Based order of twin pairs on 0: immediate; 1: delayed; .: missing info,
								such that the order of disconcordant twins is always 0-1.
								Renamed "ability" dummies into "achievement"
								Generated a dummy indicating whether pair is MZ or DZ.
								Selected only twins born between 1986 and 1999 with known zygosity.
								Renamed "hetero_class" into "delay"
	2020-09-07	a.knigge@uu.nl	Changed scoring of educational attainment (VMBO-g equal score as VMBO-t);
								Changed "achievement"-dummies slightly based on the official staatscourant table;
								Calculated centering based on cases actually analyzed;
								Kept more variables in for descriptive analyses.
	2021-01-26	a.knigge@uu.nl	Kept those pairs in where one has missing on cito but the other has not.
								Gave "pseudo-missing codes" to missing definition variables data.
*/

*	Save the data
save "`path1'NTR_P2_1_`VersionSF1'", replace

log close
