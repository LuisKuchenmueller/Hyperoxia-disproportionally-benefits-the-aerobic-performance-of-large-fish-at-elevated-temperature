README description of files Kuchenmueller 2024 

Data and scripts in the repository relate to the following research experiment:
Hyperoxia disproportionally benefits the aerobic performance of large fish at elevated temperature
Luis L. Kuchenmueller*, Elizabeth Hoots, Timothy D. Clark
*corresponding author: l.kuchenmueller@deakin.edu.au

Please see the manuscript for a description of the experimental protocol and the statistical approach.

Please note that script '1_rainbow_import_processing.R' has to be run first.
Further, the working directory has to be set to the environment of the data files 'SMR.xlsx', 'MMR.xlsx', and 'Mortalities.xlsx'.

1_rainbow_import_processing.R
	Import respirometry data, calculate metabolic rates, export metabolic data
2_rainbow_statistics.R
	Confirm model assumptions, set up linear model for each metabolic parameters, extract scaling factors and scaling exponents for each metabolic parameter, combine into full data table
3_rainbow_allometry_graphs.R
	Produce Figure 2 showing standard metabolic rate (SMR), maximum metabolic rate (MMR), and absolute aerobic scope (AAS) for all treatments across the mass range on log-log axes, note axes titles and legend are added subsequently
4_rainbow_mass_standardised_graphs.R
	Set up treatment-specific power functions, mass-standardise data, produce Figure 3 showing mass-standardised SMR, MMR, and AAS for all treatments, note axes labels and axes titles are added subsequently
5_rainbow_diagnostic_graphs.R
	All 130 individual respirometry trials are plotted with SMR calculations and background respiration, calculated SMR is compared across trials
6_rainbow_mortalities.R
	Metabolic rates are calculated and respirometry trials are plotted for all 8 mortalities

Raw LabChart data
	(provided for review of how raw data was obtained from oxygen slopes)
	See the manuscript for methods
	The raw data files are named and described as follows:
		1. date in yyyymmdd
		2. oxygen treatment
		3. temperature treatment
		4. fish size class
		5. 'TIMER' if oxygen was samples at 1 Hz causing latency in pump timer to slopes, consequently slopes were extracted per hand instead of macro

SMR.xlsx
	Extracted respirometry data from raw labchart oxygen slopes
	Each row contains one oxygen slope of one "Phase" of one respirometry chamber
		Trial		trial identifier; individuals which were experimented simultaneously share the same number
		ID_chamber	identifies the respirometry tray (A-D) and chamber (1-4); NOTE "volume_ch" has to be referenced to distinguish between respirometer volume
		ID_fish		individual fish identificator (1-130)
		mass		individual mass (g)
		volume_ch	respirometer volume (ml)
		DO_unit		dissolved oxygen unit; in this experiment always (mg O2)
		O2		oxygen treatment; 'normoxia'= 100% air saturation or 'hyperoxia'= 150% air saturation
		Temp_class	temperature treatment class; '17' or '21' or '25' (degrees celsius)  
		dateTime	mean date and time of each "Phase" (dd/mm/yyyy hh:mm:ss)
		Phase		continous numerical identifier of each extracted oxygen slope per 'Trial'
		Temp		mean temperature sensor measurements for each 'Phase' (degrees celsius)
		slope_wBR	mean chamber oxygen slope (mg O2) for each 'Phase'
		BRSlope		linearly regressed and time-matched slope (mg O2) of mean background oxygen slopes before and after each 'Trial'

MMR.xlsx
	Extracted post-chase respirometry data from raw labchart oxygen slopes
	Each row contains one oxygen slope of one individual in one respirometry chamber
		Trial		trial identifier; individuals which were experimented simultaneously share the same number
		ID_chamber	identifies the respirometry tray (A-D) and chamber (1-4); NOTE "volume_ch" has to be referenced to distinguish between respirometer volume
		ID_fish		individual fish identificator (1-130)
		mass		individual mass (g)
		Maturity	individual sexual maturity; 'N' (No) describes gonads of stage 3 and below, 'Y' (Yes) describes gonads greater than stage 3
		volume_ch	respirometer volume (ml)
		DO_unit		dissolved oxygen unit; in this experiment always (mg O2)
		O2		oxygen treatment; 'normoxia'= 100% air saturation or 'hyperoxia'= 150% air saturation
		Temp_class	temperature treatment class; '17' or '21' or '25' (degrees celsius) 
		dateTime	mean date and time of each oxygen slope (dd/mm/yyyy hh:mm:ss)
		Temp		mean temperature sensor measurements for each oxygen slope (degrees celsius)
		slope_wBR	mean chamber oxygen slope (mg O2)
		BRSlope		linearly regressed and time-matched slope (mg O2) of mean background oxygen slopes before and after each 'Trial'

Mortalities.xlsx
	Extracted respirometry data from raw labchart oxygen slopes of mortalities
	Each row contains one oxygen slope of one "Phase" of one respirometry chamber
		ID_chamber	identifies the respirometry tray (A-D) and chamber (1-4); NOTE "volume_ch" has to be referenced to distinguish between respirometer volume
		ID_fish		individual fish identifier (131-138)
		mass		individual mass (g)
		volume_ch	respirometer volume (ml)
		DO_unit		dissolved oxygen unit; in this experiment always (mg O2)
		O2		oxygen treatment; 'normoxia'= 100% air saturation or 'hyperoxia'= 150% air saturation
		Temp_class	temperature treatment class; '17' or '21' or '25' (degrees celsius)  
		dateTime	mean date and time of each "Phase" (dd/mm/yyyy hh:mm:ss)
		ch_time		time passed since the individual was introduced into the chamber (hours)
		Temp		mean temperature sensor measurements for each oxygen slope (degrees celsius)
		slope_wBR	mean chamber oxygen slope (mg O2) for each oxygen slope
		BRSlope		linearly regressed and time-matched slope (mg O2) of mean background oxygen slopes before and after each "Trial"
		
	