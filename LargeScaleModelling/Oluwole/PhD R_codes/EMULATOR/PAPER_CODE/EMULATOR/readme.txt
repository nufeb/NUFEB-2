###### LPJmL CROP YIELD EMULATOR FOR THE PAPER "Emulating impacts of ####climate change on crop yields"

DESCRIPTION
This version of LPJmL crop emulator produces crop yield change (gC/m^2) and its baseline values for 5 different crops: Temperate cereal, rice, maize, oil and groundnuts both rainfed and irrigated. All the outputs are decadal values. The emulator provides its outputs as a “HIGH” resolution. The HIGH resolution is a gridded scale of 0.5 by 0.5 degrees. The emulator is driven by ClimGEN climate data. The output file is a Netcdf file.
 
All required data and R scripts for the emulator are in the directories rainfed and irrigated subfolders. All the outputs from running the emulator go into the subfolder “output_results”. The main R source code is "LPJMLEM.R". There are some options to select and set in the configuration file before the emulator can be run to produce desire results.

INPUT DATA:
The emulator code requires four climate data (temperature, precipitation, cloud cover, wetday frequency). The data are provided separately as a Netcdf format high resolution on 0.5 by 0.5 degree resolution from 2001-2100. 


CONFIGURATION FILE
This emulator can run for various configuration options in order to change the default options. There are 5 basic options namely (time-step, RCP, CO2 fertilization, crop management levels and crop types. Any desire output results can be set in the configuration file "configuration_option" before running the main code. 

Major options to set in the configuration file:
(i) Select "j" from options (j=c(0,1,2,3,4,5,6,7,8); decadal levels

(ii) Select "man" from (man=c(1,2,3,4,5,6,7); corresponds to crop management levels. 1 represents poorly managed and 7 represents well-managed option.

(iii) Select the "rcp" by uncomment the line; RCP scenarios. 
NOTE: If option “other” is selected for the RCP. It represents an arbitrary emission scenario and then another 100 year of CO2 data has to be replaced in the “CLIMGEN” subfolder as a “.txt” file format. See the sample data already in the folder.

(iv) Select the "fertilization" by uncomment the line. This option represents the CO2 fertilization level. This option can either be 0 and 1. Emulator assumes there is full (100%) CO2 level if this option is set to equal 1. There is no CO2 fertilization if set equal to 0. It can also be set to any real number between 0 and 1 to indicate intermediate level of CO2.

(v)Select "type" by uncomment the line; this option represents crop types, it can either be a rainfed or irrigated crop.

RUN THE CODE
The main source code is "LPJMLEM.R". The emulator can be run from either the Window or Linux platform. 
 (A) Put the four climate data in the subfolder “CLIMGEN”. These are precipitation, temperature, wetdays frequency and cloud cover.
 (B) Set the configuration file.
 (C) Run the emulator: 
(i) From Window machine:
 At R console………
> setwd(“path to/EMU”)
> source(“LPJMLEM.R”,echo=TRUE)

(ii) To run from Linux terminal:
$ cd pathto/EMU
$ R CMD BATCH LPJMLEM.R result.txt 

(D)Requirements: The emulator requires installation of ncdf, foreach, doParalell and iterator packages in addition to the R standard base packages.


(E)Obtain the results from subfolder “output_results”.

OUTPUT RESULTS
For each run of the emulator, the Netcdf output data produces both the change in yield and its respective baseline value. 
The Netcdf are self-explanatory for each corresponding results.

NOTE: The emulator outputs correspond to the following time-steps for combinations of any option set above for different "j"
j=1 change in crop yield between (2005-2014) and (2015-2024)
j=2  ............................(2005-2014) and (2025-2034)
j=3  ............................(2005-2014) and (2035-2044)
j=4  ............................(2005-2014) and (2045-2054)
j=5  ............................(2005-2014) and (2055-2064)
j=6  ............................(2005-2014) and (2065-2074)
j=7  ............................(2005-2014) and (2075-2084)
j=8  ............................(2005-2014) and (2085-2094)
j==0, ALL 8 decade results at once.


CONTACT: 
oluwole.oyebamiji@open.ac.uk or 
wolemi2@yahoo.com 
for further details.
