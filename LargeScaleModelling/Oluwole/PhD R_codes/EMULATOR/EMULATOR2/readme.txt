LPJmL EMULATOR with version "LPJML_crop-yield_code_v3f" called LPJmL with management_irrigated_rainfed_co2_on_off emulator.
This version has option to choose between rainfed/irrigated. It has option to change CO2 fertilization levels to any desirable values. The output data will either be a "csv" or netcdf depending on what is specified in the configuration options.

########################################################################################################## 
The source code is "LPJML_crop-yield_code_v3f.R".

 All data for the emulator are in the directories rainfed and irrigated folder. This emulator can run for various configuration options (time-step, RCP,CO2 fertilization, crop management levels, coupling type and crop type. Any desire output result has to be set in the configuration file "configuration_option". The "path" to the climate input for running the emulators has to be indicated in the configuration file. All the output data will go into the folder "output_result".
The emulators output the results in 2 different formats which has to be specified in the configuation_option "coupling" which can either be "GEMINI-E3" or "other", if "coupling==other" then the emulator will produce the output in a gridded format on 0.5 by 0.5 degree resolutions. Otherwise, the emulator will futher aggregates the result to country levels as an excel (.csv) file. 

PLEASE NOTE:
Outputs for these emulators are the change in crop-yield relative to the baseline period (2005-2014). The 4 outputs are "tropical cereal", "pulses", "temperate root" and "soybean".

Major options to set in the configuration file for running the emulator
##################################################################
(i)  select "man" from (man=c(1,2,3,4,5,6,7); corresponds to management levels
(ii) select "j" from options (j=c(0,1,2,3,4,5,6,7,8); decadal levels
(iii) select "GCM" by uncomment the line in configuration_file; climate models (THIS OPTION IS NO LONGER APPLICABLE)
(iv) select the "rcp" by uncomment the line in configuration_file; RCP scenarios
(v)  select the "coupling" by uncomment the line in configuration_file; this represent output format
(vi) select the "fertilization" by uncomment the line in configuration_file; CO2 fertilization level
(vii) select "type" by uncomment the line in configuration_file; crop types: rainfed or irrigated

NOTE: this emulator will produce output correspond to the time-step below by selecting "j"
j=1 change in crop yield between (2005-2014) and (2015-2024)
j=2  ............................(2005-2014) and (2025-2034)
j=3  ............................(2005-2014) and (2035-2044)
j=4  ............................(2005-2014) and (2045-2054)
j=5  ............................(2005-2014) and (2055-2064)
j=6  ............................(2005-2014) and (2065-2074)
j=7  ............................(2005-2014) and (2075-2084)
j=8  ............................(2005-2014) and (2085-2094)
NOTE: j==0 implies output all the 8 time points at once



