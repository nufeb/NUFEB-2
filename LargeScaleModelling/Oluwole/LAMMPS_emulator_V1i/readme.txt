############LAMMPS EMULATOR version V1i####################
(A) DESCRIPTION
A coupled-emulator which produces rates of change of nutrient concentrations namely 
(Sub,O2,NO2,NO3,NH4) and rates of change of the following characterized outputs
(rate of change of nutrient concentration and Biomass for each species) respectively.
For instance, dS/dt=(S_2-S_1)/(T_2-T_1)

#9 outputs
(i) The rate of change of SUB/Carbon substrate concentration (kg/m^3)
(ii) The rate of change of O2 concentration (kg/m^3)
(iii) The rate of change of NO2 concentration (kg/m^3)
(iv) The rate of change of NO3 concentration (kg/m^3)
(v) The rate of change of NH4 concentration (kg/m^3)
(vi) The rate of change of biomass for HET (kg)
(vii) The rate of change of biomass for AOB (kg)
(viii) The rate of change of biomass for NOB (kg)
(ix) The rate of change of biomass for EPS (kg)

(B) All required data and C++ scripts for the emulator are in the directory "LAMMPS_emulator_V1i". 
The main c++ source code is "emu.cpp". 

(C) Input file description:
A sample "input.txt" file is placed in the "input" subfolder which is a required input parameter for running the emulator.
This contain 9 different values
The input values [1-5] represent "SUB","O2","NO2","NO3","NH4" 
local nutrient substrates (kg/m^3) respectively. 
The input values [6-9] represent current values of Biomass 
concentration for each species.

(D) RUN the code:
The main source code is "emu.cpp". The emulator can be run from either the Window or Linux platform. 
Put the input file "input.txt" to be run in the subfolder "input".
 These files contain the parameter values and initial conditions.
To run any case, the input parameters must be converted to the 
units stated above to get correct results.

To run from Linux terminal:
$ cd path to /LAMMPS_emulator_V1i
/** compile .cpp file
$ g++ -O3 -I/usr/include/eigen3/ emu.cpp -o emu
/** run the code as
$ ./emu

(E) Obtain the simulation results from subfolder "output".
(F) The valid ranges of the nutrient concentrations (kg/m^3) and other parameters:
["SUB","O2","NO2","NO3", "NH4","Biomass_HET",
"Biomass_AOB","Biomass_NOB","Biomass_EPS"]  
min=[0.003593248 0.001885759  2.622689e-119 5.293362e-118 4.234222e-116
 1.046317e-14  4.196982e-15    2.840420e-15  2.840420e-15]
max=[0.044765000 0.014969000  8.995006e-03  8.995000e-03  8.995000e-03 
 2.275782e-10  2.159382e-10    4.261543e-14  4.271532e-14]
  
NOTE: This emulator will produce output for any specified time-steps and any parameter combinations. 
CONTACT: oluwole.oyebamiji@ncl.ac.uk or wolemi2@yahoo.com for further details.