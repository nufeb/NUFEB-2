##Configuration options for the LPJML emulator
# select "j" from options (j=c(0,1,2,3,4,5,6,7,8)
# select the "rcp" by uncomment the line
# select the "coupling" type  by uncomment the line
# select CO2 fertilization level (ON or OFF)
# select "mannagement type" from options (man=c(1,2,3,4,5,6,7)
# select "crop types" from options (irrigated or rainfed)
# NOTE: j==0 implies output all the 8 decadal change once
#Requred: ClimGEN netcdf climate data to run the emulator in the format below 
##############indicate path to climGEN climate data
path="E:\\clim\\management_irrigated_rainfed_co2_on_off\\had2"

j=8   ###choose time point

out="current"
#out="future"
###########choose RCP
#rcp="RCP3PD"
#rcp="RCP45"
rcp="RCP6"
#rcp="RCP85"
#rcp="other"
###########################CO2 fertilization
fertilization=1   ###runif(n,0,1)## any number betwen 0 to 1
###0==OFF; 1==ON
###########################management levels
#man=1
#man=2
#man=3
#man=4
man=5
#man=6
#man=7
##########################choose crop types
type="irrigated"
#type="rainfed"
###############################END