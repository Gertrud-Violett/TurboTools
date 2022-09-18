# TurboTools
Turbo Compressor Design Tools
MIT License.
***

## ROSS
ROSS Sample Scripts for Tutorial
Refer to blog for detailed usage
***

## Power Matching
Turbo Compressor Matching Tool

### Inputs
Contains input .ini and toml data
#### cmap_n.csv
Compressor Map for nth stage: e.b. cmap_1.csv for 1st stage\

#### tmap_n.csv
Turbine Map for nth stage: e.b. cmap_1.csv for 1st stage\

#### setup.ini
Setting file for whole calculation\

[Settings]\
Initial rpm: Used to set initail rpm. Range must be within cmap min./max\
Analysis Name: any name can be used.\
Input filename: .toml input file name\
Matching Mode (Bypass or AFR): Matching mode, Bypass mode sets bypass turbine flow. AFR automatically sets required turbine inlet temp AFR.\

[Environment Variables]\
humidity [0-1]: use 0 as default, dry air\


#### params_n.toml
Calculation parameters for nth stage: e.g. params_1.toml for 1st stage.
Sample values are shown

[Compressor]\
AmbientPressure = 100 #Inlet Pressure [kPaA] Total Pressure\
AmbientTemperature = 298 #Inlet Temperature in K\
Weight = 0.1 #kg\
Inertia = 0.1 #kgm2\
gamma = 1.395 #gas constant\
Flow = 0.8 #kg/s Flow Target\
PRC = 5.0 #Pressure Ratio Target\
BleedRatio = 0.1 #Ratio of compressor flow bleed\

[Combustor]\
CombustionTemp_T1t = 1473 #Combustion Temp in K\
CombustordeltaP= 0.1 #Ratio Inlet Pres.loss\
AirFuelRatio = 40 #kg/kg\
FuelLatentHeatLHV = 38 #MJ/kg for lower fuel heating value\
TempEfficiency = 1.0 #Actual Temp vs Theoretical Temp\

[Turbine]\
Weight = 0.1 #kg\
Inertia = 0.1 #kgm2\
BackPressure = 110 #Turbine outlet Initial backpressure [kPaA] Static Pressure: Overwritten by calcualtion. Delta can be used a dynamic pressure\
BleedRatio = 0.0 #Ratio of turbine flow bleed\

[Shaft]\
Weight = 0.1 #kg\
Inertia = 0.1 #kgm2\
MechLoss = 0.05 #Mechanical Power Loss %\

[Calculation]\
Step = 100 #rpm resolution for calculation convergence[rpm]\
PowRes = 3 #Target Power Resolution for calculation convergence[kW]\
PRCRes = 0.01 #PRC Resolution for calculation convergence\
WtRes = 0.1 #Wt* Resolution for calculation convergence\

### Outputs
Contains output toml data. Sample values are shown\
#### Result_AFR_Analysis name.toml\
Output result for AFR mode\

#### Result_Vypass_Analysis name.toml
Output result for Bypass mode\
[Compressor]\
"Power[kW]" = 195.002  #Power used by compressor \
"Corrected Flow rate Wc*" = 0.824  #Compressor corrected flow rate\
"Speed RPM" = 49800.0  #rpm at matching speed\
"Outlet Pres[kPaA]" = 499.89  #Compressor outlet pressure\
"Outlet Temp[degC]" = 271.439  #Compressor outlet temperature\
"Pres Ratio PRC" = 4.999  #Compressor pressure ratio\
"Comp Efficiency" = 0.698   #Compressor stage adiabatic efficiency\
"Air Flow[kg/s]" = 0.8   #Actual air flow kg.s\
"Turbine Pressure Ratio P2c/P1t" = 1.111   #Compressor outlet / Turbine Inlet pressure ratio\

[Turbine]\
"Enthalpy[kW]" = 300.959   #Turbine Total Enthalpy, including bleed.\
"Turbine Power[kW]" = 205.266   #Power used for turbine\
"Speed RPM" = 49800.0   #rpm\
"Inlet Pres.[kPa]" = 449.901  #Turbine inlet pressure\
"Inlet Temp[C]" = 962.212   #Turbine Inlet temperature at adjuested AFR\
"Outlet Pres[kPa]Tot" = 160.726   #Turbine outlet pressure\
"Outlet Temp[degC]" = 745.253   #Turbine oulet temperature\
"Pres Ratio PRT" = 2.799   #turbine expansion ratio\
"Turb Efficiency" = 0.713   #Turinbe stage adiabatic efficincy\
"Total Flow[kg/s]" = 0.733  #Turbine stage total flow\
"Turbine Flow[kg/s]" = 0.738  #Flow to turbine wheel
"Bypass Ratio" = 0.008\

["Shaft&Combustor"]\
"Fuel Flow [kg/s]" = 0.013   #Fuel flow after adjusted AFR\
"Turbine Required Flow[kg/s]" = 0.727   #Required flow for turbine stage power\
"Corrected Flow rate[kg/s]" = 0.344  #Corrected turbine flow rate\
"Theoretical T1t[degC]" = 962.212\
"Mechanical Loss [kW]" = 10.263\
AFR = 54.651\

