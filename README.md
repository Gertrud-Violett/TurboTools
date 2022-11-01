# TurboTools
Turbo Compressor Design Tools
MIT License.
***

## ROSS
ROSS Sample Scripts for Tutorial
Refer to blog for detailed usage
***

## Power Matching
### TransientSingleStage: 
Single Stage Transient Matching Tool
* MAIN_Constant.py: 
Main script for Constant flow matching with multi data and AFR/Bypass mode
syntax: MAIN_Constant.py setupfilename matchingmode (AFR or Bypass)

* MAIN_Transient.py:
Main script for transient flow matching, single data only

* StageCalc.py:
Classes for each turbo component calculation. Execute from MAIN

#### Inputs: 
* cmap_n.csv: compressor map for nth stage
  * NC*: Compressor rpm
  * PRC: Pressure ratio
  * WC*: Flow rate [kg/s]
  * EtaC: Compressor Adiabatic Efficiency
* tmap_n.csv: turbine map for nth stage
  * NT*: Turbine rpm
  * PRT: Pressure ratio
  * WT*: Flow rate [kg/s]
  * EtaT: Turbine Adiabatic Efficiency
* params.toml: Calculation conditions
* target.csv: target file for csv mode multi data batch
* setup.ini: Calculation conditions to be specified by MAIN syntax
  * Input filename: filename for parameters file: default is 'params'
  * Target filename: filename for matching target: default is 'target'
  * Input mode: Select csv for multi data batch or toml for single data
  * Data numbers: must be same as number or rows in target file
  * Compressor Map Name: default is cmap
  * Turbine Map Name: default is tmap

#### Outputs:
* Calculation result Output files
