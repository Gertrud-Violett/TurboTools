# TurboTools
Turbo Tools Public Repo cloned from NHI_TurboTools <br>
Copyright (c) makkiblog.com <br>
MIT License <br>
coding: utf-8 <br>

***

### history
history for outdated code for each tool.<br>


<br><br>
***

## CompressorDesign
Centrifugal Compressor Design Tool
* MAIN.py:
  Main Calculation: <br>
   syntax: MAIN.py parameters_filename Studyname <br>
   output: Result_Studyname.toml <br>
* CentrifugalCompressor.py:<br>
   Centrifugal compressor wheel calculation. Execute from MAIN
* bladegen.py:<br>
   2D blade generation tool. <br>
    syntax: bladegen.py parameters_filename<br>
    output: bladeprofile_Studyname.csv<br>
* RadialDiffuser.py:<br>
   Diffuser design tool (vaned and vaneless). Execute from MAIN
* outletduct.py:<br>
   outlet duct desing tool. Execute from MAIN
#### Input
Input toml and csv file
#### Output
Output png and csv file


<br><br>
***


## Compressor Map Plot
Simple Centrifugal Compressor Mapping Tool
* CompressorMapPlotvCV_v3.ipynb:<br>
  Main Jupyter file to plot compressor map
* CompressorMapPlotCV_v3.py:<br>
  Same python file
#### Data
Input and output datas


<br><br>
***

## DuctDesignTool
Sub-sonic nozzle, turning vane and Diffuser calculation tool <br>
### TurbAnnularDiffuser
Simple 1D Turbine Outlet Annular Diffuser Calculation tool <br>
* MAIN.py: <br>
  Main calculation <br>
  Syntax: MAIN.py parameter_filename geometry_filename
* calc.py: <br>
  Calculation script. Execute from MAIN.py

#### Input
Input toml and csv file
#### Output
Output png and csv file


<br><br>
***



## LifeAnalysis
Low Cycle Fatigue Calculation Tool
* LCF_Calc.py:<br>
Syntax: LCF_Calc.py DamageMap DutyCycle RatedSpeed InputMultiplyer LowPassFilter <br>
Main Calculation script<br>


#### Input
Calcualtion Input files
#### Output
Output files

<br><br>
***


## NoiseAnalysis
Noise Analysis Tools
* Noise_FFT_Analyser.ipynb<br>
  wav file noise analysis tool
### FreqCut
  Plotted result file

<br><br>
***

## PerformanceTestPlotting
Performance Test Plotting Tools<br>
* GS_Data_Analysis.ipynb:<br>
  Gas Stand data analysis tool including thrust and LPT test.

* GeneralPlot.ipynb:<br>
  General Plot tool can be used for x-y axis transposed datas.

### Data
Input and result datas


<br><br>
***

## Power Matching

### EngineMapMaker:
Pre-Processing Tools to create Engine data for ICE<br>
* MAIN_EngineMap.py<br>
Main Script for Engine data creation<br>
syntax: python MAIN_EngineMap.py basemapfile targetmapfile<br>

Following Tools are used:<br>
* T1t_Interpolate.py <br>
Turbine inlet temperature calcualtor. Refer to TurbineTemp_Calculator for full tool.<br>
* MAIN_Constant_Map.py<br>
Matching tool. Refer to TransientSingleStage for full tool <br>  
* CompressorMapPlotCV_v3.py<br>
Compressor Map resolution improvement tool. Refer to CompressorMapPlot for full tool.<br>

### TransientSingleStage: 
Single Stage Transient Matching Tool<br>
* MAIN_Constant.py: <br>
Main script for Constant flow matching with multi data and AFR/Bypass mode<br>
syntax: MAIN_Constant.py setupfilename matchingmode (AFR or Bypass)

* MAIN_Transient.py:<br>
Main script for transient flow matching, single data only: WIP

* StageCalc.py:<br>
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


<br><br>
***

## ROSS
Customized ROSS Script for Turbo Compressor. Refer to ROSS for details.

<br><br>
***

## Thrust Balance
Thrust Calculation Tool<br>
* ThrustCalculator.py:<br>
  Main script<br>
  syntax:python ThrustCalculator.py<br>
### input
* input.csv:<br>
  Input measured data: pressure, temps etc.
* settings.ini:<br>
  Geometry Settings on impeller and turbine
### output
* geom.json:<br>
  output json file for geometry based on settings.ini

<br><br>
***


## Turbine Map Plot
Simple Turinbe Map Plotting Tool. Notebook format.<br>
* TurbineMapPlot_FixedGeometry.ipynb:<br>
  Main for fixed geometry and wastegated turbine

* TurbineMapPlot_VariableGeometry.ipynb:<br>
  Main for variable geometry turbine

<br><br>
***

## TurbineTemp_Calculator
Engine Combustion Temperature Calculation Tool
* T1t_Interpolate.py:<br>
  MAIN Calculation<br>
  syntax: T1t_Interpolate.py mapfile targetfile

* ShowMap.py:<br>
  T1t Map Plotting Tool Execute from MAIN

* curvefit.py:<br>
  scipy curvefit tool Execute from MAIN
### Input
* mapfile.csv:<br>
  Engine AFR and T1t Map
* targetfile.csv:<br>
  Target operating point
### Output
* Results

<br><br>
***
