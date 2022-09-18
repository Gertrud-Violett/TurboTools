"""
-*- coding: utf-8 -*-
Copyright (c) makkiblog.com
MIT License 
-*- SingleStage Power Matching Tool for Turbocompressor -*-

-*- SYNTAX USAGE -*-
>python SingleStage_Matching.py
vvvCODEvvv
"""

import sys
import os
import math
import configparser
import numpy as np
import matplotlib.pyplot as plt
import importlib
import pandas
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
import toml
import csv
from StageCalc import Stage

#import setting file
if __name__ == '__main__':
    setting_file = './Inputs/setup.ini'
    setting = configparser.ConfigParser()
    setting.optionxform = str  #Case sensitve option
    setting.read(setting_file, encoding='utf8')
    initrpm = setting.getfloat("Settings", "Initial rpm")
    humidity = setting.getfloat("Environment Variables", "humidity [0-1]") #not used yet
    matchingmode = setting.get("Settings", "Matching Mode (Bypass or AFR)") # matching mode Bypass or AFR correction
    analysisname = setting.get("Settings", "Analysis Name")
    params_file = setting.get("Settings", "Input filename")


#Method for one time iteration matching calculation for flow bypass mode
def MatchingBypassMode(Stgno,Nc):
    Stg = Stage(Stgno,params_file)
    print("Flow",Stg.C_Wt)
    Nt = Nc
    StgEff = 1 - Stg.S_MechLoss
    BackPres = Stg.T_Pout
    Pin_Tot = Stg.P_Amb
    Cont = 1
    while Cont == 1: #match compressor
        PwC,Wcstar,Nc,P2c,T2c,PRC,EtaC = Stg.CompStage(1,Stg.T_Amb,Pin_Tot,Nc,Stg.C_flow,Stg.C_PRC) 
        Cont,Nc,Nt = Stg.PowBal(PRC,Nc, Nt, Stg.S_Step)
    W_F,Wt,P1t,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,Stg.CB_AFR,Stg.C_BLDR,Stg.T_BLDR)
    T1t_calc = T1tth*Stg.CB_Tempeff
    PwT,W_reqd,W_Tstar,PwAct,Nt,P2t,T2t,PRT,BPR,EtaT = Stg.TurbStage(1,StgEff,PwC,T1t_calc,P1t,Nt,Wt,BackPres)
    Wt_act = Wt
    Enth_act = PwT
    AFR_act = Stg.CB_AFR
    while Enth_act > PwAct + Stg.S_PowRes:
        Enth_act,W_reqd,W_Tstar,PwAct,Nt,P2t,T2t,PRT,BPR,EtaT = Stg.TurbStage(1,StgEff,PwC,T1t_calc,P1t,Nt,Wt_act,BackPres)
        Wt_act -= (Wt-W_reqd)*BPR*Stg.Wtres
    print('Matched power[kW] %.1f %.1f:' % (Enth_act, PwAct))
    BPR = (Wt - W_reqd)/Wt
    
    dict_c = {"Power[kW]":PwC,"Corrected Flow rate Wc*":Wcstar, "Speed RPM":Nc,"Outlet Pres[kPaA]":P2c,"Outlet Temp[degC]":T2c-273.15,"Pres Ratio PRC":PRC,"Comp Efficiency":EtaC,"Air Flow[kg/s]":Stg.C_flow,"Turbine Pressure Ratio P2c/P1t":P2c/P1t}
    dict_t = {"Enthalpy[kW]":PwT,"Turbine Power[kW]":PwAct, "Speed RPM":Nc,"Inlet Pres.[kPa]":P1t, "Inlet Temp[C]":T1t_calc-273.15,"Outlet Pres[kPa]Tot":P2t, "Outlet Temp[degC]":T2t-273.15, "Pres Ratio PRT":PRT, "Turb Efficiency":EtaT,"Total Flow[kg/s]":Wt,"Turbine Flow[kg/s]":Wt_act,"Bypass Ratio":BPR}
    dict_s = {"Fuel Flow [kg/s]":W_F,"Turbine Required Flow[kg/s]":W_reqd,"Corrected Flow rate[kg/s]":W_Tstar,"Theoretical T1t[degC]":T1tth-273.15,"Mechanical Loss [kW]":(1-StgEff)*PwAct,"AFR": AFR_act}
    return dict_c,dict_t,dict_s


#Method for one time iteration matching calculation for T1t correction mode
def MatchingAFRMode(Stgno,Nc):
    Stg = Stage(Stgno,params_file)
    print("Flow",Stg.C_Wt)
    Nt = Nc
    StgEff = 1 - Stg.S_MechLoss
    BackPres = Stg.T_Pout
    Pin_Tot = Stg.P_Amb
    Cont = 1
    while Cont == 1: #match compressor
        PwC,Wcstar,Nc,P2c,T2c,PRC,EtaC = Stg.CompStage(1,Stg.T_Amb,Pin_Tot,Nc,Stg.C_flow,Stg.C_PRC) 
        Cont,Nc,Nt = Stg.PowBal(PRC,Nc, Nt, Stg.S_Step)
    W_F,Wt,P1t,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,Stg.CB_AFR,Stg.C_BLDR,Stg.T_BLDR)
    T1t_calc = T1tth*Stg.CB_Tempeff
    PwT,W_reqd,W_Tstar,PwAct,Nt,P2t,T2t,PRT,BPR,EtaT = Stg.TurbStage(1,StgEff,PwC,T1t_calc,P1t,Nt,Wt,BackPres)
    Wt_act = Wt
    Enth_act = PwT
    AFR_act = Stg.CB_AFR
    print(AFR_act)
    while Enth_act > PwAct + Stg.S_PowRes:
        Enth_act,W_reqd,W_Tstar,PwAct,Nt,P2t,T2t,PRT,BPR,EtaT = Stg.TurbStage(1,StgEff,PwC,T1t_calc,P1t,Nt,Wt_act,BackPres)
        AFR_act += (Wt-W_reqd)/Wt*BPR*Stg.Wtres*100
        W_F,Wt,P1t,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,AFR_act,Stg.C_BLDR,Stg.T_BLDR)
        T1t_calc = T1tth*Stg.CB_Tempeff

    print('Matched power[kW] %.1f %.1f:' % (Enth_act, PwAct))
    BPR = (Wt - W_reqd)/Wt
    
    dict_c = {"Power[kW]":PwC,"Corrected Flow rate Wc*":Wcstar, "Speed RPM":Nc,"Outlet Pres[kPaA]":P2c,"Outlet Temp[degC]":T2c-273.15,"Pres Ratio PRC":PRC,"Comp Efficiency":EtaC,"Air Flow[kg/s]":Stg.C_flow,"Turbine Pressure Ratio P2c/P1t":P2c/P1t}
    dict_t = {"Enthalpy[kW]":PwT,"Turbine Power[kW]":PwAct, "Speed RPM":Nc,"Inlet Pres.[kPa]":P1t, "Inlet Temp[C]":T1t_calc-273.15,"Outlet Pres[kPa]Tot":P2t, "Outlet Temp[degC]":T2t-273.15, "Pres Ratio PRT":PRT, "Turb Efficiency":EtaT,"Total Flow[kg/s]":Wt,"Turbine Flow[kg/s]":Wt_act,"Bypass Ratio":BPR}
    dict_s = {"Fuel Flow [kg/s]":W_F,"Turbine Required Flow[kg/s]":W_reqd,"Corrected Flow rate[kg/s]":W_Tstar,"Theoretical T1t[degC]":T1tth-273.15,"Mechanical Loss [kW]":(1-StgEff)*PwAct,"AFR": AFR_act}
    #list_other = [W_reqd,WT_Map,PwAct,Nt,BPR,Cont,Nc,Nt,W_F,P1t,]
    return dict_c,dict_t,dict_s



#Run Matching
if matchingmode == "Bypass":
    dict_c,dict_t,dict_s = MatchingBypassMode(1,initrpm)
elif matchingmode == "AFR":
    dict_c,dict_t,dict_s = MatchingAFRMode(1,initrpm)
else:
    print("Matching mode must be selected")
    exit()

#Export toml results
for k,v in dict_c.items():
    dict_c[k] = float("{:.3f}".format(float(v)))
for k,v in dict_t.items():
    dict_t[k] = float("{:.3f}".format(float(v)))
for k,v in dict_s.items():
    dict_s[k] = float("{:.3f}".format(float(v)))

dict = {"Compressor":dict_c,
"Turbine":dict_t,
"Shaft&Combustor":dict_s}

toml_string = toml.dumps(dict)
output_file = "./Outputs/Result_"+matchingmode+"_"+analysisname+".toml"
with open(output_file,"w") as output:
    toml.dump(dict,output)