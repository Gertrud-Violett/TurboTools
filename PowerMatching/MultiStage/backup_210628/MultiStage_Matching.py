# -*- coding: utf-8 -*-
"""
Created on Thu Jun-17 2021 by GV

"""
import sys
import os
import math
import configparser
import numpy as np
import matplotlib.pyplot as plt
import importlib
import pandas
from scipy import interpolate
from scipy import optimize
#from scipy.optimize import newton
#import xlrd
import json
#import csv
from GasTable import GasProp


if __name__ == '__main__':
    setting_file = 'MultiStage_Matching_input.ini'
    #stage_params = 'stage_params.txt'
    stage_json = open('stage_params.json','r')
    stage_params = json.load(stage_json)
    
class Stage:
    def __init__(self, setting_file, reload = False):
        self.setting_file = setting_file
        self.setting = configparser.ConfigParser()
        self.setting.optionxform = str  # 大文字小文字を区別するおまじない
        self.setting.read(setting_file, encoding='utf8')
        setting = self.setting

        self.stage_params = stage_params
        self.stagepar = configparser.ConfigParser()
        self.stagepar.optionxform = str  # 大文字小文字を区別するおまじない
        self.stagepar.read(stage_params, encoding='utf8')
        stagepar = self.stagepar

        #Read Setting File (Global Variables)
        self.name = setting.get("Settings", "Analysis Name")
        self.P_Amb = setting.getfloat("Environment Variables", "Ambient Pressure [kPaA]")
        self.T_Amb = setting.getfloat("Environment Variables", "Ambient Temperature [K]")
        self.T_Comb = setting.getfloat("System Settings", "Combustion Temp (T1t) [K]") 
        self.Pdelta_Comb = setting.getfloat("System Settings", "Combustor delta P [kPa]")*10**3
        self.AFR = setting.getfloat("System Settings", "Air Fuel Ratio [kg/kg]") 
        self.HHV = setting.getfloat("System Settings", "Fuel Latent Heat HHV [MJ/kg]")*10**6
        
        self.Q1c_abs = setting.getfloat("Stage Settings", "Inlet Air Flow [kg/s]")
        self.C_stg_qty = setting.getfloat("Stage Settings", "No of Compressor Stages")
        self.T_stg_qty = setting.getfloat("Stage Settings", "No of Turbine Stages")

        self.Spool_Qty = setting.getfloat("Stage Settings", "Shaft Spool Qty")

        #Spool Inertia
        #Future Implementation
        
        #Read Setting File (Each Stage Variables): PR,Eta,deltaT,deltaP
        #Define Parameters: CsetX: Compressor Stage X, TsetX: Turbine Stage X
        #PRC_X: Compressor Stage Pressure Ratio for Stage X
        #EtaC_X: Compressor Stage Efficiency for Stage X 
        #Inertia_Y: Spool Rotating Group Total Inertia
        #SplNo: Defines Spool No.for Stage X:Support up to 4 spools
        #deltaT: Defines Intercooler Temperature drop between stages (Compressor Only)
        #deltaP: Defines Intercooler Pressure drop between stages (Compressor Only)
        #Inertia: Spool Total Inertia including all rotors and shaft
        self.Cset1 = stage_params["Compressor1"]
        self.Cset2 = stage_params["Compressor2"]
        self.Cset3 = stage_params["Compressor3"]
        self.Tset1 = stage_params["Turbine1"]
        self.Tset2 = stage_params["Turbine2"]

        #Define initial stage efficiency
        self.EtaC1 = 0.9
        self.EtaT1 = 0.9


    def CompStage(self,stage,Tin,Pin,PRC,Cp,Cv,dT,dP,Q):
        #Read Compressor Map
        df = pandas.read_csv("./WheelMap/%s_CompMap_%s.csv" % (self.name, stage), sep =",")
        PaiC = df["PRC"]
        WC = df["WC"]
        N_c = df["NC"]
        index = ((PaiC - PRC)**2+(WC - Q)**2).idxmin()
        self.EtaC = df.iloc[index,3]
        #self.rho = 1.293*(Pres/101.325)*(Temp/298.15)**0.5
        self.gamma = Cp/Cv
        self.Pout = Pin*PRC - dP
        self.Toutadb = Tin*PRC**((self.gamma-1)/self.gamma)
        self.Tout = self.Toutadb/self.EtaC - dT

        return self.Tout,self.EtaC,Cp,Cv


    def TurbStage(self):
        self.T2c = 1+1 



#display function global Parameters
    def display(self): #CLI上に結果を表示
        print("")
        print("Analysis Name :\t\t%s " % (self.name))
        print("Ambient Temp :\t\t%.2f " % (self.P_Amb))
        print("")
        print (self.Cset1,self.EtaC1)

    def print(self): #.outファイルに結果を出力
        with open("StageResult.out","w") as output:
            print("Analysis Name :\t\t%s " % (self.name),file=output)
            print("Ambient Temp :\t\t%.2f " % (self.P_Amb),file=output)



#COMPRESSOR SIDE Each Stage Calculation 
#1st stage only ambient settings, initial calculation
C_1 = Stage(setting_file)
T1in = C_1.T_Amb 
P1in = C_1.P_Amb 

#obtain from GasTable.py
Cp1 = 1.006
Cv1 = 0.7171
gamma1 = Cp1/Cv1
rho = 1.293*(P1in/101.325)*(T1in/298.15) 
C_1.CompStage(1,T1in,P1in,float(C_1.Cset1["PRC"]),Cp1,Cv1,float(C_1.Cset1["deltaT"]),float(C_1.Cset1["deltaP"]),0.223)
#C_1.CompStage(C_1.T_Amb,C_1.P_Amb,C_1.Cset1["PRC"],C_1.Cset1["EtaC"],1.293,1.4,1,C_1.Cset1["deltaT"],C_1.Cset1["deltaP"])
P2in = C_1.Pout
T2in = C_1.Tout
print(C_1.EtaC)
#print(C_1.EtaC,rho, gamma1, P2in, T2in,C_1.Cset1["PRC"])

        

#COMBUSTOR Calculation

#TURBINE SIDE Each Stage Calculation



#Output results
plt.close("all")
plt.ion()
pout = Stage(setting_file)
#pout.display()
#pout.print()

with open("StageResult.out","a") as output:
    print("P2in :\t\t%.2f " % (P2in),file=output)




#Adiabatic Specific Heat Ratio Table


"""
#Obtain from Refprop

        coolant = Coolant()
        coolant.load_axis('./TableData/Methane_table/axes.csv')
        coolant.load_table('./TableData/Methane_table/density.csv', 'rho')

        coolant.load_table('./TableData/Methane_table/temperature.csv', 'temperature', rowAxis='pressure',colAxis='enthalpy')
       
        rho_f = coolant.rho(T_fup, P_fup*1e6)
        gamma_f = coolant.gamma(T_fup, P_fup*1e6)

        #flow velocity
        mdot_o = mdot_os * (1 - mr_Tert_o)
        self.LOx_vo1 = mdot_o*self.LOx_PA_ratio/rho_o/self.Area_LOx1

        #Pressure loss delta p
        self.deltap_o =  mdot_o**2/(2*rho_o*self.Area_LOx**2*LOx_Cd**2)/10**6        
        
   
"""

"""

def Turbopomp_power_calculation(g,dot_m_F,dot_m_O,P_GG_center,dot_m_GG,T_GG,Cp_GG,P_OT,gamma_GG,dH_F,dH_O,eta_pf,eta_po,eta_t):
    
    dH_t=Cp_GG*T_GG*(1-(P_OT/P_GG_center)**((gamma_GG-1)/gamma_GG))
    
    L_PO=g*dot_m_O*dH_O/eta_po
    L_PF=g*dot_m_F*dH_F/eta_pf
    L_t=dot_m_GG*eta_t*dH_t

    return L_PF,L_PO,L_t

"""
