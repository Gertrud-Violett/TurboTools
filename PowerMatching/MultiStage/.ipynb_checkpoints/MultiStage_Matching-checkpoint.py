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
from GasTable import AirProp
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
        self.humidity = setting.getfloat("Environment Variables", "humidity [0-1]")
        self.T_Comb = setting.getfloat("System Settings", "Combustion Temp (T1t) [K]") 
        self.Pdelta_Comb = setting.getfloat("System Settings", "Combustor delta P [%]")/100
        self.AFR = setting.getfloat("System Settings", "Air Fuel Ratio [kg/kg]") 
        self.HHV = setting.getfloat("System Settings", "Fuel Latent Heat HHV [MJ/kg]")*10**6
        
        self.Q1c_abs = setting.getfloat("Stage Settings", "Inlet Air Flow [kg/s]")
        self.C_stg_qty = setting.getfloat("Stage Settings", "No of Compressor Stages")
        self.T_stg_qty = setting.getfloat("Stage Settings", "No of Turbine Stages")
        
        self.Q1t_abs = self.Q1c_abs*1.05

        self.Spool_Qty = setting.getfloat("Stage Settings", "Shaft Spool Qty")

        #Spool Definition. N1>N2>N3>Combustor>N3>N2>N1
        #Spool Inertia
        #Future Implementation
        
        #Read Setting File (Each Stage Variables): PR,Eta,deltaT,deltaP
        #Define Parameters: CsetX: Compressor Stage X, TsetX: Turbine Stage X
        #PRC_X: Compressor Stage Pressure Ratio for Stage X
        #EtaC_X: Compressor Stage Efficiency for Stage X 
        #Inertia_Y: Spool Rotating Group Total Inertia
        #SplNo: Defines Spool No.for Stage X:Support up to 4 spools
        #deltaT: Defines Intercooler Temperature drop between stages (Compressor Only) or combustor temp rise in delta K
        #deltaP: Defines Intercooler Pressure drop between stages (Compressor Only) or combustor pressure loss in delta kPa
        #Inertia: Spool Total Inertia including all rotors and shaft
        #W_Pw: kW Power loss in each shaft
        self.Combustor = stage_params["Combustor"]
        self.Cset1 = stage_params["Compressor1"]
        self.Cset2 = stage_params["Compressor2"]
        self.Cset3 = stage_params["Compressor3"]
        self.Tset1 = stage_params["Turbine1"]
        self.Tset2 = stage_params["Turbine2"]
        self.Sp1 = stage_params["Spool1"]


    def CompStage(self,stage,Tin,Pin,Nc,Cp,gamma,rho,dT,dP,Q):
        #Read Compressor Map
        df = pandas.read_csv("./WheelMap/%s_CompMap_%s.csv" % (self.name, stage), sep =",")
        #PaiC = df["PRC"]
        self.gamma = gamma
        WC = df["WC"]
        N_c = df["NC"]
        #index = ((WC - Q)**2+(N_c - Nc)**2).idxmin()
        #self.EtaC = df.iloc[index,3]/100
        #self.PRC = df.iloc[index,2]
        self.EtaC = 0.83
        self.PRC = 15
        self.PoutC = Pin*self.PRC - dP
        self.TadbC = Tin*(self.PRC**((self.gamma-1)/self.gamma)-1)
        self.ToutC = self.TadbC/self.EtaC - dT + Tin
        self.PowC = Cp*self.TadbC/self.EtaC
        
        return self.PoutC,self.ToutC,self.EtaC,self.PowC,Cp


    def TurbStage(self,stage,Tin,Pin,Nt,Cp,gammaT,rho,Q):
        #Read Turbine Map
        df = pandas.read_csv("./WheelMap/%s_TurbMap_%s.csv" % (self.name, stage), sep =",")
        self.gammaT = gammaT
        WT = df["WT"]
        N_t = df["NT"]
        #index = ((WT - Q)**2+(N_t - Nt)**2).idxmin()
        #self.EtaT = df.iloc[index,3]/100
        #self.PRT = df.iloc[index,2]
        self.EtaT = 0.85
        self.PRT = 10
        self.PoutT = Pin/self.PRT
        self.TadbT = Tin*(1-self.PRT**((self.gammaT-1)/self.gammaT))
        self.ToutT = self.TadbT*self.EtaT + Tin
        self.PowT = -Cp*self.TadbT/self.EtaT
        
        return self.PoutT,self.ToutT,self.EtaT,self.PowT,Cp




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



#COMPRESSOR SIDE Each Stage Calculation =====================================================================
#1st stage only ambient settings, initial calculation
C_1 = Stage(setting_file)
T1in = C_1.T_Amb 
P1in = C_1.P_Amb 

#input: Temp,Pres,N2_MF,O2_MF,CO2_MF,humidity,Ar_MF
#output: rho, Cv, Cp, gamma,P,Pi, H2O_MF, Ps
C_1_Prop = AirProp(T1in,P1in,0.78,0.21,0.0007,C_1.humidity,0.0093)
rho1 = C_1_Prop[0]
Cp1 = C_1_Prop[2]
gamma1 = C_1_Prop[3]
Nc1_init = float(C_1.Sp1["Ncinit"])
dT1C = float(C_1.Cset1["deltaT"])
dP1C = float(C_1.Cset1["deltaP"])
print(Cp1)

#input (stage,Tin,Pin,Nt,Cp,gamma,rho,dT,dP,Q)
#output PoutT,ToutT,EtaT,PowT
C_1_out = C_1.CompStage(1,T1in,P1in,Nc1_init,Cp1,gamma1,rho1,dT1C,dP1C,C_1.Q1c_abs)
P2in = C_1.PoutC
T2in = C_1.ToutC
#print(C_1.EtaC)
print(C_1_out)


        

#COMBUSTOR Calculation =====================================================================
CB = Stage(setting_file)
PCBin = P2in
TCBin = T2in

TCBout = CB.T_Comb
PCBout = PCBin*(1 - CB.Pdelta_Comb) 
MF = (float(CB.Combustor["N2MF"]),float(CB.Combustor["O2MF"]),float(CB.Combustor["CO2MF"]),float(CB.Combustor["COMF"]),float(CB.Combustor["H2OMF"]),float(CB.Combustor["ARMF"]))
print(TCBout, PCBout)


#TURBINE SIDE Each Stage Calculation =====================================================================
T_1 = Stage(setting_file)
T_1_Prop = GasProp(TCBout,PCBout,MF[0],MF[1],MF[2],MF[3],MF[4],MF[5],) 
rhoT1 = T_1_Prop[0]
CpT1 = T_1_Prop[2]
gammaT1 = T_1_Prop[3]
Nc1T_init = float(T_1.Sp1["Ncinit"])

#input (stage,Tin,Pin,Nt,Cp,gamma,rho,dT,dP,Q)
#output PoutT,ToutT,EtaT,PowT
T_1_out = T_1.TurbStage(1,TCBout,PCBout,Nc1T_init,CpT1,gammaT1,rhoT1,T_1.Q1t_abs)
print(T_1_out)




#Output results
plt.close("all")
plt.ion()
pout = Stage(setting_file)
#pout.display()
pout.print()

#with open("StageResult.out","a") as output:
#    print("P2in :\t\t%.2f " % (P2in),file=output)




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
