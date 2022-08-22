# -*- coding: utf-8 -*-
"""
MIT License

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
from scipy import optimize
import json
import toml
import csv


if __name__ == '__main__':
    setting_file = './Inputs/setup.ini'
    setting = configparser.ConfigParser()
    setting.optionxform = str  # 大文字小文字を区別
    setting.read(setting_file, encoding='utf8')
    params_file = './Inputs/params.toml'
    with open(params_file) as inputf:
        params = toml.load(inputf)
    #cmap = np.loadtxt("./Inputs/cmap.csv", delimiter=',', skiprows=1) 
    #tmap = np.loadtxt("./Inputs/tmap.csv", delimiter=',', skiprows=1) 

#Class for each element calculation method
class Stage:
    def __init__(self, stgno, reload = False):  

        #Read Setting File (Global Variables)
        self.name = setting.get("Settings", "Analysis Name")
        self.P_Amb = setting.getfloat("Environment Variables", "Ambient Pressure [kPaA]")
        self.T_Amb = setting.getfloat("Environment Variables", "Ambient Temperature [K]")
        self.humidity = setting.getfloat("Environment Variables", "humidity [0-1]")
    

        #Read Setting File (Each Stage Variables) toml file
        self.C_Wt=params['Compressor']['Weight']
        self.C_I=params['Compressor']['Inertia']
        self.C_Cp=params['Compressor']['Cp']
        self.C_gam=params['Compressor']['gamma']
        self.C_flow=params['Compressor']['Flow']
        self.C_PRC=params['Compressor']['PRC']
        self.C_BLDR=params['Compressor']['BleedRatio']
        
        self.T_Wt=params['Turbine']['Weight']
        self.T_I=params['Turbine']['Inertia']
        self.T_Cp=params['Turbine']['Cp']
        self.T_gam=params['Turbine']['gamma']
        self.T_Pout=params['Turbine']['BackPressure']
        self.T_BLDR=params['Turbine']['BleedRatio']

        
        self.S_Wt=params['Shaft']['Weight']
        self.S_I=params['Shaft']['Inertia']

        self.CB_t=params['Combustor']['CombustionTemp_T1t']+273.15
        self.CB_dp=params['Combustor']['CombustordeltaP']
        self.CB_AFR=params['Combustor']['AirFuelRatio']
        self.CB_LHV=params['Combustor']['FuelLatentHeatLHV']

        self.S_Step=params['Calculation']['Step']
        self.S_PowRes=params['Calculation']['PowRes']

        
        """
        self.T_Comb = setting.getfloat("System Settings", "Combustion Temp (T1t) [K]") 
        self.Pdelta_Comb = setting.getfloat("System Settings", "Combustor delta P [%]")/100
        self.AFR = setting.getfloat("System Settings", "Air Fuel Ratio [kg/kg]") 
        self.HHV = setting.getfloat("System Settings", "Fuel Latent Heat HHV [MJ/kg]")*10**6
        self.Q1c_abs = setting.getfloat("Stage Settings", "Initial Inlet Air Flow [kg/s]")    
        self.Q1t_abs = self.Q1c_abs*(1+1/self.AFR)
        """


    def CompStage(self,stgno,Tin,Pin,N_C,W_C,PRC_Target):
        #Read Compressor Map
        df = pandas.read_csv("./Inputs/cmap_%s.csv" % (stgno), sep =",")
        f_PRC = interp2d(x = df['NC'], y = df['WC'], z = df['PRC'])
        f_EtaC = interp2d(x = df['NC'], y = df['WC'], z = df['EtaC'])
        PRC = f_PRC(N_C,W_C)
        EtaC = f_EtaC(N_C,W_C)
        Pout = Pin*PRC
        Tadb = Tin*(PRC**((self.C_gam-1)/self.C_gam)-1) #deltaT
        Tout = Tadb/EtaC + Tin
        Power = self.C_Cp*Tadb/EtaC*W_C
        print('NC',N_C,W_C,'____Tout',Tout)
        print('PwC',Power,EtaC,PRC)
        return Power, N_C, Pout,Tout
        #return self.PoutC,self.ToutC,self.EtaC,self.PowC,Cp


    def Combustor(self,W_C,Tin,Pin,AFR,C_BLDR,T_BLDR):
        W_F = W_C*(1-C_BLDR)/AFR
        W_T = (W_C+W_F)*(1-T_BLDR)
        deltaT = self.CB_LHV*10**6*W_F/(self.C_Cp*10**3*W_T)
        Tout = Tin + deltaT
        Pout = Pin*(1-self.CB_dp)
        return W_T,Pout,Tout


    def TurbStage(self,stgno,Tin,Pin,N_T,W_T):
        #Read Turbine Map
        df = pandas.read_csv("./Inputs/tmap_%s.csv" % (stgno), sep =",")
        f_PRT = interp2d(x = df['NT'], y = df['WT'], z = df['PRT'])
        f_EtaT = interp2d(x = df['NT'], y = df['WT'], z = df['EtaT'])
        PRT = f_PRT(N_T,W_T)
        EtaT = f_EtaT(N_T,W_T)
        Pout = Pin/PRT
        Tadb = Tin*(1-PRT**((self.T_gam-1)/self.T_gam))  #deltaT
        Tout = Tadb*EtaT + Tin
        Power = -self.T_Cp*Tadb*EtaT*W_T
        print('PwT',Power)
        return Power, N_T, Pout,Tout
    
    def PowBal(self,PwC,PwT,Nc,Nt,Step,PowRes):
        if np.abs(PwC - PwT) < PowRes:
            return 0,Nc,Nt
        if PwC > PwT:
            Nc = Nc - Step
            Nt = Nt + Step
            return 1,Nc,Nt
        elif PwT > PwC:
            Nc = Nc + Step
            Nt = Nt - Step
            return 1,Nc,Nt
        else:
            print("Error:Power does not converge, change Step or PowRes")
            

#Method for one time iteration matching calculation
def Matching(Stgno,Nc):
    Stg = Stage(Stgno)
    print(Stg.C_Wt)
    Nt = Nc
    Cont = 1
    while Cont == 1:
        PwC,Nc,P2c,T2c = Stg.CompStage(1,Stg.T_Amb,Stg.P_Amb,Nc,Stg.C_flow,Stg.C_PRC) 
        Wt,P1t,T1tth = Stg.Combustor(1,T2c,P2c,Stg.CB_AFR,Stg.C_BLDR,Stg.T_BLDR)
        PwT,Nt,P2t,T2t = Stg.TurbStage(1,Stg.CB_t,P1t,Nt,Wt) #Pin
        Cont,Nc,Nt = Stg.PowBal(PwC, PwT, Nc, Nt, Stg.S_Step,Stg.S_PowRes)
    print('Power kW',PwC,PwT,'RPM',Nc,'T2c[degC]',T2c-273.15,'T2t',T2t-273.15)
    print('Wt',Wt, 'P1t[kPa]',P1t,'T1t_th[degC]',T1tth-273.15)  

Matching(1,30000)


"""





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





#Output results
plt.close("all")
plt.ion()
pout = Stage(setting_file)
#pout.display()
pout.print()

#with open("StageResult.out","a") as output:
#    print("P2in :\t\t%.2f " % (P2in),file=output)




#Adiabatic Specific Heat Ratio Table



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
