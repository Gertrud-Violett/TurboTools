"""
-*- coding: utf-8 -*-
Copyright (c) makkiblog.com
MIT License 
-*- SingleStage Power Matching Tool for Turbocompressor -*-
stage calc subroutine
v0.2 Improved convergence

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
import warnings
warnings.filterwarnings("ignore")

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

#Class for each stage element calculation methods
class Stage:
    def __init__(self, stgno, params_file, reload = False):  
        params_set = './Inputs/' + params_file + '_' + str(stgno) + '.toml'
        with open(params_set) as inputf:
            params = toml.load(inputf)
        #Read Setting File (Each Stage Variables) toml file
        self.P_Amb = params['Compressor']['AmbientPressure']
        self.T_Amb = params['Compressor']['AmbientTemperature']
        self.C_Wt=params['Compressor']['Weight']
        self.C_I=params['Compressor']['Inertia']
        self.C_Cp=(0.2208+51.48*10**-6*(self.T_Amb))*360/86
        self.C_gam=params['Compressor']['gamma']
        self.C_flow=params['Compressor']['Flow']
        self.C_PRC=params['Compressor']['PRC']
        self.C_BLDR=params['Compressor']['BleedRatio']

        self.CB_t=params['Combustor']['CombustionTemp_T1t']
        self.CB_dp=params['Combustor']['CombustordeltaP']
        self.CB_AFR=params['Combustor']['AirFuelRatio']
        self.CB_LHV=params['Combustor']['FuelLatentHeatLHV']
        self.CB_Tempeff=params['Combustor']['TempEfficiency']

        self.T_Wt=params['Turbine']['Weight']
        self.T_I=params['Turbine']['Inertia']
        self.T_Cp=(0.2208+51.48*10**-6*(self.CB_t+273.15))*360/86
        self.T_gam=1.42-84.6*10**-6*(self.CB_t+273.15) 
        self.T_Pout=params['Turbine']['BackPressure']
        self.T_BLDR=params['Turbine']['BleedRatio']
        
        self.S_Wt=params['Shaft']['Weight']
        self.S_I=params['Shaft']['Inertia']
        self.S_MechLoss=params['Shaft']['MechLoss']


        self.S_Step=params['Calculation']['Step']
        self.S_PowRes=params['Calculation']['PowRes']
        self.S_PRCRes=params['Calculation']['PRCRes']
        self.Wtres=params['Calculation']['WtRes']
        self.T1tres=self.Wtres

    #def showmap(self,stgno):

    def CompStage(self,stgno,Tin,Pin,N_C,W_C,PRC_Target):
        #Read Compressor Map
        df = pandas.read_csv("./Inputs/cmap_%s.csv" % (stgno), sep =",")
        W_Cstar = W_C*(Tin/288.15)**0.5*101.325/Pin  #Corrected flow rate
        N_Cstar = N_C/(Tin/288.15)**0.5 #Corrected Speed
        if N_Cstar > max(df['NC*']):
            print ("ERROR:Initial rpm too large")
            exit()
        elif N_Cstar < min(df['NC*']):
            print ("ERROR:Initial rpm too small")
            exit()
        else:
            f_PRC = interp2d(x = df['NC*'], y = df['WC*'], z = df['PRC'])
            f_EtaC = interp2d(x = df['NC*'], y = df['WC*'], z = df['EtaC'])
            PRC = f_PRC(N_Cstar,W_Cstar)
            EtaC = f_EtaC(N_Cstar,W_Cstar)
            Pout = Pin*PRC
            Tadb = Tin*(PRC**((self.C_gam-1)/self.C_gam)-1) #deltaT
            Tout = Tadb/EtaC + Tin
            Power = self.C_Cp*Tadb/EtaC*W_C
            print('COMPRESSOR','PwC[kW]:%.1f|' % Power,'NC[rpm]:%.0f|' % N_C,'NC*:%.0f|' % N_Cstar,'WC[kg/s]:%.2f|' % W_C,'WC*:%.2f|' % W_Cstar,'EtaC:%.3f|' % EtaC,'PRC:%.3f|' % PRC)
            return Power,W_Cstar,N_C,Pout,Tout,PRC,EtaC

    def Combustor(self,W_C,Tin,Pin,AFR,C_BLDR,T_BLDR):
        W_Comb = W_C*(1-C_BLDR)
        W_F = W_Comb/AFR
        W_T = (W_Comb+W_F)*(1-T_BLDR)
        deltaT = self.CB_LHV*10**6*W_F/(self.C_Cp*10**3*W_T)
        Tout = Tin + deltaT
        Pout = Pin*(1-self.CB_dp)
        print('COMBUSTOR','P2c[kPa]:%.1f|'% Pin,'P1t[kPa]:%.1f|' % Pout)
        return W_F,W_T,Pout,Tout

    def TurbStage(self,stgno,EtaM,PwReq,Tin,Pin,N_T,W_T,Pout):
        #Read Turbine Map
        W_Tstar = W_T*(Tin/288.15)**0.5*101.325/Pin #Corrected flow rate
        N_Tstar = N_T/(Tin/288.15)**0.5 #Corrected Speed
        df = pandas.read_csv("./Inputs/tmap_%s.csv" % (stgno), sep =",")
        f_PRT = interp2d(x = df['NT*'], y = df['WT*'], z = df['PRT'])
        f_EtaT = interp2d(x = df['NT*'], y = df['WT*'], z = df['EtaT'])
                        
        PRT = f_PRT(N_Tstar,W_Tstar)
        EtaT = f_EtaT(N_Tstar,W_Tstar)
        Tadb = Tin*(1-PRT**((self.T_gam-1)/self.T_gam))  #deltaT
        Tout = Tadb*EtaT + Tin
        Enth = -self.T_Cp*Tadb*EtaT*W_T
        PwAct = PwReq/EtaM
        W_reqd = -PwAct/(self.T_Cp*Tadb*EtaT)
        W_reqdstar = W_reqd*(Tin/288.15)**0.5*101.325/Pin
        #W_Tinit = W_reqdstar
        BypassRatio = (W_T - W_reqd)/W_T
        Poutact = Pin/PRT
        if BypassRatio < 0:
            print("ERROR: Bypass Ratio must be larger than 0. RESULTS WILL HAVE ERROR")
            print('TURBINE','Tot. Enthal[kW]:%.1f|' % Enth,'PwT Reqd[kW]:%.1f|' % PwAct,'NT[rpm]:%.0f|' % N_T,'NT*:%.0f|' % N_Tstar,'W_Tot[kg/s]:%.2f|' % W_T,'W_Tot*:%.2f|' % W_Tstar,'W_Reqd Turb[kg/s]:%.2f|' % W_reqd,'W_Req Turb*:%.2f|' % W_reqdstar,'EtaT:%.3f|' % EtaT,'PRT:%.3f|' % PRT,'BPR:%.2f|' % BypassRatio,'P2t Actual[kPaA]:%.1f|' % Poutact)
            return Enth,W_reqd,W_Tstar,PwAct,N_T,Poutact,Tout,PRT,BypassRatio,EtaT
        else:
            print('TURBINE','Tot. Enthal[kW]:%.1f|' % Enth,'PwT Reqd[kW]:%.1f|' % PwAct,'NT[rpm]:%.0f|' % N_T,'NT*:%.0f|' % N_Tstar,'W_Tot[kg/s]:%.2f|' % W_T,'W_Tot*:%.2f|' % W_Tstar,'W_Reqd Turb[kg/s]:%.2f|' % W_reqd,'W_Req Turb*:%.2f|' % W_reqdstar,'EtaT:%.3f|' % EtaT,'PRT:%.3f|' % PRT,'BPR:%.2f|' % BypassRatio,'P2t Actual[kPaA]:%.1f|' % Poutact)
            return Enth,W_reqd,W_Tstar,PwAct,N_T,Poutact,Tout,PRT,BypassRatio,EtaT

    def PowBal(self,PRC,Nc,Nt,Step): #Compressor side power matching
        if PRC > self.C_PRC+self.S_PRCRes:
            Nc =  Nc + Nc*1/2*((self.C_PRC/PRC)**(1/3) - 1) 
            return 1,Nc,Nt
        elif PRC < self.C_PRC-self.S_PRCRes:
            Nc = Nc + Nc*1/2*((self.C_PRC/PRC)**(1/3) - 1) 
            return 1,Nc,Nt
        else:
            Nt = Nc
            return 0,Nc,Nt

    
    def PowBalRot(self,PRC,Nc,Nt,Step): #Option for unbalanced power and N rot turbine sliding mode
        if PRC > self.C_PRC+self.S_PRCRes:
            Nc =  Nc + Nc*1/2*((self.C_PRC/PRC)**(1/3) - 1) 
            return 1,Nc,Nt
        elif PRC < self.C_PRC-self.S_PRCRes:
            Nc = Nc + Nc*1/2*((self.C_PRC/PRC)**(1/3) - 1) 
            return 1,Nc,Nt
        else:
            Nt = Nc
            return 0,Nc,Nt

            """
            if np.abs(Nc - Nt) < Step:
                return 0,Nc,Nt
            if Nc > Nt - Step:
                Nt = Nt+1/2*(Nc-Nt)
                return 1,Nc,Nt
            elif Nt > Nc - Step:
                Nt = Nt-1/2*(Nc-Nt)
                return 1,Nc,Nt
            else:
                print("Error:Power does not converge, change Step or PowRes")
               """ 
