# -*- coding: utf-8 -*-
"""
MIT License
SingleStage Power Matching Tool for Turbocompressor
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


if __name__ == '__main__':
    setting_file = './Inputs/setup.ini'
    setting = configparser.ConfigParser()
    setting.optionxform = str  #Case sensitve option
    setting.read(setting_file, encoding='utf8')
    params_file = './Inputs/params.toml'
    with open(params_file) as inputf:
        params = toml.load(inputf)

#Class for each stage element calculation methods
class Stage:
    def __init__(self, stgno, reload = False):  

        #Read Setting File (Global Variables)
        self.name = setting.get("Settings", "Analysis Name")
        self.P_Amb = setting.getfloat("Environment Variables", "Ambient Pressure [kPaA]")
        self.T_Amb = setting.getfloat("Environment Variables", "Ambient Temperature [K]")
        self.humidity = setting.getfloat("Environment Variables", "humidity [0-1]") #not used yet

        #Read Setting File (Each Stage Variables) toml file
        self.C_Wt=params['Compressor']['Weight']
        self.C_I=params['Compressor']['Inertia']
        self.C_Cp=(0.2208+51.48*10**-6*(self.T_Amb))*360/86
        self.C_gam=params['Compressor']['gamma']
        self.C_flow=params['Compressor']['Flow']
        self.C_PRC=params['Compressor']['PRC']
        self.C_BLDR=params['Compressor']['BleedRatio']

        self.CB_t=params['Combustor']['CombustionTemp_T1t']+273.15
        self.CB_dp=params['Combustor']['CombustordeltaP']
        self.CB_AFR=params['Combustor']['AirFuelRatio']
        self.CB_LHV=params['Combustor']['FuelLatentHeatLHV']

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


    def CompStage(self,stgno,Tin,Pin,N_C,W_C,PRC_Target):
        #Read Compressor Map
        df = pandas.read_csv("./Inputs/cmap_%s.csv" % (stgno), sep =",")
        if N_C > max(df['NC']):
            print ("ERROR:Initial rpm too large")
            return 0
        elif N_C < min(df['NC']):
            print ("ERROR:Initial rpm too small")
            return 0
        else:
            f_PRC = interp2d(x = df['NC'], y = df['WC'], z = df['PRC'])
            f_EtaC = interp2d(x = df['NC'], y = df['WC'], z = df['EtaC'])
            PRC = f_PRC(N_C,W_C)
            EtaC = f_EtaC(N_C,W_C)
            Pout = Pin*PRC
            Tadb = Tin*(PRC**((self.C_gam-1)/self.C_gam)-1) #deltaT
            Tout = Tadb/EtaC + Tin
            Power = self.C_Cp*Tadb/EtaC*W_C
            print('NC [RPM]',N_C,'WC[kg/s]',W_C)
            print('PwC',Power,'EtaC',EtaC,'PRC',PRC)
            return Power,N_C,Pout,Tout,PRC,EtaC

    def Combustor(self,W_C,Tin,Pin,AFR,C_BLDR,T_BLDR):
        W_F = W_C*(1-C_BLDR)/AFR
        W_T = (W_C+W_F)*(1-T_BLDR)
        deltaT = self.CB_LHV*10**6*W_F/(self.C_Cp*10**3*W_T)
        Tout = Tin + deltaT
        Pout = Pin*(1-self.CB_dp)
        print('P2c',Pin,'P1t',Pout)
        return W_T,Pout,Tout

    def TurbStage(self,stgno,EtaM,PwReq,Tin,Pin,N_T,W_T,Pout):
        #Read Turbine Map
        df = pandas.read_csv("./Inputs/tmap_%s.csv" % (stgno), sep =",")
        f_WT = interp2d(x = df['NT'], y = df['PRT'], z = df['WT'])
        f_EtaT = interp2d(x = df['NT'], y = df['PRT'], z = df['EtaT'])
        PRT = Pin/Pout
        #WT = f_WT(N_T,PRT)
        EtaT = f_EtaT(N_T,PRT)
        Tadb = Tin*(1-PRT**((self.T_gam-1)/self.T_gam))  #deltaT
        Tout = Tadb*EtaT + Tin
        Power = -self.T_Cp*Tadb*EtaT*W_T*EtaM
        W_reqd = -PwReq/(self.T_Cp*Tadb*EtaT*EtaM)
        #EtaT_w = f_EtaT(N_T,PRT)

        W_bypass = W_T - W_reqd
        BypassRatio = W_bypass/W_T
        if BypassRatio < 0:
            print("ERROR: Bypass Ratio must be larger than 0")
            print('PwT',Power, 'EtaT',EtaT,'PRT',PRT,'BPR',BypassRatio,'W_T',W_T,'NT',N_T,'W_Reqd',W_reqd)
            return Power,W_reqd,N_T,Pout,Tout,PRT,BypassRatio,EtaT
    
        else:
            print('PwT',Power, 'EtaT',EtaT,'PRT',PRT,'BPR',BypassRatio,'W_T',W_T,'NT',N_T,'W_Reqd',W_reqd)
            return Power,W_reqd,N_T, Pout,Tout,PRT,BypassRatio,EtaT
    
    def PowBal(self,PRC,PwC,PwT,Nc,Nt,Step,PowRes):
        if PRC > self.C_PRC:
            Nc = Nc - Step
            return 1,Nc,Nt
        elif PRC < self.C_PRC-Step/10000:
            Nc = Nc + Step
            return 1,Nc,Nt
        else:
            if np.abs(Nc - Nt) < Step:
                return 0,Nc,Nt
            if Nc > Nt:
                Nt = Nt + Step
                return 1,Nc,Nt
            elif Nt > Nc:
                Nt = Nt - Step
                return 1,Nc,Nt
            else:
                print("Error:Power does not converge, change Step or PowRes")
                
#Method for one time iteration matching calculation
def Matching(Stgno,Nc):
    Stg = Stage(Stgno)
    print("Flow",Stg.C_Wt)
    Nt = Nc
    StgEff = 1 - Stg.S_MechLoss
    BackPres = Stg.T_Pout
    Cont = 1
    while Cont == 1:
        print("")
        PwC,Nc,P2c,T2c,PRC,EtaC = Stg.CompStage(1,Stg.T_Amb,Stg.P_Amb,Nc,Stg.C_flow,Stg.C_PRC) 
        Wt,P1t,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,Stg.CB_AFR,Stg.C_BLDR,Stg.T_BLDR)
        PwT,W_reqd,Nt,P2t,T2t,PRT,BPR,EtaT = Stg.TurbStage(1,StgEff,PwC,Stg.CB_t,P1t,Nt,Wt,BackPres) #Pin
        Cont,Nc,Nt = Stg.PowBal(PRC,PwC, PwT, Nc, Nt, Stg.S_Step,Stg.S_PowRes)
    print("")
    print('Power [kW]: C',PwC,'T',PwT,'RPM',Nc,Nt,'P1t[kPa]',P1t,'P2t',P2t,'P2c',P2c,'T2c[degC]',T2c-273.15,'T2t',T2t-273.15)
    print('Wt[kg/s]',Wt,'T1t_th[degC]',T1tth-273.15,'PRC',PRC,'PRT',PRT,'Turbine BPR',BPR,"Turbine Pressure Ratio",float(P2c/P1t)) 
    
    data = {
        "Compressor":{
        "Power[kW]":float(PwC),
        "rpm":float(Nc),
        "Outlet Pres":float(P2c),
        "Outlet Temp":float(T2c-273.15),
        "PRC":float(PRC),
        "Efficiency":float(EtaC),
        "Air Flow[kg/s]":float(Stg.C_flow),
        "Turbine Pressure Ratio":float(P2c/P1t)
        },
        "Turbine":{
        "Power[kW]":float(PwT),
        "rpm":float(Nc),
        "Inlet Pres.[kPa]":float(P1t),
        "Inlet Temp[K]":float(Stg.CB_t),
        "Outlet Pres[kPa]":float(P2t),
        "Outlet Temp[K]":float(T2t-273.15),
        "PRT":float(PRT),
        "Efficiency":float(EtaT),
        "Turbine Flow[kg/s]":float(Wt),
        "Bypass Ratio":float(BPR),
        },
    }

    toml_string = toml.dumps(data)
    output_file = "Result_" +Stg.name+".toml"
    with open(output_file,"w") as output:
        toml.dump(data,output)

Matching(1,60000)




"""

def PowBal(self,PRC,PwC,PwT,Nc,Nt,Step,PowRes,MechLoss):
    dPow = PwT - PwC
    dML = PwT*MechLoss
    print('dPow',dPow,'dML',dML)
    if PRC > self.C_PRC:
        Nc = Nc - Step
        return 1,Nc,Nt,dPow,dML
    elif PRC < self.C_PRC-Step/1000:
        Nc = Nc + Step
        return 1,Nc,Nt,dPow,dML
    else:
        if np.abs(Nc - Nt) < Step+1:
            return 0,Nc,Nt,dPow,dML
        if Nc > Nt:
            Nt = Nt + Step
            return 1,Nc,Nt,dPow,dML
        elif Nt > Nc:
            Nt = Nt - Step
            return 1,Nc,Nt,dPow,dML
        else:
            print("Error:Power does not converge, change Step or PowRes")
            

Option
Run = 1
dML = StgML - 2
while Run == 1:
    if dML < StgML:
        BackPres = BackPres*0.999
        Cont = 1
        while Cont ==1:
            PwC,Nc,P2c,T2c,PRC = Stg.CompStage(1,Stg.T_Amb,Stg.P_Amb,Nc,Stg.C_flow,Stg.C_PRC) 
            Wt,P1t,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,Stg.CB_AFR,Stg.C_BLDR,Stg.T_BLDR)
            PwT,Nt,P2t,T2t,PRT,BPR = Stg.TurbStage(1,Stg.CB_t,P1t,Nt,Wt,BackPres) #Pin
            Cont,Nc,Nt,dPow,dML = Stg.PowBal(PRC,PwC, PwT, Nc, Nt, Stg.S_Step,Stg.S_PowRes,StgML)
        Run = 1

    elif dML > StgML+1:
        BackPres = BackPres*1.001
        Cont = 1
        while Cont ==1:
            PwC,Nc,P2c,T2c,PRC = Stg.CompStage(1,Stg.T_Amb,Stg.P_Amb,Nc,Stg.C_flow,Stg.C_PRC) 
            Wt,P1t,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,Stg.CB_AFR,Stg.C_BLDR,Stg.T_BLDR)
            PwT,Nt,P2t,T2t,PRT,BPR = Stg.TurbStage(1,Stg.CB_t,P1t,Nt,Wt,BackPres) #Pin
            Cont,Nc,Nt,dPow,dML = Stg.PowBal(PRC,PwC, PwT, Nc, Nt, Stg.S_Step,Stg.S_PowRes,StgML)
        Run = 1
    else:
        Cont = 1
        while Cont ==1:
            PwC,Nc,P2c,T2c,PRC = Stg.CompStage(1,Stg.T_Amb,Stg.P_Amb,Nc,Stg.C_flow,Stg.C_PRC) 
            Wt,P1t,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,Stg.CB_AFR,Stg.C_BLDR,Stg.T_BLDR)
            PwT,Nt,P2t,T2t,PRT,BPR = Stg.TurbStage(1,Stg.CB_t,P1t,Nt,Wt,BackPres) #Pin
            Cont,Nc,Nt,dPow,dML = Stg.PowBal(PRC,PwC, PwT, Nc, Nt, Stg.S_Step,Stg.S_PowRes,StgML)
        print("")
        print('Power [kW]: C',PwC,'T',PwT,'RPM',Nc,Nt,'P1t[kPa]',P1t,'P2t',P2t,'P2c',P2c,'T2c[degC]',T2c-273.15,'T2t',T2t-273.15)
        print('Wt[kg/s]',Wt,'T1t_th[degC]',T1tth-273.15,'PRC',PRC,'PRT',PRT,'Turbine BPR',BPR,'Mech Power Loss[kW]',dPow)  
        Run = 0
"""
