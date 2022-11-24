"""
===-*- TurboCompressor Matching Tool -*-===
StageCalc Module

=====-*- General -*-=====
Copyright (c) makkiblog.com
MIT License 
coding: utf-8

===-*- VERSION -*-===
v0.1 Initial version
v0.2 Improved convergence
v1.0 Initial Release Working version
v1.1 Updated convergence parameter
v1.2 Fixed T2t error and updated turbine convergence

vvvCODEvvv
"""

import configparser
import numpy as np
import pandas
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
import tomli
import warnings
warnings.filterwarnings("ignore")

#if __name__ == '__main__':

#Class for each stage element calculation methods
class Stage:
    def __init__(self, stgno, settingfilename,params_file, target_file, input_mode,reload = False):
        setting_file = './Inputs/'+settingfilename+'.ini'
        setting = configparser.ConfigParser()
        setting.optionxform = str  #Case sensitve option
        setting.read(setting_file, encoding='utf8')
        initrpm = setting.getfloat("Settings", "Initial rpm")
        humidity = setting.getfloat("Environment Variables", "humidity [0-1]") #not used yet
        #matchingmode = setting.get("Settings", "Matching Mode (Bypass or AFR)") # matching mode Bypass or AFR correction
        analysisname = setting.get("Settings", "Analysis Name")
        params_file = setting.get("Settings", "Input filename")
        target_file = setting.get("Settings", "Target filename")
        input_mode = setting.get("Settings", "Input mode (csv or toml)")


        params_set = './Inputs/' + params_file + '.toml'
        with open(params_set, 'rb') as paramobj:
            params = tomli.load(paramobj)
        print("Input Geometry:\n",params)

        self.C_Wt=params['Compressor']['Weight']
        self.C_I=params['Compressor']['Inertia']
        self.C_gam=params['Compressor']['gamma']
        self.CB_LHV=params['Combustor']['FuelLatentHeatLHV']
        self.T_Wt=params['Turbine']['Weight']
        self.T_I=params['Turbine']['Inertia']
        self.S_Wt=params['Shaft']['Weight']
        self.S_I=params['Shaft']['Inertia']
        self.S_MechLoss=params['Shaft']['MechLoss']
        self.S_Step=params['Calculation']['Step']
        self.S_PowRes=params['Calculation']['PowRes']
        self.S_PRCRes=params['Calculation']['PRCRes']
        self.Wtres=params['Calculation']['WtRes']

        if input_mode == 'toml':
            target_set = './Inputs/' + target_file + '_' + str(stgno) + '.toml'
            with open(target_set, 'rb') as targobj:
                target = tomli.load(targobj)
            print("Stage Target Parameters:\n",target)            
        
            #Read Setting File (Each Stage Variables) toml file
            self.P_Amb = target['Compressor']['AmbientPressure']
            self.T_Amb = target['Compressor']['AmbientTemperature']
            self.C_Cp=(0.2208+51.48*10**-6*(self.T_Amb))*360/86
            self.C_flow=target['Compressor']['Flow']
            self.C_PRC=target['Compressor']['PRC']
            self.C_BLDR=target['Compressor']['BleedRatio']
            self.T_BLDR=target['Turbine']['BleedRatio']    

            self.CB_t=target['Combustor']['CombustionTemp_T1t']
            self.CB_dp=target['Combustor']['CombustordeltaP']
            self.CB_AFR=target['Combustor']['AirFuelRatio']        
            self.CB_Tempeff=target['Combustor']['TempEfficiency']
            self.T_Pout=target['Turbine']['BackPressure']

            self.T_Cp=(0.2208+51.48*10**-6*(self.CB_t+273.15))*360/86
            self.T_gam=1.42-84.6*10**-6*(self.CB_t+273.15) 

        elif input_mode == 'csv':
            target_set = './Inputs/' + target_file + '.csv'  #Read all cases at once
            with open(target_set, 'rb') as targobj:
                df_target = pandas.read_csv(targobj, sep =",",  index_col=0)
                index = str(stgno)
            print("Stage Target Parameters:\n",df_target.loc[index,:])
            self.P_Amb = float(df_target.loc[index,'Pc_in'])
            self.T_Amb = float(df_target.loc[index,'Tc_in'])
            self.C_Cp=(0.2208+51.48*10**-6*(self.T_Amb))*360/86
            self.C_flow = float(df_target.loc[index,'Wc_in'])
            self.C_PRC = float(df_target.loc[index,'PRC'])
            self.C_BLDR=float(df_target.loc[index,'Bleed_C'])
            self.T_BLDR=float(df_target.loc[index,'Bleed_T'])

            self.CB_t = float(df_target.loc[index,'T1t_th'])
            self.CB_dp = float(df_target.loc[index,'CB_dp'])
            self.CB_AFR = float(df_target.loc[index,'AFR']  )      
            self.CB_Tempeff = float(df_target.loc[index,'TempEff'])
            self.T_Pout = float(df_target.loc[index,'Pt_out'])

            self.T_Cp=(0.2208+51.48*10**-6*(self.CB_t+273.15))*360/86
            self.T_gam=1.42-84.6*10**-6*(self.CB_t+273.15) 

        else:
            print('setup.ini ERROR: Input mode must be either csv or toml')
            exit()

    #def showmap(self,stgno):

    def CompStage(self,stgno,cmap,Tin,Pin,N_C,W_C,PRC_Target):
        #Read Compressor Map
        df = pandas.read_csv("./Inputs/%s_%s.csv" % (cmap,1), sep =",")  #for single -stage
        #df = pandas.read_csv("./Inputs/%s_%s.csv" % (cmap,stgno), sep =",") #for multi-stage maps
        W_Cstar = W_C*(Tin/288.15)**0.5*101.325/Pin  #Corrected flow rate
        N_Cstar = N_C/(Tin/288.15)**0.5 #Corrected Speed
        if N_Cstar > max(df['NC*']):
            print ("ERROR:Initial rpm too large", N_Cstar)
            exit()
        elif N_Cstar < min(df['NC*']):
            print ("ERROR:Initial rpm too small", N_Cstar)
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
            #CLI Output for Debug
            #print('COMPRESSOR','PwC[kW]:%.1f|' % Power,'NC[rpm]:%.0f|' % N_C,'NC*:%.0f|' % N_Cstar,'WC[kg/s]:%.2f|' % W_C,'WC*:%.2f|' % W_Cstar,'EtaC:%.3f|' % EtaC,'PRC:%.3f|' % PRC)
            return Power,W_Cstar,N_C,Pout,Tout,PRC,EtaC

    def Combustor(self,W_C,Tin,Pin,AFR,C_BLDR,T_BLDR):
        W_Comb = W_C*(1-C_BLDR)
        W_F = W_Comb/AFR
        W_T = (W_Comb+W_F)*(1-T_BLDR)
        deltaT = self.CB_LHV*10**6*W_F/(self.C_Cp*10**3*W_T)
        Tout = Tin + deltaT
        Pout = Pin*(1-self.CB_dp)
        #CLI Output for Debug
        #print('COMBUSTOR','P2c[kPa]:%.1f|'% Pin,'P1t[kPa]:%.1f|' % Pout,'T1t[K]:%.1f|' % Tout)
        return W_F,W_T,Pout,Tout

    def TurbStage(self,stgno,tmap,EtaM,PwReq,Tin,Pininit,N_T,W_T,Pout):
        #Read Turbine Map
        W_Tstarinit = W_T*(Tin/288.15)**0.5*101.325/Pininit #Corrected flow rate
        N_Tstar = N_T/(Tin/288.15)**0.5 #Corrected Speed
        df = pandas.read_csv("./Inputs/%s_%s.csv" % (tmap,1), sep =",") #for single-stage
        #df = pandas.read_csv("./Inputs/%s_%s.csv" % (tmap,stgno), sep =",") #for multi-stage maps
                
        f_EtaT = interp2d(x = df['NT*'], y = df['PRT'], z = df['EtaT'])
        f_WT = interp2d(x = df['NT*'], y = df['PRT'], z = df['WT*'])
        #PRT = f_PRT(N_Tstar,W_Tstar)
        PRT = Pininit/Pout
        W_Tstar = f_WT(N_Tstar,PRT)
        EtaT = max(f_EtaT(N_Tstar,PRT),min(df['EtaT']))  #Avoid Extrapolation and divergence

        deltaW = W_Tstarinit - W_Tstar
        while abs(deltaW) > self.Wtres*min(df['WT*']):  #find WT* from PRT due to sensitivity issue
            PRT += deltaW/W_Tstar*self.Wtres
            Pin = PRT*Pout
            W_Tstar = f_WT(N_Tstar,PRT)
            EtaT = max(f_EtaT(N_Tstar,PRT),min(df['EtaT'])) #Avoid Extrapolation and divergence
            deltaW = W_Tstarinit - W_Tstar
            if PRT > max(df['PRT']): #Avoid Extrapolation and divergence
                PRT = max(df['PRT'])
                print('WARNING PRT Exceed limits, Max. PRT Used:',PRT)
                W_Tstar = max(df['WT*'])
                EtaT = min(df['EtaT'])
                break
            #print('deltaW Map:',deltaW,'deltaPRT:',deltaW/min(df['WT*']),'PRT:',PRT,'W_Tstar:',W_Tstar)
        Pin = PRT*Pout
        W_Tcalc = W_Tstar/((Tin/288.15)**0.5*101.325/Pin)
        Tadb = Tin*(1-PRT**((self.T_gam-1)/self.T_gam))  #deltaT
        Tout = Tadb + Tin
        Enth = -self.T_Cp*Tadb*EtaT*W_Tcalc
        PwAct = PwReq/EtaM
        W_reqd = -PwAct/(self.T_Cp*Tadb*EtaT)
        W_reqdstar = W_reqd*(Tin/288.15)**0.5*101.325/Pin
        BypassRatio = (W_Tcalc - W_reqd)/W_Tcalc
        Pini = Pout*PRT
        
        #print('W_Tcalc',W_Tcalc,'W_reqd',W_reqd,'BPR',BypassRatio,'PwAct',PwAct,'Enth',Enth,'EtaT',EtaT,'Tadb',Tadb)
        
        #CLI Output for Debug
        #print('TURBINE','Tot. Enthal[kW]:%.1f|' % Enth,'PwT Reqd[kW]:%.1f|' % PwAct,'NT[rpm]:%.0f|' % N_T,'NT*:%.0f|' % N_Tstar,'W_Tot[kg/s]:%.2f|' % W_T,'W_Tot*:%.2f|' % W_Tstar,'W_Reqd Turb[kg/s]:%.2f|' % W_reqd,'W_Req Turb*:%.2f|' % W_reqdstar,'EtaT:%.3f|' % EtaT,'PRT:%.3f|' % PRT,'BPR:%.2f|' % BypassRatio,'P1t Actual[kPaA]:%.1f|' % Pini,'P2t Actual[kPaA]:%.1f|' % Pout)
        if BypassRatio < 0:
            #print("ERROR: Bypass Ratio must be larger than 0.  BPR:",BypassRatio)
            return Enth,W_reqd,W_Tcalc,W_Tstar,PwAct,N_T,Pout,Pini,Tout,PRT,BypassRatio,EtaT
        else:
            return Enth,W_reqd,W_Tcalc,W_Tstar,PwAct,N_T,Pout,Pini,Tout,PRT,BypassRatio,EtaT

    def PowBal(self,PRC,Nc,Nt,Step): #Compressor side power matching
        #Convergence by compressor scaling law
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

            #Optional Nt convergence mode for slip coupling
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
