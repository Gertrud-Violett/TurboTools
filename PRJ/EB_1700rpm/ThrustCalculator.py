# -*- coding: utf-8 -*-
"""
2 Plane Thrust Calculator
MIT License 2022 K.T.
"""
import sys
import os
import json
import csv
import math
import configparser
import numpy as np
import matplotlib.pyplot as plt
import importlib
import pandas as pd

from scipy import interpolate
from scipy import optimize

#Gas Constants
kappa = 1.4
R = 287.058


#Setup (Visual)
import seaborn as sns
sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 1.5})
sns.set_style('whitegrid')

plt.rcParams["figure.figsize"] = (16,9)
plt.rcParams['figure.facecolor'] = 'white'


#Read Result file
def csv_df(input):
    return pd.read_csv(input, skiprows=2, header=0,index_col=False)

#Read Setting file
def readset(filename): #Syntax:"setting filename.ini"

    setting_file = filename
    setting = configparser.ConfigParser()
    setting.optionxform = str #Case Sensitive
    setup = setting.read(setting_file, encoding='utf8')
    
    name = setting.get("Settings", "Analysis Name")
    csvfile = setting.get("Settings", "CSV file")
    pressure = setting.get("Settings", "Gauge or Abs[kPa]")
    TipOp = setting.get("Settings", "Tip Pressure (Yes/No)")
    Shroud = setting.get("Settings", "Shroud Pressure (Yes/No)")
    C_TipDia = setting.getfloat("Compressor Geometry", "Tip Diameter[mm]")
    C_BDDia = setting.getfloat("Compressor Geometry", "Backdisk Diameter[mm]")
    C_SealDia = setting.getfloat("Compressor Geometry", "BackSeal Diameter[mm]")
    C_BPisDia = setting.getfloat("Compressor Geometry", "Balance Piston Diameter[mm]")
    C_InducDia = setting.getfloat("Compressor Geometry", "Inducer Diameter[mm]")
    C_NoseDia = setting.getfloat("Compressor Geometry", "Nose Diameter[mm]")
    C_H_ExitBlade = setting.getfloat("Compressor Geometry", "Exducer Blade Height[mm]")
    C_BladeNo = setting.getfloat("Compressor Geometry", "No. of Blades[mm]")
    C_beta = setting.getfloat("Compressor Geometry", "Blade Exit Angle Beta[deg]")

    T_TipDia = setting.getfloat("Turbine Geometry", "Tip Diameter[mm]")
    T_BDDia = setting.getfloat("Turbine Geometry", "Backdisk Diameter[mm]")
    T_SealDia = setting.getfloat("Turbine Geometry", "BackSeal Diameter[mm]")
    T_BPisDia = setting.getfloat("Turbine Geometry", "Balance Piston Diameter[mm]")
    T_ExducDia = setting.getfloat("Turbine Geometry", "Exducer Diameter[mm]")
    T_NoseDia = setting.getfloat("Turbine Geometry", "Nose Diameter[mm]")
    T_H_InducBlade = setting.getfloat("Turbine Geometry", "Inducer Blade Height[mm]")
    T_BladeNo = setting.getfloat("Turbine Geometry", "No. of Blades[mm]")

    Setting_Dict={'name' : name, 'csv' : csv}
    C_Dict={'C_TipDia': C_TipDia, 'C_BDDia' : C_BDDia, 'C_SealDia' : C_SealDia, 'C_BPisDia' : C_BPisDia, 'C_InducDia' : C_InducDia, 'C_NoseDia' : C_NoseDia, 'C_H_ExitBlade' : C_H_ExitBlade, 'C_BladeNo' : C_BladeNo, 'C_beta': C_beta}
    T_Dict={'T_TipDia' : T_TipDia, 'T_BDDia' : T_BDDia, 'T_SealDia' : T_SealDia, 'T_BPisDia' : T_BPisDia, 'T_ExducDia' : T_ExducDia, 'T_NoseDia' : T_NoseDia, 'T_H_InducBlade' : T_H_InducBlade, 'T_BladeNo' : T_BladeNo}
    return(name,C_Dict,T_Dict,csvfile,pressure,TipOp,Shroud)

#Read and create pd
name, C_Spec, T_Spec, csvfile, pressure,TipOp,Shroud= readset("settings.ini")
res_df = csv_df(csvfile)
df = res_df #create df copy for calculation

#Update data is Absolute pressure is used
while True:
    if pressure =="Abs":
        df['Pamb'] = 0
        print("Absolute Pressure")
        break
    elif pressure=="Gauge":
        print("Gauge Pressure")
        break
    else:
        print("Specify Measurement Units in Gauge or Abs")
        break

#Calculate and update properties for Wheels [metric units]
C_Spec['A_Ind']=C_Spec['C_InducDia']**2*np.pi/4*10**-6 #Includes nose
C_Spec['A_Nose']=C_Spec['C_NoseDia']**2*np.pi/4*10**-6 
C_Spec['A_Exd']=(C_Spec['C_TipDia']**2 - C_Spec['C_BDDia']**2)*np.pi/4*10**-6
C_Spec['A_ExitRad']=np.pi/2*(C_Spec['C_TipDia']+C_Spec['C_BDDia'])*np.sqrt((C_Spec['C_TipDia']-C_Spec['C_BDDia'])**2/4+C_Spec['C_H_ExitBlade']**2)*10**-6
C_Spec['A_BD']=(C_Spec['C_BDDia']**2 - C_Spec['C_SealDia']**2)*np.pi/4*10**-6
C_Spec['A_Sh']=(C_Spec['C_TipDia']**2 - C_Spec['C_InducDia']**2)*np.pi/4*10**-6
C_Spec['W_Exd']=(C_Spec['C_TipDia'] + C_Spec['C_BDDia'])*np.pi/2*C_Spec['C_H_ExitBlade']*10**-6

T_Spec['A_Ind']=(T_Spec['T_TipDia']**2 - T_Spec['T_BDDia']**2)*np.pi/4*10**-6
T_Spec['A_Exd']=T_Spec['T_ExducDia']**2*np.pi/4*10**-6
T_Spec['A_BD']=(T_Spec['T_BDDia']**2 - T_Spec['T_SealDia']**2)*np.pi/4*10**-6
T_Spec['A_Sh']=(T_Spec['T_TipDia']**2 - T_Spec['T_ExducDia']**2)*np.pi/4*10**-6


#Calculations & append to pandas df
df['P1Cabs']=df['Pamb']+df['P1C']
df['P2Cabs']=df['Pamb']+df['P2C']
df['P1Tabs']=df['Pamb']+df['P1T']
df['P2Tabs']=df['Pamb']+df['P2T']
df['PbCinabs']=df['Pamb']+df['PbCin']
df['PbCoutabs']=df['Pamb']+df['PbCout']
df['PbTinabs']=df['Pamb']+df['PbTin']
df['PbToutabs']=df['Pamb']+df['PbTout']
df['T1Cabs']=df['T1C']+273.15
df['T2Cabs']=df['T2C']+273.15
df['T1Tabs']=df['T1T']+273.15
df['T2Tabs']=df['T2T']+273.15
df['PRC']=df['P2Cabs']/df['P1Cabs']
df['PRT']=df['P2Tabs']/df['P1Tabs']

#Gas Calcs
df['rho1c']=df['P1Cabs']/df['T1Cabs']/R*10**3
df['rho2c']=df['P2Cabs']/df['T2Cabs']/R*10**3 #??
df['Sonic2c']=np.sqrt(df['T2Cabs']*R*kappa)
df['AirFlowkgs']=df['AirFlow']/3600

#Compressor Specs
df['etaC']=df['PRC']/((df['T2Cabs']/df['T1Cabs'])**(kappa/(kappa-1)))
df['v1c(C1)']=df['AirFlow']/3600/(C_Spec['A_Ind']-C_Spec['A_Nose'])/df['rho1c']
df['U2']=np.pi*df['TCrpm']/30*C_Spec['C_TipDia']/2*10**-3
df['slip']=1/(1+3.6/((1-(C_Spec['C_InducDia']/C_Spec['C_TipDia'])**2)*C_Spec['C_BladeNo']))
df['Mach']=0.3
df['Qt']=df['AirFlow']+df['FuelFlow']

print(df['P1Tabs'])

#Derive Mach no. and Compressor Pressures calculation
def CompCalc(M0,T2Cabs,T2Tabs,PS1C,PS2C,PS1T,PS2T,PRT,slip,flow,U2,Qt):
    deltaM = 0.9
    Mach_I=M0
    while deltaM > 0.001:
        TS2C=T2Cabs/(1+0.4/2*Mach_I**2)
        Sonic2c=np.sqrt(TS2C*R*kappa)
        rhoS2c=PS2C*10**3/TS2C/R
        Cr2=flow/3600/rhoS2c/C_Spec['A_ExitRad']
        phi2=Cr2/U2
        wt2=Cr2*np.tan(C_Spec['C_beta']*np.pi/180)
        Ct2=(U2-wt2)*slip
        v2c=np.sqrt(Ct2**2+Cr2**2) #Exducer Velocity
        Mach_E=v2c/Sonic2c
        deltaM = np.abs(Mach_E-Mach_I)
        Mach_I=Mach_E
    PT2C=PS2C+rhoS2c/2*v2c**2*10**-3
    PtC=PT2C/((1+(kappa-1)/2*Mach_E**2)**(kappa/(kappa-1)))
    PtT=PS1T*(1-1.1*(1-PRT**(0.33/1.33)))
    PshC=(PS1C+PtC)/2
    PshT=(PS2T+PtT)/2
    rhoS2t=PS2T/T2Tabs/R
    v2t=Qt/3600/T_Spec['A_Exd']

#   visc=0.000001458/(110.4+T2Cabs)/rhoS2c*T2Cabs**(3/2)
    return Mach_E,v2c,v2t,Cr2,wt2,Ct2,phi2,PtC,PtT,PshC,PshT

i=0
while i <= len(df.index)-1:
    df.loc[i,'Mach'], df.loc[i,'v2c(C2)'], df.loc[i,'v2t'],df.loc[i,'cr2'], df.loc[i,'wt2'],df.loc[i,'ct2'],df.loc[i,'phi2'], df.loc[i,'PtC_Cal'], df.loc[i,'PtT_Cal'], df.loc[i,'PshC_Cal'], df.loc[i,'PshT_Cal'] = CompCalc(df.loc[i,'Mach'],df.loc[i,'T2Cabs'],df.loc[i,'T2Tabs'],df.loc[i,'P1Cabs'],df.loc[i,'P2Cabs'],df.loc[i,'P1Tabs'],df.loc[i,'P2Tabs'],df.loc[i,'PRT'],df.loc[i,'slip'],df.loc[i,'AirFlow'],df.loc[i,'U2'],df.loc[i,'Qt'])
    i+=1

#Replace Tip and Shroud measurement with calculated values
while True:
    if TipOp =="No":
        df['PtCabs'] = df['PtC_Cal']
        df['PtTabs'] = df['PtT_Cal']
        print("Calculated Tip Pressure Used")
        if Shroud =="No":
            df['PshCabs'] = (df['PtCabs']+df['P1Cabs'])/2
            df['PshTabs'] = (df['PtTabs']+df['P2Tabs'])/2
            print("Calculated Shroud Pressure Used")
            break
        elif Shroud == "Yes":
            df['PshCabs'] = df['PshC']
            df['PshTabs'] = df['PshT']
            print("Measured Shroud Pressure Used")
            break
        else:
            print("Specify Yes or No")
            break
    elif TipOp=="Yes":
        df['PtCabs'] = df['Pamb']+df['PtC']
        df['PtTabs'] = df['Pamb']+df['PtT']
        print("Measured Tip Pressure Used")
        if Shroud =="No":
            df['PshCabs'] = (df['PtCabs']+df['P1Cabs'])/2
            df['PshTabs'] = (df['PtTabs']+df['P2Tabs'])/2
            print("Calculated Shroud Pressure Used")
            break
        elif Shroud == "Yes":
            df['PshCabs'] = df['PshC']
            df['PshTabs'] = df['PshT']
            print("Measured Shroud Pressure Used")
            break
        else:
            print("Specify Yes or No")
            break
    else:
        print("Specify Yes or No")
        break


#Thrust Calculation
df['RESULTS']="|==>"
df['Th1C_Stat_Inducer']=df['P1C']*10**3*C_Spec['A_Ind']*-1  #CompPush
df['Th1C_Dyn_Inducer']=df['v1c(C1)']*df['AirFlow']/3600*-1  #CompPush
df['ThshC_Shroud']=(df['PshCabs']-df['Pamb'])*10**3*C_Spec['A_Sh']*-1  #CompPush
df['ThbC_Backdisk']=((df['PbCinabs']+df['PbCoutabs'])/2-df['Pamb'])*10**3*C_Spec['A_BD']  #CompPull
df['Th2C_Dyn_Exducer']=(df['PtCabs']-df['Pamb'])*10**3*C_Spec['A_Exd']  #CompPull

df['Th2T_Stat_Exducer']=df['P2T']*10**3*T_Spec['A_Exd']  #CompPull
df['Th2T_Dyn_Exducer']=df['v2t']*df['Qt']/3600 #CompPull
df['ThshT_Shroud']=(df['PshTabs']-df['Pamb'])*10**3*T_Spec['A_Sh'] #CompPull
df['ThbT_BackDisk']=((df['PbTinabs']+df['PbToutabs'])/2-df['Pamb'])*10**3*T_Spec['A_BD']*-1 #CompPush
df['Th1T_Inducer']=(df['PtTabs']-df['Pamb'])*10**3*T_Spec['A_Ind']*-1  #CompPush
df['TotalThrust(CompPull+)']=df['Th1C_Stat_Inducer']+df['Th1C_Dyn_Inducer']+df['ThshC_Shroud']+df['ThbC_Backdisk']+df['Th2C_Dyn_Exducer']+df['Th2T_Stat_Exducer']+df['Th2T_Dyn_Exducer']+df['ThshT_Shroud']+df['ThbT_BackDisk']+df['Th1T_Inducer']

df['deltP2cP1t']=-(df['P2C']-df['P1T'])
df['deltP2tP1t']=-(df['P2T']-df['P1T'])
#print(df)
#print(C_Spec)
#print(T_Spec)
outdata=[{"Comp.":(C_Spec)}, {"Turb.":(T_Spec)}]
df.to_csv(name + '_res.csv')
with open(name + 'geom.json', 'w') as fp:
    json.dump(outdata, fp)

#plot results
plt.clf()
plt.figure()
df.plot(x='TCrpm', y=['P1Tabs', 'P2Tabs','P1Cabs','P2Cabs','PtCabs','PtTabs','deltP2cP1t'])
plt.yticks(np.linspace(0,600,20))
plt.ylabel('Pressure [kPaA]')
plt.xlabel('rpm')
plt.savefig(name + '_RPM-Pressure.png')
#plt.show()

plt.clf()
plt.figure()
df.plot(x='Time', y=['P1Tabs', 'P2Tabs','P1Cabs','P2Cabs','PtCabs','PtTabs','deltP2cP1t'])
plt.yticks(np.linspace(0,600,20))
plt.ylabel('Pressure [kPaA]')
plt.xlabel('Time[s]')
plt.savefig(name + '_Time-Pressure.png')
#plt.show()

plt.clf()
plt.figure()
df.plot(x='TCrpm', y=['TotalThrust(CompPull+)'])
plt.yticks(np.linspace(-100,250,20))
plt.ylabel('Thrust [N]')
plt.xlabel('rpm')
plt.savefig(name + '_RPM-Thrust.png')
#plt.show()

plt.clf()
plt.figure()
df.plot(x='Time', y=['TotalThrust(CompPull+)'])
plt.yticks(np.linspace(-100,250,20))
plt.ylabel('Thrust [N]')
plt.xlabel('Time[s]')
plt.savefig(name + '_Time-Thrust.png')
#plt.show()

plt.clf()
plt.figure()
df.plot(x='Time', y=['Poil'])
#plt.yticks(np.linspace(-100,250,20))
plt.ylabel('Oil Pressure [bar]')
plt.xlabel('Time[s]')
plt.savefig(name + '_Time-OilPres.png')
#plt.show()

plt.clf()
plt.figure()
df.plot(x='Time', y=['Toil'])
#plt.yticks(np.linspace(-100,250,20))
plt.ylabel('Oil Temp [degC]')
plt.xlabel('Time[s]')
plt.savefig(name + '_Time-OilTemp.png')
#plt.show()

plt.clf()
plt.figure()
df.plot(x='Time', y=['T1C','T1T','T2C','T2T'])
#plt.yticks(np.linspace(-100,250,20))
plt.ylabel('Oil Temp [degC]')
plt.xlabel('Time[s]')
plt.savefig(name + '_Time-GasTemp.png')
#plt.show()