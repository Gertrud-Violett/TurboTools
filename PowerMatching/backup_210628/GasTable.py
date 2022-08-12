"""
Gas Properties

"""

import sys
import os
os.environ['RPPREFIX'] = r'.'
import math
import numpy as np
import pandas
import matplotlib.pyplot as plt
import importlib
from scipy import interpolate
from scipy import optimize
from scipy.optimize import newton
import configparser

#from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
#%matplotlib inline
#%config InlineBackend.figure_formats = {'png', 'retina'}
plt.rcParams['figure.figsize'] = [8, 6]


#########WIP for Refprop function################

#RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])

# 明示的にパスを指定する。ライブラリのインスタンス化する。
#RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
# 次にREFPROPのインスタンスにルートディレクトリがどこか指示。 
# このルートディレクトリには、少なくとも次のものを含める必要があります。
# A）REFPRP64.DLL（32ビットウィンドウの場合はREFPROP.dll、LinuxまたはOSXの場合はそれぞれlibrefprop.soまたはlibrefprop.dylib）
# B）FLUIDSフォルダ（大文字と小文字を区別）
# C）MIXTURESフォルダ（大文字と小文字を区別）
#RP.SETPATHdll(os.environ['RPPREFIX'])
#MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum
#MASS_BASE_SI = RP.GETENUMdll(0, "MASS BASE SI").iEnum

#p_Pa = 101325
#Q = 0.0
#r = RP.REFPROPdll("Water","PQ","T",MOLAR_BASE_SI,0,0,p_Pa,Q,[1.0])
#r.Output[0]

#https://nbviewer.jupyter.org/gist/ina111/a4d9507eef905c5aeb11fffd42d32a48

#by absolute 
#r = RP.REFPROPdll(os.path.join(os.environ['RPPREFIX'],"FLUIDS","WATER.FLD"),"PQ","T",MOLAR_BASE_SI,0,0,p_Pa,Q,[1.0])





#Output Gas Properties
#Density rho,Cp,Cv,Sound Speed Mach, Thermal Conductivity kappa, Viscosity mu
#NIST Table: 0:Temperature (C)    1:Pressure (MPa)  2:Density (kg/m3) 3:Volume (m3/kg) 4:Internal Energy (kJ/kg) 5:Enthalpy (kJ/kg) 6:Entropy (J/g*K) 7:Cv (J/g*K) 8:Cp (J/g*K)  9:Sound Spd. (m/s) 10:Joule-Thomson (K/MPa)   11:Viscosity (Pa*s)    12:Therm. Cond. (W/m*K)    13:Phase    

def N2_GasProp(Temp,Pres):
    if Pres < 200:
        P = 0.1
    elif Pres < 300:
        P = 0.2
    elif Pres < 600:
        P = 0.5
    elif Pres < 1300:
        P = 1.0
    elif Pres < 1600:
        P = 1.5
    elif Pres < 2200:
        P = 2.0
    elif Pres < 3300:
        P = 3.0
    else:
        print("Pressure Out of Table Range")
    df = pandas.read_csv("./NIST_Table/N2_%sMPa.txt" % (P), sep ="\t")
    T = df["Temperature (C)"]
    index = abs(T - Temp).idxmin()
    Tp = df.iloc[index,0]
    rho = df.iloc[index,2]
    Cv = df.iloc[index,7]
    Cp = df.iloc[index,8]
    mu = df.iloc[index,11]
    gamma = Cp/Cv
    return Tp,rho,Cv,Cp,gamma,mu 

A = N2_GasProp(400,0.1)
print(A)


def O2_GasProp(Temp,Pres):
    if Pres < 200:
        P = 0.1
    elif Pres < 300:
        P = 0.2
    elif Pres < 600:
        P = 0.5
    elif Pres < 1300:
        P = 1.0
    elif Pres < 1600:
        P = 1.5
    elif Pres < 2200:
        P = 2.0
    elif Pres < 3300:
        P = 3.0
    else:
        print("Pressure Out of Table Range")
    df = pandas.read_csv("./NIST_Table/O2_%sMPa.txt" % (P), sep ="\t")
    T = df["Temperature (C)"]
    index = abs(T - Temp).idxmin()
    Tp = df.iloc[index,0]
    rho = df.iloc[index,2]
    Cv = df.iloc[index,7]
    Cp = df.iloc[index,8]
    mu = df.iloc[index,11]
    gamma = Cp/Cv
    return Tp,rho,Cv,Cp,gamma,mu  


def CO2_GasProp(Temp,Pres):
    hoge = Temp+Pres 
    return hoge

def H2O_GasProp(Temp,Pres):
    hoge = Temp+Pres

def Ar_GasProp(Temp,Pres):
    hoge = Temp+Pres


#Calculate Mixture properties
#Mixture Ratio:MR Must total 1. Others_MR compensate for missing MR as Nitrogen
R = 8.314510


def AirProp(Temp,Pres,N2_MF,O2_MF,CO2_MF,humidity,Ar_MF):
    Ps = humidity*math.exp(34.494 - 4924.99/(Temp-273.15+237.1))/((Temp-273.15+105)**1.57)
    H2O_MF = Ps/(Pres*10**3)
    Others_MF = 1-N2_MF-O2_MF-CO2_MF-Ar_MF
    N2_MF+=Others_MF
    O2 = O2_GasProp(Temp,Pres)

    #https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml Air Saturation Pressure Formula for Water
    return H2O_MF, Ps  
 #   return rho,Cp,Cv,gamma,mu,H2O_MF,Ps

def GasProp(Temp,Pres,N2_MF,O2_MF,CO2_MF,H20_MF,Ar_MF):
    Others_MF = 1- N2_MF - O2_MF - CO2_MF - H2O_MF - Ar_MF
    gamma = gamma_N2*N2_MF+gamma_O2*O2_MF+gamma_N2*Others_MF
    rho = 1 / (N2_MF/rho_N2 +O2_MF/rho_O2 + CO2_MF/rho_CO2 + H2O_MF/rho_H2O + Ar_MF/rho_Ar + Others_MF/rho_N2)
    return(rho,Cp,Cv,gamma,kappa,mu) 



print(AirProp(300,100,0.78,0.21,0,1,0))