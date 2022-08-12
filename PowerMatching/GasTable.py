"""
Gas Properties Calculator via NIST Gas Properties Table

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
plt.rcParams['figure.figsize'] = [8, 6]


#input file
setting_file = 'MultiStage_Matching_input.ini'
setting = configparser.ConfigParser()
setting.optionxform = str  # 大文字小文字を区別するおまじない
inputparams = setting.read(setting_file, encoding='utf8')

name = setting.get("Settings", "Analysis Name")
T_Amb = setting.getfloat("Environment Variables", "Ambient Pressure [kPaA]")
P_Amb = setting.getfloat("Environment Variables", "Ambient Temperature [K]")

#Output Gas Properties
#INPUT SYNTAX: Temp:K, Pres: kPa
#NIST Table Parameters: 0:Temperature (C)    1:Pressure (MPa)  2:Density (kg/m3) 3:Volume (m3/kg) 4:Internal Energy (kJ/kg) 5:Enthalpy (kJ/kg) 6:Entropy (J/g*K) 7:Cv (J/g*K) 8:Cp (J/g*K)  9:Sound Spd. (m/s) 10:Joule-Thomson (K/MPa)   11:Viscosity (Pa*s)    12:Therm. Cond. (W/m*K)    13:Phase (Name)
#WARNING FOR USAGE, LIMIT: Above 1700degC is calculated as 1700degC!

def N2_GasProp(Temp,Pres):
    if Pres < 200:
        P = 0.1 
        Pi = 0.2
    elif Pres < 500:
        P = 0.2 
        Pi = 0.5
    elif Pres < 1000:
        P = 0.5 
        Pi = 1.0
    elif Pres < 1500:
        P = 1.0 
        Pi = 1.5
    elif Pres < 2000:
        P = 1.5 
        Pi = 2.0
    elif Pres < 3000:
        P = 2.0 
        Pi = 3.0
    elif Pres < 3300:
        P = 3.0
    else:
        print("Pressure Out of Table Range")
    df = pandas.read_csv("./NIST_Table/N2_%sMPa.txt" % (P), sep ="\t")
    dfi = pandas.read_csv("./NIST_Table/N2_%sMPa.txt" % (Pi), sep ="\t")
    T = df["Temperature (C)"]
    Ti = dfi["Temperature (C)"]
    index = abs(T - (Temp-273.15)).idxmin()
    indexi = abs(Ti - (Temp-273.15)).idxmin()
    Tp = df.iloc[index,0]
    #Linear Approximation
    rho = (Pres/P/10**3-1)*(dfi.iloc[indexi,2]-df.iloc[index,2])+df.iloc[index,2]
    Cv = df.iloc[index,7]
    Cp = df.iloc[index,8]
    mu = df.iloc[index,11]
    gamma = Cp/Cv
    return Tp,rho,Cv,Cp,gamma,mu ,P,Pi

def O2_GasProp(Temp,Pres):
    if Pres < 200:
        P = 0.1 
        Pi = 0.2
    elif Pres < 500:
        P = 0.2 
        Pi = 0.5
    elif Pres < 1000:
        P = 0.5 
        Pi = 1.0
    elif Pres < 1500:
        P = 1.0 
        Pi = 1.5
    elif Pres < 2000:
        P = 1.5 
        Pi = 2.0
    elif Pres < 3000:
        P = 2.0 
        Pi = 3.0
    elif Pres < 3300:
        P = 3.0
    else:
        print("Pressure Out of Table Range")
    df = pandas.read_csv("./NIST_Table/O2_%sMPa.txt" % (P), sep ="\t")
    dfi = pandas.read_csv("./NIST_Table/O2_%sMPa.txt" % (Pi), sep ="\t")
    T = df["Temperature (C)"]
    Ti = dfi["Temperature (C)"]
    index = abs(T - (Temp-273.15)).idxmin()
    indexi = abs(Ti - (Temp-273.15)).idxmin()
    Tp = df.iloc[index,0]
    #Linear Approximation
    rho = (Pres/P/10**3-1)*(dfi.iloc[indexi,2]-df.iloc[index,2])+df.iloc[index,2]
    Cv = df.iloc[index,7]
    Cp = df.iloc[index,8]
    mu = df.iloc[index,11]
    gamma = Cp/Cv
    return Tp,rho,Cv,Cp,gamma,mu  


def CO2_GasProp(Temp,Pres):
    if Pres < 200:
        P = 0.1 
        Pi = 0.2
    elif Pres < 500:
        P = 0.2 
        Pi = 0.5
    elif Pres < 1000:
        P = 0.5 
        Pi = 1.0
    elif Pres < 1500:
        P = 1.0 
        Pi = 1.5
    elif Pres < 2000:
        P = 1.5 
        Pi = 2.0
    elif Pres < 3000:
        P = 2.0 
        Pi = 3.0
    elif Pres < 3300:
        P = 3.0
    else:
        print("Pressure Out of Table Range")
    df = pandas.read_csv("./NIST_Table/CO2_%sMPa.txt" % (P), sep ="\t")
    dfi = pandas.read_csv("./NIST_Table/CO2_%sMPa.txt" % (Pi), sep ="\t")
    T = df["Temperature (C)"]
    Ti = dfi["Temperature (C)"]
    index = abs(T - (Temp-273.15)).idxmin()
    indexi = abs(Ti - (Temp-273.15)).idxmin()
    Tp = df.iloc[index,0]
    #Linear Approximation
    rho = (Pres/P/10**3-1)*(dfi.iloc[indexi,2]-df.iloc[index,2])+df.iloc[index,2]
    Cv = df.iloc[index,7]
    Cp = df.iloc[index,8]
    mu = df.iloc[index,11]
    gamma = Cp/Cv
    return Tp,rho,Cv,Cp,gamma,mu 

def H2O_GasProp(Temp,Pres):
    if Pres < 200:
        P = 0.1 
        Pi = 0.2
    elif Pres < 500:
        P = 0.2 
        Pi = 0.5
    elif Pres < 1000:
        P = 0.5 
        Pi = 1.0
    elif Pres < 1500:
        P = 1.0 
        Pi = 1.5
    elif Pres < 2000:
        P = 1.5 
        Pi = 2.0
    elif Pres < 3000:
        P = 2.0 
        Pi = 3.0
    elif Pres < 3300:
        P = 3.0
    else:
        print("Pressure Out of Table Range")
    df = pandas.read_csv("./NIST_Table/H2O_%sMPa.txt" % (P), sep ="\t")
    dfi = pandas.read_csv("./NIST_Table/H2O_%sMPa.txt" % (P), sep ="\t")
    T = df["Temperature (C)"]
    Ti = dfi["Temperature (C)"]
    index = abs(T - (Temp-273.15)).idxmin()
    indexi = abs(Ti - (Temp-273.15)).idxmin()
    Tp = df.iloc[index,0]
    #Linear Approximation
    rho = (Pres/P/10**3-1)*(dfi.iloc[indexi,2]-df.iloc[index,2])+df.iloc[index,2]
    Cv = df.iloc[index,7]
    Cp = df.iloc[index,8]
    mu = df.iloc[index,11]
    gamma = Cp/Cv
    return Tp,rho,Cv,Cp,gamma,mu 

def Ar_GasProp(Temp,Pres):
    if Pres < 200:
        P = 0.1 
        Pi = 0.2
    elif Pres < 500:
        P = 0.2 
        Pi = 0.5
    elif Pres < 1000:
        P = 0.5 
        Pi = 1.0
    elif Pres < 1500:
        P = 1.0 
        Pi = 1.5
    elif Pres < 2000:
        P = 1.5 
        Pi = 2.0
    elif Pres < 3000:
        P = 2.0 
        Pi = 3.0
    elif Pres < 3300:
        P = 3.0
    else:
        print("Pressure Out of Table Range")
    df = pandas.read_csv("./NIST_Table/Ar_%sMPa.txt" % (P), sep ="\t")
    dfi = pandas.read_csv("./NIST_Table/Ar_%sMPa.txt" % (Pi), sep ="\t")
    T = df["Temperature (C)"]
    Ti = dfi["Temperature (C)"]
    index = abs(T - (Temp-273.15)).idxmin()
    indexi = abs(Ti - (Temp-273.15)).idxmin()
    Tp = df.iloc[index,0]
    #Linear Approximation
    rho = (Pres/P/10**3-1)*(dfi.iloc[indexi,2]-df.iloc[index,2])+df.iloc[index,2]
    Cv = df.iloc[index,7]
    Cp = df.iloc[index,8]
    mu = df.iloc[index,11]
    gamma = Cp/Cv
    return Tp,rho,Cv,Cp,gamma,mu 

def CO_GasProp(Temp,Pres):
    if Pres < 200:
        P = 0.1 
        Pi = 0.2
    elif Pres < 500:
        P = 0.2 
        Pi = 0.5
    elif Pres < 1000:
        P = 0.5 
        Pi = 1.0
    elif Pres < 1500:
        P = 1.0 
        Pi = 1.5
    elif Pres < 2000:
        P = 1.5 
        Pi = 2.0
    elif Pres < 3000:
        P = 2.0 
        Pi = 3.0
    elif Pres < 3300:
        P = 3.0
    else:
        print("Pressure Out of Table Range")
    df = pandas.read_csv("./NIST_Table/CO_%sMPa.txt" % (P), sep ="\t")
    dfi = pandas.read_csv("./NIST_Table/CO_%sMPa.txt" % (Pi), sep ="\t")
    T = df["Temperature (C)"]
    Ti = dfi["Temperature (C)"]
    index = abs(T - (Temp-273.15)).idxmin()
    indexi = abs(Ti - (Temp-273.15)).idxmin()
    Tp = df.iloc[index,0]
    #Linear Approximation
    rho = (Pres/P/10**3-1)*(dfi.iloc[indexi,2]-df.iloc[index,2])+df.iloc[index,2]
    Cv = df.iloc[index,7]
    Cp = df.iloc[index,8]
    mu = df.iloc[index,11]
    gamma = Cp/Cv
    return Tp,rho,Cv,Cp,gamma,mu 


#Calculate Mixture properties
#Mixture Ratio:MR Must total 1. Others_MR compensate for missing MR as Nitrogen
R = 8.314510


def AirProp(Temp,Pres,N2_MF,O2_MF,CO2_MF,humidity,Ar_MF):
    #Calculate water content H2O_MF based on humidity and steam saturation
    #https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice
    #H2O_MF: % Water Vapor in Weight, Ps: Molecuar Pressure Fraction of Water Vapor
    Ps = humidity*math.exp(34.494 - 4924.99/(Temp-273.15+237.1))/((Temp-273.15+105)**1.57)
    H2O_MF = Ps/(Pres*10**3)
    Others_MF = 1-N2_MF-O2_MF-CO2_MF-Ar_MF-H2O_MF
    N2_MF+=Others_MF
    O2 = O2_GasProp(Temp,Pres)
    N2 = N2_GasProp(Temp,Pres)
    H2O = H2O_GasProp(Temp,Pres)
    CO2 = CO2_GasProp(Temp,Pres)
    Ar = Ar_GasProp(Temp,Pres)
    rho = 1/(O2_MF/O2[1]+N2_MF/N2[1]+H2O_MF/H2O[1]+CO2_MF/CO2[1]+Ar_MF/Ar[1])
    Cv = O2[2]*O2_MF+N2[2]*N2_MF+H2O[2]*H2O_MF+CO2[2]*CO2_MF+Ar[2]*Ar_MF
    Cp = O2[3]*O2_MF+N2[3]*N2_MF+H2O[3]*H2O_MF+CO2[3]*CO2_MF+Ar[3]*Ar_MF
    gamma = Cp/Cv
    P = N2[6]
    Pi =  N2[7]
    return rho, Cv, Cp, gamma, P,Pi,H2O_MF,Ps
 #   return rho,Cp,Cv,gamma,mu,H2O_MF,Ps

def GasProp(Temp,Pres,N2_MF,O2_MF,CO2_MF,CO_MF,H2O_MF,Ar_MF):
    Others_MF = 1- N2_MF - O2_MF - CO2_MF - H2O_MF - Ar_MF
    N2_MF+=Others_MF
    O2 = O2_GasProp(Temp,Pres)
    N2 = N2_GasProp(Temp,Pres)
    H2O = H2O_GasProp(Temp,Pres)
    CO2 = CO2_GasProp(Temp,Pres)
    CO = CO_GasProp(Temp,Pres)
    Ar = Ar_GasProp(Temp,Pres)
    rho = 1/(O2_MF/O2[1]+N2_MF/N2[1]+H2O_MF/H2O[1]+CO2_MF/CO2[1]+Ar_MF/Ar[1]+CO_MF/CO[1])
    Cv = O2[2]*O2_MF+N2[2]*N2_MF+H2O[2]*H2O_MF+CO2[2]*CO2_MF+Ar[2]*Ar_MF+CO[2]*CO_MF
    Cp = O2[3]*O2_MF+N2[3]*N2_MF+H2O[3]*H2O_MF+CO2[3]*CO2_MF+Ar[3]*Ar_MF+CO[3]*CO_MF
    gamma = Cp/Cv
    P = N2[6]
    Pi =  N2[7]
    return rho, Cv, Cp, gamma, P, Pi

humidity = 0
#print(AirProp(288,101.3,0.78,0.21,0.0007,humidity,0.0093))
#print(GasProp(288,500,0.66,0.10,0.12,0,0.11,0.01))


#########WIP for Refprop function################
#from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
#%matplotlib inline
#%config InlineBackend.figure_formats = {'png', 'retina'}

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

