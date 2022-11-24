"""
===-*- Turbine Inlet Temperature Correlation Tool -*-===
curvefit Module for Turbine Inlet Temperature Interpolation

=====-*- General -*-=====
Copyright (c) makkiblog.com
MIT License 
coding: utf-8

===-*- VERSION -*-===
v0.1 Initial version
v1.0 Released

Reference: https://rikei-fufu.com/2020/07/05/post-3270-fitting/
vvvCODEvvv
"""

#===========================INIT=====================================
import sys
import math
import configparser
import numpy as np
import matplotlib.pyplot as plt, mpld3
import importlib
import pandas as pd
import xarray as xr
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy import stats
#from scipy.interpolate import LinearNDInterpolator as lNDI
from mpl_toolkits.mplot3d import Axes3D


#=====================Scipy Curvefit================================
#Curve fit Method
def fitting_curve(df):
    def t1tcurve(X,C0,C1,C2,D1,D2,E1,E2,F1): # 2D Gaussian
        x, y, bet, gam = X #rpm, AFR,P1t
        z = C0+C1*x**2+C2*x+D1*y**2+D2*y+E1*bet**2+E2*bet+F1*gam
        return z
 
    x_obs = df['rpm']
    y_obs = df['AFR']
    bet_obs = df['P1t']
    gam_obs = df['P1c']
    z_obs = df['T1t']
 
    #fittingのメイン計算部分
    popt, pcov = curve_fit(t1tcurve, (x_obs, y_obs, bet_obs, gam_obs), z_obs) 
    perr = np.sqrt(np.diag(pcov)) 
 
    
    #Chi2 contingency
    o = z_obs
    e = t1tcurve((x_obs, y_obs, bet_obs, gam_obs), popt[0], popt[1], popt[2], popt[3], popt[4],popt[5],popt[6],popt[7]) 
    chi2 = stats.chisquare(o, f_exp = e, ddof=4) #カイ自乗計算のメイン部分。chi2には[カイ二乗, p値]の2つが出力される。
 
    #R2 calc
    residuals =  o - e #残渣
    rss = np.sum(residuals**2)      #残差平方和: residual sum of squares = rss
    tss = np.sum((o-np.mean(o))**2) #全平方和: total sum of squares = tss
    r_squared = 1 - (rss / tss)     #決定係数R^2
    statistics_numbers = {  "X-squared": format(chi2[0], '.3f'),
                            "p-value": format(chi2[1], '.5f'),
                            "R^2": format(r_squared, '.4f')}

    #Print results
    print("==========Fit Result============")
    print("T1t = C0 + C1*rpm^2 + C2*rpm + D1*AFR^2 + D2*AFR + E1*P1t^2 + E2*P1t + F1*P1c")
    print("Constants:",popt)
    print('Error:',perr)   
    print("X-squared, p-value, R^2:", statistics_numbers)
    print("=============EOF=============") 
    return popt,perr,statistics_numbers

if __name__=="__main__":
    Constants,Error,Stats = fitting_curve(df)
