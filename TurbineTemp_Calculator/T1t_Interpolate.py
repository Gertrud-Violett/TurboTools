"""
===-*- Turbine Inlet Temperature Correlation Tool -*-===
MAIN Module for Turbine Inlet Temperature Interpolation

=====-*- General -*-=====
Copyright (c) makkiblog.com
MIT License 
coding: utf-8

======-*- SUBROUTINE STRUCTURE -*-======
T1t_Interpolate.py
|_ShowMap.py
|_curvefit.py

===-*- SYNTAX USAGE -*-===
>python T1t_Interpolate.py mapfile targetfile

===-*- VERSION -*-===
v0.1 Initial version
v1.0 Released
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
import plotly
from plotly.subplots import make_subplots
import csv
import pickle
import ShowMap
import argparse
from curvefit import fitting_curve

parser=argparse.ArgumentParser()
parser.add_argument('input file', nargs='*', default=[1], help='input mapfile and target filename (without .csv)')
args=parser.parse_args()

#Read syntax and files
if __name__ == '__main__':
    args = sys.argv
    mapfile = args[1]
    target = args[2]

targetfile =  ("./Input/" + target + ".csv")
mapfile = ("./Input/" + mapfile + ".csv")
df = pd.read_csv(mapfile, sep =",")
tar_df = pd.read_csv(targetfile, sep =",")
skip = 2 #AFR line skip option for plotting 

#obtain keys from dataframe
x = df['rpm']
y = df['AFR']
z = df['T1t']
#alp = df['fuel g/rev']
bet = df['P2c']

#===================INTERPOLATION DEFINITION========================
#Target Cubic Interpolation with P2c (beta) mode: best correlation based on combustion equivalence
T1t_corr = interpolate.interp2d(x = bet, y = y, z = z, kind='cubic', fill_value=None)
def T1t(P2c,AFR,df):
    T1t_calc = T1t_corr(P2c,AFR)
    return(float(T1t_calc))

def Target_T1t(calc_df):
    index = calc_df.index
    P2c = calc_df['P2c']
    AFR = calc_df['AFR']
    T1t_res = []
    for i in index:
        T1t_res.append(T1t(P2c[i],AFR[i],df))
    calc_df.insert(6,"T1t interp",T1t_res)
    return(calc_df)

#4D interpolation XArray method and csv export
def Target_T1t_4D(calc_df2):
    rpm = calc_df2['rpm']
    AFR = calc_df2['AFR']
    P2c = calc_df2['P2c']
    df_rows = pd.DataFrame(df).set_index(["rpm", "AFR", "P2c"])
    df_rows = df_rows[~df_rows.index.duplicated()]
    ds = xr.Dataset.from_dataframe(df_rows)
    print(ds)
    dataarr = ds.interp(rpm=rpm,AFR=AFR,P2c=P2c)
    return(dataarr)

#=====================Scipy Curvefit================================
#Curve fit Method For Altitude MODE
def Target_T1t_curvefit(df,target_df):
    Cs,Error,Stats = fitting_curve(df)
    def T1t_fit(rpm,AFR,P1t,P1c,Cs):
        T1t = Cs[0]+Cs[1]*rpm**2+Cs[2]*rpm+Cs[3]*AFR**2+Cs[4]*AFR+Cs[5]*P1t**2+Cs[6]*P1t+Cs[7]*P1c
        ERRmax = np.sqrt(Error[0]**2+Error[1]**2*rpm**2+Error[2]**2*rpm+Error[3]**2*AFR**2+Error[4]**2*AFR+Error[5]**2*P1t**2+Error[6]**2*P1t+Error[7]*P1c)
        T1tmax = T1t + ERRmax
        T1tmin = T1t - ERRmax
        return T1t,T1tmax,T1tmin
    index = target_df.index
    rpm = target_df['rpm']
    AFR = target_df['AFR']
    P1t = target_df['P1t']
    P1c = target_df['P1c']
    T1t_curvefit = []
    T1t_fitmax = []
    T1t_fitmin = []
    for i in index:
        T1t,T1tmax,T1tmin = T1t_fit(rpm[i],AFR[i],P1t[i],P1c[i],Cs)
        T1t_curvefit.append(T1t)
        T1t_fitmax.append(T1tmax)
        T1t_fitmin.append(T1tmin)
    target_df.insert(8,"T1t curvefit",T1t_curvefit)
    target_df.insert(9,"T1t fitmax",T1t_fitmax)
    target_df.insert(10,"T1t fitmin",T1t_fitmin)
    return(target_df,Stats)

#================TARGET CALCULATION=================
#Correlation result for target
#2D AFR and P2c Method
result_df = Target_T1t(tar_df)
print("======Target Results=======")
print('result2D',result_df)
result_df.to_csv('./Output/T1t_result2D_' + target + "_.csv")

#4D Linear Method
result_ds = Target_T1t_4D(tar_df)
result_4dpd = result_ds.to_dataframe()
print('result4D',result_4dpd)
#result_4dpd.to_csv('./Output/T1t_result4D_' + target + "_.csv")

#scipy curvefit method
result_curvefit,Stats = Target_T1t_curvefit(df,tar_df)
print('curevefit',result_curvefit)
result_curvefit.to_csv('./Output/T1t_result_curvefit_' + target + "_.csv")


#==================MESHGRID PLOTTING===================
#Plot map data to meshgrid
def correlatemap2d(xin,yin,zin,name):
    xx = np.linspace(np.min(xin),np.max(xin))
    yy = np.linspace(np.min(yin),np.max(yin))
    xx1,yy1 = np.meshgrid(xx,yy)
    T1t =  interpolate.griddata((xin,yin),zin,(xx1.ravel(),yy1.ravel()))
    T1t_arr =  interpolate.griddata((xin,yin),zin,(xx1,yy1))
    yy1_list = np.append([0],yy1[:,0])
    T1t_arr = np.append(np.array([xx1[0]]),T1t_arr,axis=0)
    T1t_arr = np.append(np.array([yy1_list]).transpose(),T1t_arr,axis=1)
    np.savetxt("./Output/" + name + "_.csv", T1t_arr, delimiter=",")
    return(xx1,yy1,T1t,T1t_arr)

#Plot 2D map data to gridmap
xx1,yy1,T1t_1,T1t_1_arr = correlatemap2d(x,y,z,"T1t-1_rpm-AFR_" + target)
#alpha,yy2,T1t_2,T1t_2_arr = correlatemap2d(alp,y,z,"T1t-2_fuel-AFR_" + target)
#beta,yy3,T1t_3,T1t_3_arr = correlatemap2d(bet,y,z,"T1t-3_P2c-AFR_" + target)
xx4,yy4,T1t_4,T1t_4_arr = correlatemap2d(x,y,z,"T1t-4_rpm_curvefit")


#plotting 2D Method
def plot2d(dataset,title,xaxis,xliml,xlimh,skip):
    fig2d = plt.figure(figsize=(16,16))
    ax1 = fig2d.add_subplot(1,1,1)
    ax1.set_title(title)
    ax1.set_xlabel(xaxis)
    ax1.set_xlim(xliml,xlimh)
    ax1.set_ylim(np.min(z),np.max(z))
    ax1.set_ylabel("T1t")
    ncol = round(len(dataset[:,0])/skip)
    colors = plt.cm.plasma(np.linspace(0,1,ncol))
    i = 1
    while i < ncol:
        ax1.plot(dataset[0,:],dataset[i,:],label=( "%.2f" % dataset[i,0]),color=colors[i])
        i += 1
    ax1.legend()
    plt.savefig("./Output/T1t_2d_" + target + title + "_.png", dpi=300)
    

#==============PLOT RESULTS FOR INTERP====================
#plotting 2D
plot2d(T1t_1_arr,"T1t-1_rpm-AFR_","[rpm]",np.min(x),np.max(x),skip)
#plot2d(T1t_2_arr,"T1t-2_fuel-AFR_","fuel rate[g/rev]",np.min(alp),np.max(alp),skip)
#plot2d(T1t_3_arr,"T1t-3_P2c-AFR_","P2c [kPa]",np.min(bet),np.max(bet),skip)

#plotting 3D
fig = plt.figure(figsize=(16,16))

ax1 = fig.add_subplot(1,1,1, projection='3d')
ax1.set_title("T1t_rpm-AFR")
ax1.set_xlabel("rpm")
ax1.set_ylabel("AFR")
ax1.set_zlabel("T1t[degC]")
ax1.set_zlim(min(df['T1t']),max(df['T1t']))
ax1.view_init(elev=10., azim=100)
ax1.plot_trisurf(xx1.ravel(), yy1.ravel(), T1t_1, cmap=plt.cm.Spectral_r)
rpminit = result_df.at[1,'rpm']
if math.isnan(rpminit) == False:
    ax1.scatter(result_df['rpm'],result_df['AFR'],result_df['T1t interp'],color='Black')
else:
    pass

"""
#Optional P2c and fuelflow interpolation
ax2 = fig.add_subplot(2,2,2, projection='3d')
ax2.set_title("T1t_fuelflow-AFR")
ax2.set_xlabel("fuel rate [g/rev]")
ax2.set_ylabel("AFR")
ax2.set_zlabel("T1t[degC]")
ax2.set_zlim(400,750)
ax2.view_init(elev=10., azim=140)
ax2.plot_trisurf(alpha.ravel(), yy2.ravel(), T1t_2, cmap=plt.cm.Spectral_r)
fuelinit = result_df.at[1,'fuel rate']
if math.isnan(fuelinit) == False:
    ax2.scatter(result_df['fuel rate'],result_df['AFR'],result_df['T1t interp'],color='Black')
else:
    pass

ax3 = fig.add_subplot(2,2,3, projection='3d')
#ax3 = fig.add_subplot(1,1,1, projection='3d')
ax3.set_title("T1t_P2c-AFR")
ax3.set_xlabel("P2c [kPa]")
ax3.set_ylabel("AFR")
ax3.set_zlabel("T1t[degC]")
ax3.set_zlim(400,750)
ax3.view_init(elev=10., azim=100)
ax3.plot_trisurf(beta.ravel(), yy3.ravel(), T1t_3, cmap=plt.cm.Spectral_r)
P2cinit = result_df.at[1,'P2c']
if math.isnan(P2cinit) == False:
    ax3.scatter(result_df['P2c'],result_df['AFR'],result_df['T1t interp'],color='Black')
else:
    pass
"""

plt.tight_layout()
#pickle.dump(fig, open('./Output/T1t_map' + target + '_pickle', 'wb'))
plt.savefig("./Output/T1t_map3d_" + target + "_.png", dpi=300)
plt.clf()

#============Curvefit 3D Plot================
fig2 = plt.figure(figsize=(16,16))

ax1 = fig2.add_subplot(1,1,1, projection='3d')
ax1.set_title("T1t_rpm-AFR_Curvefit")

ax1.set_xlabel("rpm")
ax1.set_ylabel("AFR")
ax1.set_zlabel("T1t[degC]")
#ax1.set_zlim(400,750)
ax1.set_zlim(min(df['T1t']),max(df['T1t']))
ax1.view_init(elev=10., azim=100)
ax1.plot_trisurf(xx4.ravel(), yy4.ravel(), T1t_4, cmap=plt.cm.Spectral_r)
rpminit = result_curvefit.at[1,'rpm']
if math.isnan(rpminit) == False:
    ax1.scatter(result_curvefit['rpm'],result_curvefit['AFR'],result_curvefit['T1t fitmax'],label='max.',color='Red')
    ax1.scatter(result_curvefit['rpm'],result_curvefit['AFR'],result_curvefit['T1t curvefit'],label='nominal',color='Black')
    ax1.scatter(result_curvefit['rpm'],result_curvefit['AFR'],result_curvefit['T1t fitmin'],label='min.',color='Blue')
    ax1.text(2000, 16, max(df['T1t'])+50, (Stats))
else:
    pass
ax1.legend()
plt.tight_layout()
#pickle.dump(fig, open('./Output/T1t_map' + target + '_pickle', 'wb'))
plt.savefig("./Output/T1t_map3d_curvefit" + target + "_.png", dpi=300)
plt.show()

print(tar_df)
