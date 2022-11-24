"""
-*- coding: utf-8 -*-
Copyright (c) makkiblog.com
MIT License 
-*- Low Cycle Fatigue life calculation tool-*-
v0.1 Released

-*- SYNTAX USAGE -*-
>python LCF_Calc.py mapdata.csv dutycycle.csv MSRS
(MSRS: Maximum Sea Level Rated Speed)
vvvCODEvvv
"""

#===========================INIT=====================================
import sys
import math
import argparse
import configparser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
import importlib
import pandas as pd
from scipy import interpolate
import csv
import seaborn as sns
import progressbar
from time import sleep

#============PROGRAM DEFAULT SETTING (Do not Edit)=========
msrsratiomax = 1.0  #Max msrs range for sensitivity study
msrsratiomin = 0.9  #Min msrs range for sensitivity study
#LPF = 5000 #Low Pass Filter Setting = 1000rpm. Dismiss changes within this range
DefFailRate = 0.005 #Default Failure Rate B0.5 (0.5%)
LifeplotResolution = 3  #No of MSRS points to be plotted in life plot: Def.=5


#=====================Init====================
parser=argparse.ArgumentParser()
#parser.add_argument('input filename', nargs='*', default=[], help='||Syntax: Need to input:|| \n Damage map filename (without .csv) \n Duty cycle filename (without .csv) \n Max. Sea level rated Speed [rpm] \n multiplier for input and output data (seconds/hour, cycles/km etc.) \n Low pass filter settings')
parser.add_argument('Damage Map', help='Damage map filename (without .csv)')
parser.add_argument('Duty Cycle', help='Damage map filename (without .csv)')
parser.add_argument('MSRS [rpm]', help='Maximum Sea Level Rated Speed')
parser.add_argument('Multiplyer', help='multiplier for input and output data (seconds/hour = 3600, cycles/km etc.)')
parser.add_argument('LPF Value[rpm]', help='Low pass filter settings. Enter Value where min/max rpms delta is dismissed')
args=parser.parse_args()

#Read syntax and files
if __name__ == '__main__':
    args = sys.argv
    damagemap = args[1]  #No data type, string, damage map filename
    dutycycle = args[2]  #No data type, string, dutycycle filename
    msrs = args[3]   #No data type, float, Enter Max. sea level rated speed
    conv = args[4]  #No data type, float, Enter mutiplyer for input and output data: e.g. seconds to hours: 3600
    LPF = float(args[5])    #Low Pass Filter Setting = 1000rpm. Dismiss changes within this range

dutycyclefile =  ("./Input/" + dutycycle + ".csv")
mapfile = ("./Input/" + damagemap + ".csv")
df = pd.read_csv(mapfile, sep =",",  index_col=0)
tar_df = pd.read_csv(dutycyclefile, sep =",")


#================DATA CONVERSION==============
#Convert to np array
dfarr = df.to_numpy()
df1d = df.stack()
dfZarr = df1d.to_numpy()

#Create 1d list for X and Y coordinates for griddata interpolation 脳筋!スマートメソッド求む
dfX = list(np.float_(df.columns))
dfY = list(np.float_(df.index))
dfXarr = []
dfYarr = []
for row in dfY:
    dfXarr = np.append(dfXarr,dfX,axis=0)
n=0
for row in dfY:
    for column in dfX:
        dfYarr.append(dfY[n])
    n+=1


#=================DEF PLOTS===================
#Plot duty cycle
def plotdutycycle(dataset):
    fig2d = plt.figure(figsize=(10,10))
    ax = fig2d.add_subplot(1,1,1)
    ax.set_title("DutyCycle")
    ax.set_xlabel("Time[s]")
    ax.set_ylabel("rpm")
    print(dataset)
    ncol = len(dataset.iloc[0,:])
    colors = plt.cm.hsv(np.linspace(0,1,ncol))
    i = 1
    while i < ncol:
        ax.plot(dataset.iloc[:,0],dataset.iloc[:,i],label=(dataset.columns[i]),color=colors[i])
        i += 1
    ax.legend()
    plt.savefig('./Output/dutycycle_' + dutycycle + '_MSRS_' + str(msrs) + '_LPFrpm_' + str(LPF) + '_.png', dpi=300)
    #plt.show()

#Plot nparray heatmap
def plotheatmap(dataset):
    X = list(np.float_(dataset.columns))
    Y = list(np.float_(dataset.index))
    print(X)
    print(Y)
    fig = plt.figure(figsize=(10,12))
    ax = fig.add_subplot(1,1,1)
    ax.set_xticks(np.arange(0, dataset.index.max(), step=10))
    ax.set_yticks(np.arange(0, dataset.index.max(), step=10))
    ax.set_xlabel("rpm min")
    ax.set_ylabel("rpm max")
    ax.set_title("damage map")
    im = ax.pcolormesh(X,Y,dataset, norm = colors.LogNorm(), cmap= 'RdBu',)
    #im = ax.contourf(X,Y,dataset, locator=ticker.LogLocator(), cmap= 'RdBu',)
    fig.colorbar(im)
    plt.savefig('./Output/damagemap_'  + damagemap + '_.png', dpi=300)
    #plt.show()


#===============CALCULATE DAMAGE================
def damageget(X,Y,Z,x,y):
    points = X,Y
    values = Z
    z =  interpolate.griddata(points,values,(x,y), method = 'linear')
    return z


def calclife(map,data,j,duration,rpm):
    print("No.",j,"Max rpm",max(data.iloc[:,j]),"ratio",ratio)
    print("duration for duty cycle (check units):",duration)
    k = 1
    i = k
    #Life counter
    delta = data.iloc[k-1,j] - data.iloc[i,j]
    damagelist = []
    rpmmin = data.iloc[k,j]
    rpmmax = data.iloc[k,j]
    print("Start Calculation for item no. %i at %.0f rpm" % (j, float(rpm)))
    print("data length", len(data.iloc[:,0]))
    bar = progressbar.ProgressBar(maxval=len(data.iloc[:,0]-1), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    while k < len(data.iloc[:,0]-1):
        while delta >= 0:
            if i < len(data.iloc[:,0]-1):
                delta = data.iloc[i-1,j] - data.iloc[i,j]
                rpmmin = data.iloc[i-1,j]
                i+=1
                #print("P",i,delta)
            else:
                break
        while delta < 0:
            if i < len(data.iloc[:,0]-1):
                delta = data.iloc[i-1,j] - data.iloc[i,j]
                rpmmax = data.iloc[i-1,j]
                i+=1
                #print("N",i,delta)
            else:
                break
        life = damageget(dfXarr,dfYarr,dfZarr,rpmmin/10**3,rpmmax/10**3)
        if rpmmax - rpmmin < LPF:
            damage = 0
        else:
            damage=1/life
        damagelist.append(damage)
        k = i
        #print("k",k,"i",i, "max", rpmmax, "min",rpmmin,"damage", damage)
        bar.update(k)
    bar.finish()
    Totaldamage=sum(damagelist)
    Totallife = duration/Totaldamage
    print("list length",len(damagelist))
    print("total damage",Totaldamage, "life", Totallife)
    return Totaldamage, Totallife


def failurerate(defrate,deflife):
    time = []
    time.append(100)
    i = 0
    while time[i] < 10e5:
        time.append(time[i]*1.1)
        i+=1
    j=0
    calclife = []
    for times in time:
        calclife.append(1 - (1 - defrate)**((times / deflife)**3))
        j+=1
    return time, calclife


#=====================RESULTS===========================
print("Damage Map\n", df)
print("\n Duty Cycle \n", tar_df)
plotheatmap(df)
plotdutycycle(tar_df)

ncol = int(len(tar_df.iloc[0,:]))
msrsarr = np.linspace(float(msrs)*msrsratiomin,float(msrs)*msrsratiomax,num=LifeplotResolution)
result_df = pd.DataFrame(columns=tar_df.columns)
result_df.iloc[:,0] = msrsarr
#totaldamage,result_df.iloc[0,1] = calclife(df,tar_df,1,dur,msrs) #Multiplier seconds to hours, first duty cycle only

#sweep through msrs
no = 1
while no < ncol:
    k=0
    for rpms in msrsarr:
        print("MSRS",rpms)
        rpmdf = tar_df
        dur = rpmdf.iloc[:,0].max()/float(conv)
        ratio = float(rpms)/float(msrs)
        rpmdf = rpmdf.multiply(ratio)
        totaldamage,result_df.iloc[k,no] = calclife(df,rpmdf,no,dur,rpms) #Multiplier seconds to hours
        k+=1
        print("")
    no+=1
result_df = result_df.rename(columns={'Time [s]':'rpm'})
result_df.to_csv('./Output/result_' + dutycycle + '_MSRS_' + str(msrs) + '_LPFrpm_' + str(LPF) + '_.csv')

#failure rate calculation
life_df = pd.DataFrame()
res_df_drop = result_df.drop(['rpm'],axis=1)
for column in res_df_drop.columns:
    for index in res_df_drop.index:
        input = res_df_drop.loc[index,column]    
        colname = str("%.0f _" % (result_df.loc[index,'rpm'])) + '_' + column
        life_df['Time'],life_df[colname]  = failurerate(DefFailRate,input)
life_df.to_csv('./Output/life_' + dutycycle + '_MSRS_' + str(msrs) + '_LPFrpm_' + str(LPF) + '_.csv')


#=================PLOT RESULTS===============
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.set_title("Life-rpm (MSRS)")
ax.set_ylim(10e2,10e6)
ax.set_xlabel("rpm")
ax.set_ylabel("Life")
ncol = len(result_df.iloc[0,:])
ax.set_yscale('log')
ax.set_yticks(np.logspace(2, 5, num=16, base=10, endpoint=True))
l = 1
while l < ncol:
    ax.plot(result_df['rpm'],result_df.iloc[:,l],label=(result_df.columns[l]))
    l += 1
ax.legend()
plt.savefig('./Output/MSRS-Life_' + dutycycle + '_MSRS_' + str(msrs) + '_LPFrpm_' + str(LPF) + '_.png', dpi=300)
#plt.show()

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
ax.set_title("Life-FailureRate")
ax.set_xlabel("distance/hours")
ax.set_ylabel("Failure rate")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(10e2,10e5)
ax.set_ylim(10e-5,1)
ax.set_xticks(np.logspace(2, 5, num=16, base=10, endpoint=True))
ax.set_yticks(np.logspace(-5, 0, num=16, base=10, endpoint=True))
life_dfpop = life_df.drop(['Time'], axis=1)
for column in life_dfpop.columns:
    ax.plot(life_df['Time'],life_df[column],label=(column))
ax.legend()
plt.savefig('./Output/B-Life_' + dutycycle + '_MSRS_' + str(msrs) + '_LPFrpm_' + str(LPF) + '_.png', dpi=300)
