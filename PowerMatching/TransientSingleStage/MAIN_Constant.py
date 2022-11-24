"""
===-*- TurboCompressor Matching Tool -*-===
MAIN Module for Constant flow matching

=====-*- General -*-=====
Copyright (c) makkiblog.com
MIT License 
coding: utf-8

======-*- SUBROUTINE STRUCTURE -*-======
MAIN_Constant.py For Constant Flow
|_StageCalc.py  Classes for Compressor and Turbine Stages

===-*- SYNTAX USAGE -*-===
>python MAIN.py setupfilename matchingmode

===-*- VERSION -*-===
v0.1 Initial version
v0.2 Improved convergence
v1.0 Initial Release Working version
v1.1 Updated convergence parameter
v1.2 Fixed T2t error and updated turbine convergence

vvvCODEvvv
"""

import sys
import argparse
import configparser
import numpy as np
import matplotlib.pyplot as plt
import pandas
import toml
from StageCalc import Stage

#init parser
parser=argparse.ArgumentParser()
parser.add_argument('matching mode', nargs='*', default=[2],help='Enter setup filename and matching mode (Bypass or AFR)')
args=parser.parse_args()

#import setting file
if __name__ == '__main__':
    args = sys.argv
    setupfilename = args[1]
    matchingmode = args[2]
    
    setting_file = './Inputs/'+setupfilename+'.ini'
    setting = configparser.ConfigParser()
    setting.optionxform = str  #Case sensitve option
    setting.read(setting_file, encoding='utf8')
    initrpm = setting.getfloat("Settings", "Initial rpm")
    humidity = setting.getfloat("Environment Variables", "humidity [0-1]") #not used yet
#    matchingmode = setting.get("Settings", "Matching Mode (Bypass or AFR)") # matching mode Bypass or AFR correction
    analysisname = setting.get("Settings", "Analysis Name")
    params_file = setting.get("Settings", "Input filename")
    target_file = setting.get("Settings", "Target filename")
    input_mode = setting.get("Settings", "Input mode (csv or toml)") #input file csv or toml format
    data_no = int(setting.get("Settings", "Data numbers"))
    cmap = setting.get("Settings", "Compressor Map Name")
    tmap = setting.get("Settings", "Turbine Map Name")

#Method for one time iteration matching calculation for flow bypass mode
def MatchingBypassMode(Stgno,setupfilename,Nc):
    Stg = Stage(Stgno,setupfilename,params_file,target_file,input_mode)

    #match compressor
    Nt = Nc
    StgEff = 1 - Stg.S_MechLoss
    BackPres = Stg.T_Pout
    Pin_Tot = Stg.P_Amb
    Cont = 1
    while Cont == 1: 
        PwC,Wcstar,Nc,P2c,T2c,PRC,EtaC = Stg.CompStage(Stgno,cmap,Stg.T_Amb,Pin_Tot,Nc,Stg.C_flow,Stg.C_PRC) 
        Cont,Nc,Nt = Stg.PowBal(PRC,Nc, Nt, Stg.S_Step)
    
    #Calculate P1t and T1t
    W_F,Wt,Pout_comb,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,Stg.CB_AFR,Stg.C_BLDR,Stg.T_BLDR)
    T1t_calc = T1tth*Stg.CB_Tempeff
    P1t = Pout_comb
    PwT,W_reqd,W_Tcalc,W_Tstar,PwAct,Nt,P2t,P1ti1,T2t,PRT,BPR,EtaT = Stg.TurbStage(Stgno,tmap,StgEff,PwC,Stg.CB_t,P1t,Nt,Wt,BackPres)
    Wt_act = Wt
    Enth_act = PwT
    AFR_act = Stg.CB_AFR

    #Iterate P1t for given backpressure, obtain initial pressure
    
    Pgap = abs(P1t - P1ti1)
    while abs(Pgap) > Stg.S_PRCRes:
        Enth_act,W_reqd,W_Tcalc,W_Tstar,PwAct,Nt,P2t,P1ti2,T2t,PRT,BPR,EtaT = Stg.TurbStage(Stgno,tmap,StgEff,PwC,Stg.CB_t,P1ti1,Nt,Wt_act,BackPres)
        Pgap = P1ti2 - P1ti1
        if Pgap > 0:
            P1ti1 += 1/2*Pgap
            #P1ti1 += 1/2*Pgap*Stg.Wtres #Option convergence algorithm, Use for stabilization
        elif Pgap < 0:
            P1ti1 += 1/2*Pgap*Stg.S_PRCRes*BackPres
        #print("P1t Re-Calculating, Pgap = ",Pgap,"BPR=:",BPR,"Enth_act:",Enth_act,"PwAct:",PwAct)
    
    print("Calc Init",Stgno,"BPR=:",BPR,"Enth_act:",Enth_act,"PwAct:",PwAct)
    while Enth_act - PwAct > Stg.S_PowRes:
        Enth_act,W_reqd,W_Tcalc,W_Tstar,PwAct,Nt,P2t,P1ti2,T2t,PRT,BPR,EtaT = Stg.TurbStage(Stgno,tmap,StgEff,PwC,Stg.CB_t,P1ti1,Nt,Wt_act,BackPres)
        Pgap = P1ti2 - P1ti1
        P1ti1 += 1/2*Pgap*Stg.S_PRCRes*BackPres
        Wt_act -= 1/2*(W_Tcalc-W_reqd)*Stg.Wtres
        #Wt_act -= 1/2*(Wt-W_reqd)*BPR*Stg.Wtres #Option convergence algorithm, Use for stabilization
        #print("Bypass:Stage",Stgno,"Pgap:",Pgap,"Wt Re-calc: Actual", Wt_act,"Reqd:",W_reqd,'Enth:',Enth_act,'PwAct:',PwAct,'BPR:',BPR,'PRT',PRT)
    print("Converged Bypass:Stage",Stgno,"Pgap:",Pgap,"Wt Re-calc: Actual", Wt_act,"Reqd:",W_reqd,'Enth:',Enth_act,'PwAct:',PwAct,'BPR:',BPR,'PRT',PRT)

    P1t = P1ti2
    BPR = (Wt - Wt_act)/Wt
    TurbEnth = Enth_act/EtaT/StgEff
    PwT = TurbEnth/(1-BPR)
    T_eff = Stg.CB_t/T1tth

    print('===========<Matched power[kW]>=========\n >>> TOTAL:%.1f, TURBINE:%.1f:\n \n' % (Enth_act, PwAct))
    dict_c = {"Power[kW]":PwC,"Corrected Flow rate Wc*":Wcstar, "Speed RPM":Nc,"Outlet Pres[kPaA]":P2c,"Outlet Temp[degC]":T2c-273.15,"Pres Ratio PRC":PRC,"Comp Efficiency":EtaC,"Air Flow[kg/s]":Stg.C_flow,"CombustorPressure Ratio P2c/P1t":P2c/P1t}
    dict_t = {"Total Enthalpy[kW]":PwT,"Turbine Enthalpy[kW]":TurbEnth,"Turbine Power[kW]":PwAct, "Speed RPM":Nt,"Inlet Pres.[kPa]":P1t, "Inlet Temp[C]":Stg.CB_t-273.15,"Outlet Pres[kPa]Tot":P2t, "Outlet Temp[degC]":T2t-273.15, "Pres Ratio PRT":PRT, "Turb Efficiency":EtaT,"Total Flow[kg/s]":Wt,"Turbine Flow[kg/s]":Wt_act,"Bypass Ratio[BPR]":BPR}
    dict_s = {"Fuel Flow [kg/s]":W_F,"Turbine Required Flow[kg/s]":W_reqd,"Turbine Corrected Flow rate[kg/s]":W_Tstar,"Theoretical T1t[degC]":T1tth-273.15,"Mechanical Loss [kW]":(1-StgEff)*PwAct,"AFR": AFR_act,'Temp Efficiency':T_eff}
    #list_other = [WT_Map,Cont]
    return dict_c,dict_t,dict_s

#Method for one time iteration matching calculation for T1t correction mode
def MatchingAFRMode(Stgno,setupfilename,Nc):
    Stg = Stage(Stgno,setupfilename,params_file,target_file,input_mode)

    #match compressor
    Nt = Nc
    StgEff = 1 - Stg.S_MechLoss
    BackPres = Stg.T_Pout
    Pin_Tot = Stg.P_Amb
    Cont = 1
    while Cont == 1: 
        PwC,Wcstar,Nc,P2c,T2c,PRC,EtaC = Stg.CompStage(Stgno,cmap,Stg.T_Amb,Pin_Tot,Nc,Stg.C_flow,Stg.C_PRC) 
        Cont,Nc,Nt = Stg.PowBal(PRC,Nc, Nt, Stg.S_Step)

    #Calculate P1t and T1t
    W_F,Wt,Pout_comb,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,Stg.CB_AFR,Stg.C_BLDR,Stg.T_BLDR)
    T1t_calc = T1tth*Stg.CB_Tempeff
    P1t = Pout_comb
    PwT,W_reqd,W_Tcalc,W_Tstar,PwAct,Nt,P2t,P1ti1,T2t,PRT,BPR,EtaT = Stg.TurbStage(Stgno,tmap,StgEff,PwC,T1t_calc,P1t,Nt,Wt,BackPres)
    Wt_act = Wt 
    Enth_act = PwT
    AFR_act = Stg.CB_AFR

    #Iterate P1t for given backpressure, obtain initial pressure
    Pgap = abs(P1t - P1ti1)
    Enth_act,W_reqd,W_Tcalc,W_Tstar,PwAct,Nt,P2t,P1ti2,T2t,PRT,BPR,EtaT = Stg.TurbStage(Stgno,tmap,StgEff,PwC,Stg.CB_t,P1ti1,Nt,Wt_act,BackPres)
    #Settle P
    while abs(Pgap) > Stg.S_PRCRes:
        Enth_act,W_reqd,W_Tcalc,W_Tstar,PwAct,Nt,P2t,P1ti2,T2t,PRT,BPR,EtaT = Stg.TurbStage(Stgno,tmap,StgEff,PwC,Stg.CB_t,P1ti1,Nt,Wt_act,BackPres)
        Pgap = P1ti2 - P1ti1
        if Pgap > 0:
            P1ti1 += 1/2*Pgap
            #P1ti1 += 1/2*Pgap*Stg.Wtres #Option convergence algorithm, Use for stabilization
        elif Pgap < 0:
            P1ti1 += 1/2*Pgap*Stg.S_PRCRes*BackPres
        #print("P1t Re-Calculating, Pgap = ",Pgap,"BPR=:",BPR,"Enth_act:",Enth_act,"PwAct:",PwAct)
    print("Converged P1t:Stage",Stgno,"Pgap = ",Pgap,"T1t = ",T1tth,T1t_calc,"BPR=:",BPR,"Enth_act:",Enth_act,"PwAct:",PwAct,'PRT',PRT)

    #Settle T1t
#    Wt_act = W_Tcalc
    #if Enth_act < PwAct*0.1:
    #    Enth_act = 0.5*PwAct #Recovery mode for divergence
      
    while abs(Enth_act - PwAct) > Stg.S_PowRes:
        Pgap = P1ti2 - P1ti1
        P1ti1 += 10*Stg.Wtres*Pgap*abs(BPR)
        #P1ti1 += 1/2*Pgap*Stg.S_PRCRes*BackPres
        if Enth_act > PwAct:
            Stepac = 1/2*(abs(Enth_act - PwAct))**0.1/Enth_act**0.1*AFR_act*Stg.Wtres*abs(BPR)
            AFR_act += Stepac
        elif Enth_act < PwAct:
            Stepac = 1/2*(abs(Enth_act - PwAct))**0.1/Enth_act**0.1*AFR_act*Stg.Wtres
            AFR_act -= Stepac
        W_F,Wt,P1_cb,T1tth = Stg.Combustor(Stg.C_flow,T2c,P2c,AFR_act,Stg.C_BLDR,Stg.T_BLDR)
        dP_comb = P2c - P1ti1
        T1t_calc = T1tth*Stg.CB_Tempeff
        print('AFR:Tot %.2f, Act Enth %.2f, P1ti %.1f, %.1f, T1t %.1f, W_Tcalc %.3f, Wt %.3f, PRT %.3f, BPR %.3f, AFR %.2f' % (Enth_act,PwAct,P1ti1,P1ti2,T1t_calc,W_Tcalc,Wt_act,PRT,BPR,AFR_act))
        Enth_act,W_reqd,W_Tcalc,W_Tstar,PwAct,Nt,P2t,P1ti2,T2t,PRT,BPR,EtaT = Stg.TurbStage(Stgno,tmap,StgEff,PwC,T1t_calc,P1ti1,Nt,Wt_act,BackPres)    
        #print("T1t Iterate, T1t = ",T1tth,T1t_calc,'AFR:',AFR_act,"Delta AFR",Stepac,"Pgap",Pgap,"BPR=:",BPR,"Enth_act:",Enth_act,"PwAct:",PwAct,'EtaT',EtaT)
    print("Converged T1t:Stage",Stgno, "T1t = ",T1tth,T1t_calc,'AFR',AFR_act,"Pgap",Pgap,"BPR=:",BPR,"Enth_act:",Enth_act,"PwAct:",PwAct)
    
    P1t = P1ti2
    BPR = (Wt - Wt_act)/Wt
    TurbEnth = Enth_act/EtaT/StgEff
    PwT = TurbEnth/(1-BPR)
    T_eff = Stg.CB_t/T1tth

    print('=======<Matched power[kW]>=====\n >>> TOTAL:%.1f, TURBINE:%.1f:\n \n' % (Enth_act, PwAct))  
    dict_c = {"Power[kW]":PwC,"Corrected Flow rate Wc*":Wcstar, "Speed RPM":Nc,"Outlet Pres[kPaA]":P2c,"Outlet Temp[degC]":T2c-273.15,"Pres Ratio PRC":PRC,"Comp Efficiency":EtaC,"Air Flow[kg/s]":Stg.C_flow,"Turbine Pressure Ratio P2c/P1t":P2c/P1t}
    dict_t = {"Total Enthalpy[kW]":PwT,"Turbine Enthalpy[kW]":TurbEnth,"Turbine Power[kW]":PwAct, "Speed RPM":Nt,"Inlet Pres.[kPa]":P1t, "Inlet Temp[C]":T1t_calc-273.15,"Outlet Pres[kPa]Tot":P2t, "Outlet Temp[degC]":T2t-273.15, "Pres Ratio PRT":PRT, "Turb Efficiency":EtaT,"Total Flow[kg/s]":Wt,"Turbine Flow[kg/s]":Wt_act,"Bypass Ratio[BPR]":BPR}
    dict_s = {"Fuel Flow [kg/s]":W_F,"Turbine Required Flow[kg/s]":W_reqd,"Turbine Corrected Flow rate[kg/s]":W_Tstar,"Theoretical T1t[degC]":T1tth-273.15,"Mechanical Loss [kW]":(1-StgEff)*PwAct,"AFR": AFR_act,'dP_comb[kPa]':dP_comb,'Temp Efficiency':T_eff}
    #list_other = [WT_Map,Cont]
    return dict_c,dict_t,dict_s


def matching(analysisname,setupfilename,matchingmode,stage,initrpm):
    #Run Matching
    if matchingmode == "Bypass":
        dict_c,dict_t,dict_s = MatchingBypassMode(stage,setupfilename,initrpm)
    elif matchingmode == "AFR":
        dict_c,dict_t,dict_s = MatchingAFRMode(stage,setupfilename,initrpm)
    else:
        print("Matching mode must be selected (Bypass or AFR)")
        exit()

    #Export toml results
    for k,v in dict_c.items():
        dict_c[k] = float("{:.4f}".format(float(v)))
    for k,v in dict_t.items():
        dict_t[k] = float("{:.4f}".format(float(v)))
    for k,v in dict_s.items():
        dict_s[k] = float("{:.4f}".format(float(v)))

    print("==========<COMPRESSOR>==========\n",dict_c,"\n\n=============<TURBINE>===========\n",dict_t,"\n\n======<SHAFT AND COMBUSTOR>======\n",dict_s,"\n") 
    
    res_dict = {"Compressor":dict_c,
    "Turbine":dict_t,
    "Shaft":dict_s}

    #toml_string = toml.dumps(res_dict)  #toml output option
    #output_file = "./Outputs/Result_"+matchingmode+"_"+analysisname+"_"+ str(stage) +".toml"
    #with open(output_file,"w") as output:
    #    toml.dump(res_dict,output)

    #Create multiindex df for csv export
    reformed_dict = {}
    for outerKey, innerDict in res_dict.items():
        for innerKey, values in innerDict.items():
            reformed_dict[(outerKey,
                           innerKey)] = values
    res_df = pandas.Series(reformed_dict,index=reformed_dict.keys())
    res_df = res_df.rename("Stage:"+ str(stage))
    print(res_df)
    #res_df.to_csv('./Outputs/Result_'+matchingmode+'_'+analysisname+"_"+ str(stage) +'.csv')
    return res_df

#Run Matching
res_matrix = pandas.DataFrame()
i=1
while i < data_no + 1:
    res_matrix[("Stage:"+str(i))] = (matching(analysisname,setupfilename,matchingmode,i,initrpm))
    i+=1
print(res_matrix)
res_matrix.to_csv('./Outputs/Result_Summary_'+matchingmode+'_'+analysisname+'.csv')
resplt = res_matrix.transpose()
print('=======RESULT MATRIX========')
print(resplt)
print('============================')

#===========PLOT RESULTS================
fig = plt.figure(figsize=(16,16))

ax1 = fig.add_subplot(1,1,1)
ax1.set_title("Flow-Bypass Ratio")
ax1.set_xlabel("Turbine Total Flow [kg/s]")
ax1.set_ylabel("Bypass Ratio")
ax1.plot(resplt['Turbine']['Total Flow[kg/s]'],resplt['Turbine']['Bypass Ratio[BPR]'])
plt.tight_layout()
plt.savefig("./Outputs/MatchingResult" + analysisname + "_.png", dpi=300)
plt.show()




