"""
EasyPintle
ピントルインジェクタの形状から各基本諸元を吐き出すスクリプト

Primary: 1列目
Secondary: 2列目
Tertiary: 3列目 (補助穴、予備燃焼穴)

各パラメータは下記参照のこと
RP_OR_PRP_REF_用語定義_ピントルインジェクタ形状
https://drive.google.com/open?id=10SRJ_7X-CTOBotop3nZ1QJVhWbI7XsVbEQAcDDQ_6j4

3列目については、穴径、Cd値が1/2列目と著しく違うため、別枠として考え、流量比のみをみる

"""



import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import importlib
from scipy import interpolate
from scipy import optimize
from scipy.optimize import newton
import configparser
from matplotlib.font_manager import FontProperties
from coolant_table import Coolant

plt.close('all')


class Pintle:
    def __init__(self, setting_file, reload = False):
        if (reload):  #再読込の際はself.settingの値をそのまま使う
            pass
        else:
            self.setting_file = setting_file
            self.setting = configparser.ConfigParser()
            self.setting.optionxform = str  # 大文字小文字を区別するおまじない
            self.setting.read(setting_file, encoding='utf8')
        setting = self.setting

        self.name = setting.get("全体", "名前")
        #self.is_save_fig = setting.getboolean("全体", "図を保存する？")
        InjType = setting.getfloat("Type", "Circular[1] or Rectangular[2]")
        Dist_Dp = setting.getfloat("PintleDim", "Pintle Tip Diameter Dp[mm]")/10**3
        Dist_Ls = setting.getfloat("PintleDim", "Skip Distance Ls[mm]")/10**3
        Dist_Dfo = setting.getfloat("PintleDim", "Fuel channel Outer Diameter Dfo[mm]")/10**3
        Dist_Dfm = setting.getfloat("PintleDim", "Fuel manifold Inner Diameter Dfm[mm]")/10**3
        Dist_Doi = setting.getfloat("PintleDim", "LOx channel Inner Diameter Doi [mm]")/10**3
        Dist_Lo1 = setting.getfloat("PintleDim", "Primary Oxidizer slot height Lo1[mm]")/10**3
        Dist_deltao1 = setting.getfloat("PintleDim", "Primary Oxidizer slot Circumferential size deltao1[mm]")/10**3
        Dia_Lo1 =  setting.getfloat("PintleDim", "Primary Oxidizer Diameter [mm]")/10**3
        Dist_Lo2 = setting.getfloat("PintleDim", "Secondary Oxidizer slot height Lo2[mm]")/10**3
        Dist_deltao2 = setting.getfloat("PintleDim", "Secondary Oxidizer slot Circumferential size deltao2[mm]")/10**3
        Dia_Lo2 =  setting.getfloat("PintleDim", "Secondary Oxidizer Diameter [mm]")/10**3
        CN_N1 = setting.getfloat("PintleDim", "Number of Primary slots N1")
        CN_N2 = setting.getfloat("PintleDim", "Number of Secondary slots N2")
        CN_N3 = setting.getfloat("PintleDim", "Number of Tertiary slots N3")
        Dist_Lo12 = setting.getfloat("PintleDim", "Distance between Primary and Secondary Oxidizer slots Lo12[mm]")/10**3
        rho_o = setting.getfloat("OxidizerProp", "LOx Density[kg/m3]")
        mdot_os = setting.getfloat("OxidizerProp", "LOx mass flow rate [kg/s]")
        P_fup = setting.getfloat("FuelProp", "Fuel Injector Upstream Pressure[MPa]")  
        T_fup = setting.getfloat("FuelProp", "Fuel Injector Upstream Temperature[K]")
        mdot_f = setting.getfloat("FuelProp", "Fuel mass flow rate [kg/s]")
        LOx_Cd = setting.getfloat("InjectorParam", "LOx injector Cd")
        Fuel_Cd = setting.getfloat("InjectorParam", "Fuel injector Cd")
        mr_Tert_o = setting.getfloat("InjectorParam", "Tertiary slot flow ratio")
        P_fdo_guess = setting.getfloat("Hypothesis", "Fuel Injector Downstream Pressure[MPa]")

        coolant = Coolant()
        coolant.load_axis('./TableData/Methane_table/axes.csv')
        coolant.load_table('./TableData/Methane_table/density.csv', 'rho')
        coolant.load_table('./TableData/Methane_table/cp.csv', 'cp')
        coolant.load_table('./TableData/Methane_table/cv.csv', 'cv')
        coolant.load_table('./TableData/Methane_table/viscosity.csv', 'viscosity')
        coolant.load_table('./TableData/Methane_table/k.csv', 'thermal_conductivity')
        coolant.load_table('./TableData/Methane_table/prandtl.csv', 'prandtl')
        coolant.load_table('./TableData/Methane_table/gamma.csv', 'gamma')
        coolant.load_table('./TableData/Methane_table/enthalpy.csv', 'enthalpy')
        coolant.load_table('./TableData/Methane_table/temperature.csv', 'temperature', rowAxis='pressure',colAxis='enthalpy')

        rho_f = coolant.rho(T_fup, P_fup*1e6)
        gamma_f = coolant.gamma(T_fup, P_fup*1e6)

        #Injector geometry area
        if InjType == 1:
            self.Area_LOx1 = Dia_Lo1**2*np.pi/4*CN_N1
            self.Area_LOx2 = Dia_Lo2**2*np.pi/4*CN_N2
            self.BLF1 = Dia_Lo1/np.pi/Dist_Dp*CN_N1
            self.BLF2 = Dia_Lo2/np.pi/Dist_Dp*CN_N2
        elif InjType == 2:
            self.Area_LOx1 = Dist_Lo1*Dist_deltao1*CN_N1
            self.Area_LOx2 = Dist_Lo2*Dist_deltao2*CN_N2
            self.BLF1 = Dist_deltao1/np.pi/Dist_Dp*CN_N1
            self.BLF2 = Dist_deltao2/np.pi/Dist_Dp*CN_N2
        else:
            print("error")

        self.Area_LOx = self.Area_LOx1+self.Area_LOx2
        self.Area_Fuel = (Dist_Dfo**2-Dist_Dp**2)*np.pi/4
        self.LOx_PA_ratio = self.Area_LOx1/(self.Area_LOx1+self.Area_LOx2)

        #flow velocity
        mdot_o = mdot_os * (1 - mr_Tert_o)
        self.LOx_vo1 = mdot_o*self.LOx_PA_ratio/rho_o/self.Area_LOx1
        self.LOx_vo2 = mdot_o*(1-self.LOx_PA_ratio)/rho_o/self.Area_LOx2
        self.LOx_vo = self.LOx_vo1
        self.Fuel_vf = mdot_f/rho_f/self.Area_Fuel
        self.Mom_Ff = rho_f*self.Fuel_vf**2*self.Area_Fuel
        self.Mom_Fo1 = rho_o*self.LOx_vo1**2*self.Area_LOx1
        self.Mom_Fo2 = rho_o*self.LOx_vo2**2*self.Area_LOx2
        self.Mom_Fo = self.Mom_Fo1+self.Mom_Fo2
        self.TMR = self.Mom_Ff/self.Mom_Fo
        self.man_lox = mdot_o/rho_o/(Dist_Doi**2*np.pi/4)
        self.man_fuel = mdot_f/rho_f/((Dist_Dfm**2-Dist_Dp**2)*np.pi/4)

        #Other Parameters
        self.Skip_Dist_D = Dist_Ls/Dist_Dp
        self.Skip_Dist_V = Dist_Ls/self.Fuel_vf
        self.ATM_Cone = np.arctan(mdot_o/mdot_f)*180/np.pi
        self.mdot_Tert_o = mr_Tert_o*mdot_os
        self.Tert_Dia = (self.mdot_Tert_o*4/(self.LOx_vo1*rho_o*CN_N3*np.pi))**0.5*10**3
        self.rho_f = rho_f
        self.P_fup = P_fup

        #Pressure loss delta p
        self.deltap_o =  mdot_o**2/(2*rho_o*self.Area_LOx**2*LOx_Cd**2)/10**6        
        
        mdot_f_cal = self.Area_Fuel*np.sqrt(P_fup*1e6*rho_f*gamma_f)*(2/gamma_f)**((gamma_f+1)/(gamma_f-1))
        mdot_f_exp = mdot_f/Fuel_Cd

        #圧縮性(ガス、超臨界) or 非圧縮(液体)で場合分け
        #臨界点(Tm, Pm)の密度を閾値として液体か、その他か判定
        if rho_f  > 161.383 :
          self.deltap_f = mdot_f**2/(2*rho_f*self.Area_Fuel**2*Fuel_Cd**2)/10**6
          self.mdot_f_con = mdot_f
          self.P_fdo = P_fup - self.deltap_f

        else:
        #インジェクタ下流圧力 P_fdo[MPa]を求める。
            def pressure_cal(x):
            
                return self.Area_Fuel*np.sqrt(P_fup*1e6*rho_f)*np.sqrt(2*gamma_f/(gamma_f-1)*((x/P_fup)**(2/gamma_f)-(x/P_fup)**((gamma_f+1)/gamma_f))) - mdot_f_exp

            P_fdo = newton(pressure_cal, P_fdo_guess, disp=False)
            mdot_f_con = Fuel_Cd*self.Area_Fuel*np.sqrt(P_fup*1e6*rho_f)*np.sqrt(2*gamma_f/(gamma_f-1)*((P_fdo/P_fup)**(2/gamma_f)-(P_fdo/P_fup)**((gamma_f+1)/gamma_f)))
            self.P_fdo = P_fdo
            self.mdot_f_con = mdot_f_con
        
            if P_fdo > P_fup/2:
                self.deltap_f = P_fup - P_fdo
            else:
              self.deltap_f = P_fup/2
          

    def display(self):
        print("")
        print("TMR (Total Momentum Ratio) :\t\t%.2f " % (self.TMR))
        print("LOx Primary slot flow ratio:\t\t%.2f " % (self.LOx_PA_ratio))
        print("LOx Outlet Velocity vo1 [m/s] :\t\t%.1f " % (self.LOx_vo1))
        print("Fuel Outlet Velocity vf [m/s] :\t\t%.1f " % (self.Fuel_vf))
        print("")
        print("LOx Manifold Velocity [m/s] :\t\t%.1f " % (self.man_lox))
        print("Fuel Manifold Velocity [m/s] :\t\t%.1f " % (self.man_fuel))
        print("Non-Dimensional skip distance Ls/Dp :\t%.2f " % (self.Skip_Dist_D))
        print("Normalized skip distance Ls/vf[s] :\t%.5f " % (self.Skip_Dist_V))
        print("Theoretical Atomizing Cone angle[deg]:\t%.2f " % (self.ATM_Cone))
        print("")
        print("Primary slot Blockage Factor :\t\t%.2f " % (self.BLF1))
        print("Secondary slot Blockage Factor :\t%.2f " % (self.BLF2))
        print("Fuel Density [kg/m^3]:  \t\t%.2f " % (self.rho_f))
        print("Oxidizer Injector delta p [MPa]:\t%.2f " % (self.deltap_o))
        print("Tertiary slot flow rate [kg/s]:\t\t%.3f " % (self.mdot_Tert_o))
        print("Tertiary slot flow Dia [mm]:\t\t%.3f " % (self.Tert_Dia))
        print("")
        if self.rho_f > 161.383 : 
                print("Fuel is uncompressed")
                print("Fuel Injector delta p [MPa]:\t\t%.2f " % (self.deltap_f))
                print("Fuel Injector Downstream Pressure [MPa]:%.2f " % (self.P_fdo))
        else:
            
            if self.P_fdo > self.P_fup/2:
                print("Fuel is compressed and non-choked")
                print("Fuel Injector delta p [MPa]:\t\t%.2f " % (self.deltap_f))
                print("Fuel Injector Downstream Pressure [MPa]:%.2f " % (self.P_fdo))
            else:
                print("Fuel is choked")
                print("Fuel Injector delta p is less than " + str(self.P_fup/2) +"[MPa]")

    def print(self):
        with open("PintleParams.out","w") as output:
            print("TMR (Total Momentum Ratio) :\t\t%.2f " % (self.TMR),file=output)
            print("LOx Primary slot flow ratio:\t\t%.2f " % (self.LOx_PA_ratio),file=output)
            print("LOx Outlet Velocity vo1 [m/s] :\t\t%.1f " % (self.LOx_vo1),file=output)
            print("Fuel Outlet Velocity vf [m/s] :\t\t%.1f " % (self.Fuel_vf),file=output)
            print("Non-Dimensional skip distance Ls/Dp :\t%.2f " % (self.Skip_Dist_D),file=output)
            print("Normalized skip distance Ls/vf[s] :\t%.5f " % (self.Skip_Dist_V),file=output)
            print("Primary slot Blockage Factor :\t\t%.2f " % (self.BLF1),file=output)
            print("Secondary slot Blockage Factor :\t%.2f " % (self.BLF2),file=output)
            print("Oxidizer Injector delta p [MPa]:\t%.2f " % (self.deltap_o),file=output)
            print("Tertiary slot flow rate [kg/s]:\t\t%.3f " % (self.mdot_Tert_o),file=output)
            print("Tertiary slot flow Dia [mm]:\t\t%.3f " % (self.Tert_Dia),file=output)
            print("LOx Manifold Velocity [m/s] :\t\t%.1f " % (self.man_lox),file=output)
            print("Fuel Manifold Velocity [m/s] :\t\t%.1f " % (self.man_fuel),file=output)
            if self.rho_f > 161.383 : 
                print("Fuel is uncompressed",file=output)
                print("Fuel Injector delta p [MPa]:\t\t%.2f " % (self.deltap_f),file=output)
                print("Fuel Injector Downstream Pressure [MPa]:%.2f " % (self.P_fdo),file=output)
            else:
            
                if self.P_fdo > self.P_fup/2:
                    print("Fuel is compressed and non-choked",file=output)
                    print("Fuel Injector delta p [MPa]:\t\t%.2f " % (self.deltap_f),file=output)
                    print("Fuel Injector Downstream Pressure [MPa]:%.2f " % (self.P_fdo),file=output)
                else:
                    print("Fuel is choked",file=output)
                    print("Fuel Injector delta p is less than " + str(self.P_fup/2) +"[MPa]",file=output)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        setting_file = 'setting.ini'
    else:
        setting_file = sys.argv[1]
        assert os.path.exists(setting_file), "ファイルが存在しません"
    plt.close("all")
    plt.ion()
    pout = Pintle(setting_file)
    pout.display()
    pout.print()
