#Bearing and Seal definition
#MIT License: https://github.com/ross-rotordynamics/ross
#Run after Materials.py

#Initialize=================================================================
import os
from pathlib import Path
import ross as rs
import numpy as np
import ross.bearing_seal_element
import plotly.graph_objects as go

#VVV SAMPLE CODE VVV SHOWN BELOW
#BEARING DEFINITIONS========================================================
#SAMPLES====================================================================
# Constant Coefficient Bearing==============================================
"""
stfx = 1e6
stfy = 0.8e6
bearing1 = rs.BearingElement(n=0, kxx=stfx, kyy=stfy, cxx=1e3)

#print("="*55)
#print(f"Kxx coefficient: {bearing1.kxx}")
"""

# Speed Dependent Coefficient Bearing========================================

bearing5 = rs.BearingElement(
    n=0, 
    kxx=[0.5e6, 1.0e6, 2.5e6],
    kyy=[1.5e6, 2.0e6, 3.5e6],
    cxx=[0.5e3, 1.0e3, 1.5e3],
    frequency=[0, 1000, 2000],
)

#print("="*79)
#print(f"Kxx coefficient: {bearing5.kxx}")

"""

# n_link series bearing: two elements in one node============================
# Use with list for speed dependent coefficients n_link=1,2... for link
stfx = 1e6
stfy = 0.8e6
bearing3 = rs.BearingElement(n=0, kxx=stfx, kyy=stfy, cxx=1e3, n_link=1, tag="journal_bearing")
bearing4 = rs.BearingElement(n=1, kxx=1e7, kyy=1e9, cxx=10, tag="support")
"""

#BALL BEARING DEFINITIONS===================================================
#SAMPLES====================================================================
# Constant Coefficient BALL Bearing=========================================
"""
n = 0
n_balls= 8
d_balls = 0.03
fs = 500.0
alpha = np.pi / 6
tag = "ballbearing"
ballbearing = rs.BallBearingElement(
    n=n,
    n_balls=n_balls,
    d_balls=d_balls,
    fs=fs,
    alpha=alpha,
    tag=tag,
)

# Constant Coefficient ROLLER Bearing=======================================
n = 0
n_rollers= 8
l_rollers = 0.03
fs = 500.0
alpha = np.pi / 6
tag = "rollerbearing"
rollerbearing = rs.RollerBearingElement(
    n=n,
    n_rollers=n_rollers,
    l_rollers=l_rollers,
    fs=fs,
    alpha=alpha,
    tag=tag
)

# Magnetic Bearing =========================================================
#Refer to ROSS documentations

#EXCEL OPTION: Read bearing elemets from excel file=========================
file_path = Path("bearing_seal_si.xls")
bearing = rs.BearingElement.from_table(n=0, file=file_path)
"""

#Create Bearing Array=======================================================
#bearings = [bearing0, bearing1, bearing2]


#VVV WRITE CODE HEREVVV
#BEARING DEFINITIONS========================================================
stfx = 1e6
stfy = 0.8e6
#bearing0 = rs.BearingElement(n=3, kxx=stfx, kyy=stfy, cxx=1e2, n_link=None)

bearing0 = rs.BearingElement(
    n=3, 
    kxx=[1.0e6, 2.0e6, 3.0e6],
    kyy=[0.8e6, 1.8e6, 2.8e6],
    cxx=[1.0e3, 1.5e3, 2.0e3],
    frequency=[0, 1000, 2000],
    n_link=None
)

bearing1 = rs.BearingElement(
    n=16, 
    kxx=[5.0e6, 10.0e6, 15.0e6],
    kyy=[5.0e6, 10.0e6, 15.0e6],
    cxx=[1.0e3, 1.5e3, 2.0e3],
    frequency=[0, 1000, 2000],
)


#bearing2 = rs.BearingElement(n=31, kxx=stfx*100, kyy=stfy*100, cxx=0)

#bearings = [bearing0, bearing1]
#bearings = [bearing0, bearing1,bearing2]

#BEARING CHECK==============================================================
#Check Bearing specs========================================================
bearing_fig=bearing0.kxx.plot() #Enter any bearing no and attribute
#bearing_fig.show()
bearing_fig.write_html("./output/bearing_pic.html")

#SEAL DEFINITIONS===========================================================
#DEFINITION SEAL============================================================
# Constant Coefficient Seals================================================
stfx = 1e2
stfy = 0.8e2
seal1 = rs.SealElement(n=2, kxx=stfx, kyy=stfy, cxx=1e3, cyy=0.8e3)
seal2 = rs.SealElement(n=17, kxx=stfx, kyy=stfy, cxx=1e3, cyy=0.8e3)
seal3 = rs.SealElement(n=9, kxx=stfx, kyy=stfy, cxx=1e1, cyy=0.8e1)


bearings = [bearing0, bearing1, seal1, seal2, seal3]