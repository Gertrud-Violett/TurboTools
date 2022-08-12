#Shaft and Disk geometry definition
#MIT License: https://github.com/ross-rotordynamics/ross
#Run after Materials.py
#DEFAULT SETTING: Timoshenko beam Theory (rotary_inertia=True, shear_effects=True - used as default)

#Initialize and Load Materials=============================================
import os
from pathlib import Path
import ross as rs
import numpy as np
steel = rs.Material.load_material('Steel')

#VVV SAMPLE CODE VVV SHOWN BELOW
#SHAFT DEFINITIONS==========================================================
#SAMPLES====================================================================
# Cylindrical Shaft element 1===============================================
"""
L = 0.25
i_d = 0
o_d = 0.05
cy_elem = rs.ShaftElement(
    L=L,
    idl=i_d,
    odl=o_d,
    material=steel,
    shear_effects=True,
    rotary_inertia=True,
    gyroscopic=True,
)


# Conical Shaft element 1===================================================
L = 0.25
idl = 0
idr = 0
odl = 0.05
odr = 0.07
co_elem = rs.ShaftElement(
    L=L,
    idl=idl,
    idr=idr,
    odl=odl,
    odr=odr,
    material=steel,
    shear_effects=True,
    rotary_inertia=True,
    gyroscopic=True,
)


#Creating a list of shaft elements==========================================
L =   [0.20, 0.20, 0.10, 0.10, 0.20, 0.20]
i_d = [0.01,    0,    0,    0,    0, 0.01]
o_d = [0.05, 0.05, 0.06, 0.06, 0.05, 0.05]
N = len(L)  #count number of elements and add to node no.
shaft_elements = [
    rs.ShaftElement(
        L=L[i],
        idl=i_d[i],
        odl=o_d[i],
        material=steel,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for i in range(N)
]


#EXCEL OPTION: Read shaft elemets from excel file (Mass and Inertia only=====
#disk_file_path = Path("shaft_si.xls")
#list_of_disks_excel = rs.DiskElement.from_table(file=disk_file_path, sheet_name="More")
"""

#VVV WRITE CODE HEREVVV
#SHAFT=====================================================================
n = 10
shaft_elem = [
    rs.ShaftElement(
        L=0.05,
        idl=0.0,
        odl=0.1,
        material=steel,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for _ in range(n)
]

n = 10
shaft_elem2 = [
    rs.ShaftElement(
        L=0.05,
        idl=0.0,
        odl=0.12,
        material=steel,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for _ in range(n)
]

#Collar
shaft_elem3 = rs.ShaftElement(
        L=0.05,
        idl=0.12,
        odl=0.16,
        material=steel,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
        n=18
    )

shaft_elem4 = rs.ShaftElement(
        L=0.05,
        idl=0.12,
        odl=0.16,
        material=steel,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
        n=19
    )

shafts=[shaft_elem,shaft_elem2,shaft_elem3,shaft_elem4]

