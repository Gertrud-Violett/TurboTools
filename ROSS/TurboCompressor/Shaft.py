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
ChroMoly1 = rs.Material.load_material('CrMo4340')

#VVV WRITE CODE HEREVVV
#SHAFT=====================================================================
n = 7
shaft_elem = [
    rs.ShaftElement(
        L=0.00917/7,
        idl=0.0102,
        odl=0.01731,
        material=ChroMoly1,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for _ in range(n)
]

shaft_elem2 = [
    rs.ShaftElement(
        L=0.001156,
        idl=0.0,
        odl=0.01365,
        material=ChroMoly1,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
]


shaft_elem3 = [
    rs.ShaftElement(
        L=0.001156,
        idl=0.0,
        idr=0.0,
        odl=0.01365,
        odr=0.010156,
        material=ChroMoly1,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
]

n = 8
shaft_elem4 = [
    rs.ShaftElement(
        L=0.01063/8,
        idl=0.0,
        odl=0.010156,
        material=ChroMoly1,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for _ in range(n)
]

n = 10
shaft_elem5 = [
    rs.ShaftElement(
        L=0.02/10,
        idl=0.0,
        odl=0.009756,
        material=ChroMoly1,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for _ in range(n)
]

n = 8
shaft_elem6 = [
    rs.ShaftElement(
        L=0.012/9,
        idl=0.0,
        odl=0.010156,
        material=ChroMoly1,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for _ in range(n)
]

n = 5
shaft_elem7 = [
    rs.ShaftElement(
        L=0.009299/5,
        idl=0.0,
        odl=0.006529,
        material=ChroMoly1,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for _ in range(n)
]

n = 3
shaft_elem8 = [
    rs.ShaftElement(
        L=0.01342/3,
        idl=0.0,
        odl=0.006268,
        material=ChroMoly1,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for _ in range(n)
]


#Thrust Collar
ThrustCollar = rs.ShaftElement(
        L=0.0083,
        idl=0.006529,
        odl=0.0173,
        material=steel,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
        n=35
    )

shafts=[shaft_elem,shaft_elem2,shaft_elem3,shaft_elem4,shaft_elem5,shaft_elem6,shaft_elem7,shaft_elem8,ThrustCollar]

