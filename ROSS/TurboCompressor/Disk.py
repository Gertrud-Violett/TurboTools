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
I713C = rs.Material.load_material('I713C')
AC1 = rs.Material.load_material('A2618')


#VVV WRITE CODE HEREVVV
#DISK=======================================================================
#Turbine
disk0 = rs.DiskElement(
    n=0, m=0.357,   Id=90e-6,Ip=117e-6,
)

#Compressor
disk1 = rs.DiskElement(
    n=40, m=0.139,Id=42e-6,Ip=51e-6,
)

#Create Disk Array==========================================================
disks = [disk0, disk1]




