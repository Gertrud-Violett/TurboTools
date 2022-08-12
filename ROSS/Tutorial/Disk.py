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
#DISK DEFINITIONS==========================================================
#SAMPLES===================================================================
# Single disk element 1====================================================
"""
disk_elem = rs.DiskElement.from_geometry(
    n=4,  #Node no.
    material=steel,
    width=0.07,
    i_d=0.05,
    o_d=0.28
)

#print("="*76)
#print(f"Disk mass:              {disk1_Ex.m}")
#print(f"Disk polar inertia:     {disk1_Ex.Ip}")
#print(f"Disk diametral inertia: {disk1_Ex.Id}")

#Creating a list of disk elements==========================================
n_list = [2, 4]   #Node nos.
width_list = [0.7, 0.7]
i_d_list = [0.05, 0.05]
o_d_list = [0.15, 0.18]

N = len(n_list)
disk_elements = [
    rs.DiskElement.from_geometry(
        n=n_list[i],
        material=steel,
        width=width_list[i],
        i_d=i_d_list[i],
        o_d=o_d_list[i],
    )
    for i in range(N)
]

#OPTION with Id and Ip=====================================================
m_list = [85.875, 128.38]
Id_list =[3.6408, 5.5224]
Ip_list = [0.26836, 0.56007]
disk_elements_dir = [
    rs.DiskElement(
        n=n_list[i],
        m=m_list[i],
        Id=Id_list[i],
        Ip=Ip_list[i],
    )
    for i in range(N)
]

#print(disk_elements_dir)

#EXCEL OPTION: Read shaft elemets from excel file (Mass and Inertia only=====
#disk_file_path = Path("shaft_si.xls")
#list_of_disks_excel = rs.DiskElement.from_table(file=disk_file_path, sheet_name="More")

"""

#VVV WRITE CODE HEREVVV
#DISK=======================================================================
disk0 = rs.DiskElement.from_geometry(
    n=1, material=steel, width=0.20, i_d=0.10, o_d=0.23
)
disk1 = rs.DiskElement.from_geometry(
    n=20, material=steel, width=0.05, i_d=0.12, o_d=0.28
)
disk2 = rs.DiskElement.from_geometry(
    n=10, material=steel, width=0.3, i_d=0.12, o_d=0.20
)

#Create Disk Array==========================================================
disks = [disk0, disk1, disk2]




