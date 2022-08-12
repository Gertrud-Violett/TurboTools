#Rotor Creation definition
#MIT License: https://github.com/ross-rotordynamics/ross
#Run after Materials, Shaft, Disk, Bearing_Seal and PointMass
#DEFAULT SETTING: Timoshenko beam Theory (rotary_inertia=True, shear_effects=True - used as default)

#Initialize=================================================================
import os
from pathlib import Path
import ross as rs
import numpy as np
import Materials
import Shaft
import Disk
import Bearing_Seal
import PointMass

shafts=Shaft.shafts
disks=Disk.disks
bearings=Bearing_Seal.bearings
pointmass=PointMass.pointmass

#VVV SAMPLE CODE VVV SHOWN BELOW
#ROTOR DEFINITIONS==========================================================
#CREATE ROTOR===============================================================
#REGULAR ELEMENT METHOD=====================================================
"""
rotor1 = rs.Rotor(shafts, disks, bearings, pointmass)

print("Rotor total mass = ", np.round(rotor1.m, 2))
print("Rotor center of gravity =", np.round(rotor1.CG, 2))

# plotting the rotor model==================================================
rotorfig1 = rotor1.plot_rotor()
#rotorfig.update_yaxes(range=[-0.2, 0.2])
rotorfig1.show()
rotorfig1.write_html("rotor1_pic.html")

#BY SECTION METHOD==========================================================
rotor2 = rs.Rotor.from_section(
    brg_seal_data=bearings,
    disk_data=disks,
    idl_data=i_ds_data,
    leng_data=leng_data,
    odl_data=o_ds_data, 
    nel_r=4,
    material_data=steel,
)

"""

#VVV WRITE CODE HEREVVV
#ROTOR DEFINITIONS==========================================================
rotor = rs.Rotor(shafts, disks, bearings)
#rotor = rs.Rotor(shafts, disks, bearings, pointmass) #for n_link option

# plotting the rotor model==================================================
rotorfig = rotor.plot_rotor()
#rotorfig.show()
rotorfig.write_html("./output/rotor_pic.html")

#Save & Load ROTOR==========================================================
rotor.save('./output/rotor.toml')
#rotor1 = rs.Rotor.load('./output/rotor.toml')
#rotor == rotor1


