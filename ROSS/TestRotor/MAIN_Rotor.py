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
#BY EXCEL EXTERNAL FILE=====================================================
shaft_file = Path("shaft_si.xls")
shaft = rs.ShaftElement.from_table(
    file=shaft_file, sheet_type="Model", sheet_name="Model"
)

print(shaft)

file_path = Path("shaft_si.xls")
list_of_disks = rs.DiskElement.from_table(file=file_path, sheet_name="More")

bearing1 = rs.BearingElement.from_table(n=7, file="bearing_seal_si.xls", scale_factor=0.5)
bearing2 = rs.BearingElement.from_table(n=48, file="bearing_seal_si.xls", scale_factor=0.5)

bearings = [bearing1, bearing2]

rotor3 = rs.Rotor(shaft, list_of_disks, bearings)

# plotting the rotor model==================================================
rotorfig3 = rotor3.plot_rotor(nodes=2)
#rotorfig3.show()
rotorfig3.write_html("./output/rotor3_pic.html")


#Save & Load ROTOR==========================================================
rotor3.save('rotor3.toml')
#rotor3_1 = rs.Rotor.load('rotor3.toml')
#rotor3_1 == rotor3
#use excel2x


#CALCULATIONS & PLOT========================================================
#1.1.Static Analysis========================================================
"""
static = rotor3.run_static()
print(rotor3.bearing_forces_nodal)
deform = static.plot_deformation()
deform.show()
"""

#1.2.Modal Analysis=========================================================
rotor_speed = 100.0 # rad/s
modal = rotor3.run_modal(rotor_speed)
print(f"Undamped natural frequencies:\n {modal.wn}")
print()
print(f"Damped natural frequencies:\n {modal.wd}")
#print()
#print(f"Damping ratio for each mode:\n {modal.damping_ratio}")
mode = 5
modeplot1 = modal.plot_mode_3d(mode)
#modeplot1.show()
modeplot1.write_html("./output/mode3d_pic.html")

#1.3 Campbell Plot===========================================================
samples = 31
speed_range = np.linspace(0, 1100, samples)
campbell = rotor3.run_campbell(speed_range)
campbell3 = campbell.plot(harmonics=[1.0,2.0,3.0], frequency_units="RPM")
#campbell3.show()
#campbell3.write_html("./output/campbell3_pic.html")

#1.4 Frequency Response======================================================
nodenum = 43
localdof = 1  #(x=0,y=1,alpha=2,beta=3)
bodeindex = nodenum*4+localdof
speed_range = np.linspace(315, 1150, 31) # rads/s, samples
results3 = rotor3.run_freq_response(speed_range=speed_range)
#freqplot = results3.plot(inp=bodeindex, out=bodeindex,frequency_units="RPM")
#freqplot.show()

#1.5 Unbalance Response=======================================================
unbalnodes = [29, 33] #Node no.
amps = [0.003,0.002] #Amplitude
phase = [0,0] #Phase

frequency_range=np.linspace(315, 1150, 101)
results2 = rotor3.run_unbalance_response(unbalnodes, amps, phase, frequency_range)

# probe = (probe_node, probe_orientation)
probe1 = (15, 45) # node 15, orientation 45º
probe2 = (35, 45) # node 35, orientation 45º

#unbalplot = results2.plot(probe=[probe1, probe2], probe_units="degrees",frequency_units="RPM")
#unbalplot.show()

#unbalshape = results2.plot_deflected_shape(speed=649,frequency_units="RPM")
#unbalshape.show()

#1.6 Time response external force===============================================
speed = 100*2*np.pi #rad/s
time_samples = 1000
node = 26
t = np.linspace(0, 1, time_samples)
F = np.zeros((time_samples, rotor3.ndof))  #No DoF is rotor3
# component on direction x
F[:, 4 * node + 0] = 1 * np.cos(t*2*np.pi)
# component on direction y
F[:, 4 * node + 1] = 1 * np.sin(t*2*np.pi)
response3 = rotor3.run_time_response(speed, F, t)

probe1 = (3, 0)   # node 3, orientation 0° (X dir.)
probe2 = (3, 90)  # node 3, orientation 90°(Y dir.)
timeres1 = response3.plot_1d(probe=[probe1, probe2], probe_units="degree")
timeres2 = response3.plot_2d(node=node)
timeres3 = response3.plot_3d()
#timeres1.show()
#timeres2.show()
#timeres3.show()
timeres1.write_html("./output/timeres1_pic.html")
timeres2.write_html("./output/timeres2_pic.html")
timeres3.write_html("./output/timeres3_pic.html")


#1.7 UCS Undamped Critical Speed Map==============================================
stiff_range = (6, 11)  #10e6 to 10e11 N/m shown as 6,11
ucs_results = rotor3.run_ucs(stiffness_range=stiff_range, num=20, num_modes=16)
ucs_fig = ucs_results.plot()
#ucs_fig.show()