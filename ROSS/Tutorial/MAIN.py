#MAIN Calculation
#MIT License: https://github.com/ross-rotordynamics/ross
#Run after Materials, Shaft, Disk, Bearing_Seal,PointMass and Rotor
#DEFAULT SETTING: Timoshenko beam Theory (rotary_inertia=True, shear_effects=True - used as default)

#Initialize=================================================================
import os
from pathlib import Path
import ross as rs
import numpy as np
import plotly.express as px
import Rotor

rotor=Rotor.rotor

#OPTION Load ROTOR==========================================================
#rotor = rs.Rotor.load('./output/rotor.toml')


#CALCULATIONS & PLOT========================================================
#1.1.Static Analysis========================================================
#static = rotor.run_static()
#print(rotor.bearing_forces_nodal)
#deform = static.plot_deformation()
#deform.show()


#1.2.Modal Analysis=========================================================
rotor_speed = 2500.0 # rad/s
modal = rotor.run_modal(rotor_speed)
print(f"Undamped natural frequencies[rad/s, rpm]:\n {modal.wn}, {modal.wn*60/np.pi/2}")
print()
print(f"Damped natural frequencies[rad/s, rpm]:\n {modal.wd}, {modal.wd*60/np.pi/2}")
print()
print(f"Damping ratio for each mode:\n {modal.damping_ratio}")
print()
print(f"Log Decrement for each mode:\n {modal.log_dec}")
print()
print(f"Whirl Direction for each mode:\n {modal.whirl_values()}")

dampfig = px.line(x=modal.wn, y=modal.damping_ratio, title="Damping ratio", labels={'x':'rad/s', 'y':'Damping Ratio'})

logdecfig = px.line(x=modal.wn, y=modal.log_dec, title="Log Decrement Stability Map", labels={'x':'rad/s', 'y':'Log Decrement'})
logdecfig.show()
logdecfig.write_html("./output/logdecfig_pic.html")

whirlfig = px.scatter(x=modal.wn, y=modal.whirl_values(), title="Whirl Direction", labels={'x':'rad/s', 'y':'Whirl Direction: Forward:0, Backward:0.5, Mixed:1.0'})
whirlfig.show()
whirlfig.write_html("./output/whirlfig_pic.html")

mode = 1
while mode < 6:
	modeplot3 = modal.plot_mode_2d(mode)
	modeplot3d3 = modal.plot_mode_3d(mode)
	modeplot3d3.show()
	modeplot3.write_html("./output/mode2d_mode"+str(mode)+".html")
	modeplot3d3.write_html("./output/mode3d_mode"+str(mode)+".html")
	mode = mode +1



#1.3 Campbell Plot===========================================================
samples = 20
speed_range = np.linspace(0, 2500, samples)
campbell_res = rotor.run_campbell(speed_range)
campbell = campbell_res.plot(harmonics=[1.0,2.0,3.0], frequency_units="RPM")
campbell.show()
campbell.write_html("./output/campbell_pic.html")


#1.4 Frequency Response======================================================
samples = 200
nodenum = 0   #specify node number
localdof = 0  #(x=0,y=1,alpha=2,beta=3) #specify direction in each node
bodeindex = nodenum*4+localdof
speed_range = np.linspace(0, 2500, samples) # rads/s, samples
freqres = rotor.run_freq_response(speed_range=speed_range)
freqplot = freqres.plot(inp=bodeindex, out=bodeindex,frequency_units="RPM")
freqplot.show()
freqplot.write_html("./output/freqplot_pic.html")


#1.5 Unbalance Response=======================================================
unbalnodes = [1,20] #Node no.
amps = [0.0001,0.0001] #Amplitude
phase = [0,np.pi*0.5] #Phase in radians
ubspeed = 1908 #in rad/s
frequency_range=np.linspace(0, ubspeed, 100)

ubalres = rotor.run_unbalance_response(unbalnodes, amps, phase, frequency_range)

# probe = (probe_node, probe_orientation)
probe1 = (1, 45) # node 4, orientation 45º
probe2 = (20, 45) # node 9, orientation 45º

unbalplot = ubalres.plot(probe=[probe1, probe2], probe_units="degrees",frequency_units="RPM")
unbalplot.show()
unbalplot.write_html("./output/unbalplot_pic.html")

unbalshape = ubalres.plot_deflected_shape(speed=ubspeed,frequency_units="RPM")  #speed input in rad/s
unbalshape.show()
unbalshape.write_html("./output/unbalshape_pic.html")


#1.6 Time response external force===============================================
speed = 1000 #rad/s
time_samples = 1000
node = 6
t = np.linspace(0, 10, time_samples)
F = np.zeros((time_samples, rotor.ndof))  #No DoF is rotor
# component on direction x
F[:, 4 * node + 0] = 5 * np.cos(2*t)  #Defult harmonic: t*2*np.pi
# component on direction y
F[:, 4 * node + 1] = 5 * np.sin(2*t)
response = rotor.run_time_response(speed, F, t)

probe1 = (10, 0)   # node 3, orientation 0° (X dir.)
probe2 = (10, 90)  # node 3, orientation 90°(Y dir.)
timeres1 = response.plot_1d(probe=[probe1, probe2], probe_units="degree")
timeres2 = response.plot_2d(node=node)
timeres3 = response.plot_3d()
timeres1.show()
timeres2.show()
timeres3.show()
timeres1.write_html("./output/timeres1_pic.html")
timeres2.write_html("./output/timeres2_pic.html")
timeres3.write_html("./output/timeres3_pic.html")


#1.7 UCS Undamped Critical Speed Map==============================================
stiff_range = (3, 11)  #10e6 to 10e11 N/m shown as 6,11
ucs_results = rotor.run_ucs(stiffness_range=stiff_range, num=20, num_modes=16, frequency_units="RPM")
ucs_fig = ucs_results.plot()
ucs_fig.show()
ucs_fig.write_html("./output/ucs_fig.html")

