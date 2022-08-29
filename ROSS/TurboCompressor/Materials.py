#ROSS Material Definition File
#Initiate First 

#ROSS API Reference:
#https://ross.readthedocs.io/en/latest/references/api.html

import os
from pathlib import Path
import ross as rs
import numpy as np

#Define Materials===============================================
steel = rs.Material(name="Steel", rho=7810, E=211e9, G_s=81.2e9)
steel.save_material()

ChroMoly1 = rs.Material(name="CrMo4340", rho=7900, E=212e9, G_s=82.2e9)
ChroMoly1.save_material()

#Aluminum Alloy
AC1 = rs.Material(name="A2618", rho=2780, E=78e9, Poisson=0.27)
AC1.save_material()

I713C = rs.Material(name="I713C", rho=7980, E=211e9, G_s=82.0e9)
I713C.save_material()


