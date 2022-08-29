#Point Mass geometry definition
#MIT License: https://github.com/ross-rotordynamics/ross
#Run after Materials.py

#Initialize=================================================================
import os
from pathlib import Path
import ross as rs
import numpy as np


#VVV WRITE CODE HEREVVV
#POINT MASS DEFINITIONS====================================================
#Rotor Locknut
pm0 = rs.PointMass(n=42, m=0.05)
pointmass = [pm0]