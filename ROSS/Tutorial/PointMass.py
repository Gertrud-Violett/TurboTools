#Point Mass geometry definition
#MIT License: https://github.com/ross-rotordynamics/ross
#Run after Materials.py

#Initialize=================================================================
import os
from pathlib import Path
import ross as rs
import numpy as np

#VVV SAMPLE CODE VVV SHOWN BELOW
#POINT MASS DEFINITIONS====================================================
#SAMPLES===================================================================
# Single mass element 1====================================================
"""
p0 = rs.PointMass(n=0, m=2)
p1 = rs.PointMass(n=0, mx=2, my=3)
#print(p1.M())
"""

#VVV WRITE CODE HEREVVV
#POINT MASS DEFINITIONS====================================================
pm0 = rs.PointMass(n=20, m=0.1)
pointmass = [pm0]