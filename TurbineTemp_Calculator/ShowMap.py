"""
===-*- Turbine Inlet Temperature Correlation Tool -*-===
ShowMap Module for Turbine Inlet Temperature Interpolation

=====-*- General -*-=====
Copyright (c) makkiblog.com
MIT License 
coding: utf-8

===-*- VERSION -*-===
v0.1 Initial version
v1.0 Released

vvvCODEvvv
"""
import matplotlib.pyplot as plt 
import pickle


def show():
    fig = pickle.load(open("./Output/T1t_map.pickle", "rb") )
    fig.show()

