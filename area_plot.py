# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 10:45:46 2015

@author: Sierra Anderson
"""
import sys
import numpy as np
from matplotlib import pyplot as plt

MAX_DATA_POINTS = 350

def plot(profile, parameters):
    if (profile.abundance_data.shape[1] > MAX_DATA_POINTS):
        print("ERROR: Too many data points to generate area plot.")        
        sys.exit(0)
    
    profile.abundance_data.plot(kind='bar',stacked=True)
    
    

