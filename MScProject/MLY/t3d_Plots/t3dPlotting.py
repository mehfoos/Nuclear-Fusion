"""
t3d results plotting for Mehfoos MSc Dissertation
"""

### IMPORTS ###
from t3d.tools.profile_plot import plotter
import os
import numpy as np
###
    
### IMPORT AND EXTRACT DATA

# Define full path
fullPath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/t3d_Plots/"
        
# Instantiating plotter class
t3dPlot = plotter()

t3dPlot.grid_lines = True

#t3dOutputFile1 = "input_step1.log.npy"
t3dOutputFiles = sorted([f for f in os.listdir(fullPath) if f.endswith('.npy')])
print(t3dOutputFiles)

fullData = []
for iFile in t3dOutputFiles:
    data = np.load(fullPath+iFile, allow_pickle=True)
    fullData.append(data)

print('debug')

#t3dPlot.read_data(t3dOutputFile)

#t3dPlot.plot_panel()

def merge_defaultdicts(d,d1):
    for k,v in d1.items():
        if (k in d):
            d[k].update(d1[k])
        else:
            d[k] = d1[k]
    return d