"""
t3d results plotting for Mehfoos MSc Dissertation
"""

### IMPORTS ###
from t3d.tools.profile_plot import plotter
import os
import numpy as np
###
    
### IMPORT AND EXTRACT DATA

def merge_defaultdicts(d,d1): # From SO
    for k,v in d1.items():
        if (k in d):
            d[k].update(d1[k])
        else:
            d[k] = d1[k]
    return d

# Define full path
fullPath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/t3d_Plots/"
        
# Instantiating plotter class
t3dPlot = plotter()

t3dPlot.grid_lines = True

#t3dOutputFile1 = "input_step1.log.npy"
t3dOutputFiles = sorted([f for f in os.listdir(fullPath) if f.endswith('.npy')]) # From SO
print(t3dOutputFiles)

fullData = np.load(fullPath+t3dOutputFiles[0], allow_pickle=True).tolist()
for iFile in t3dOutputFiles:
    data = np.load(fullPath+iFile, allow_pickle=True).tolist()
    fullData = merge_defaultdicts(fullData,data)
    #fullData = np.concatenate([fullData,data])

print('debug')

#t3dPlot.read_data(t3dOutputFile)

#t3dPlot.plot_panel()

