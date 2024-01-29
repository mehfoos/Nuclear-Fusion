"""
t3d results plotting for Mehfoos MSc Dissertation
"""

### IMPORTS ###
from profile_plot_MY1 import plotter
#from profile_plot_MY2 import plotter
from t3d.tools.profile_plot import plotter as t3dPlotter
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from netCDF4 import Dataset
###
    

# Define full path
fullPath = "C:\\Users\\mly509\\gitRepo1\\Nuclear-Fusion\\MscProject\\MLY\\t3d_Results\\"
        
singleCase = False

ncf = Dataset(fullPath+"profiles.cdf", format="NETCDF4")

# Plot just one file
if singleCase:
    t3dPlot = t3dPlotter()
    #t3dPlot = plotter()
    t3dPlot.grid_lines = True
    t3dOutputFile = "input_step1.log.npy"
    #t3dOutputFile = "profiles.cdf"
    t3dPlot.read_data(fullPath + t3dOutputFile)
    #t3dPlot.plot_panel()
    
    fig, axs = plt.subplots(1, 1)
    t3dPlot.plot_power_balance(axs)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.2, left=0.2)
    #t3dPlot.plot_density_profiles(axs)
    
    #_,axs2 = plt.subplots(1, 1)
    #t3dPlot.plot_state_profiles(axs2, profile='ni')

    fig, axs = plt.subplots(1, 1)
    plt.plot(t3dPlot.beta)
    
    plt.show()

# Plot multiple files
if not(singleCase):
    t3dPlots = []
    t3dOutputFiles = sorted([f for f in os.listdir(fullPath) if f.endswith('.npy')]) # From SO
    print(t3dOutputFiles)

    #fullData = np.load(fullPath+t3dOutputFiles[0], allow_pickle=True).tolist()
    
    # Read in all data
    for iFile in range(len(t3dOutputFiles)):
        t3dPlots.append(plotter())
        t3dPlots[iFile].grid_lines = True
        t3dPlots[iFile].read_data(fullPath +t3dOutputFiles[iFile])
        #t3dPlots[0].plot_panel()
        #data = np.load(fullPath+iFile, allow_pickle=True).tolist()
        #fullData = merge_defaultdicts(fullData,data)
        #fullData = np.concatenate([fullData,data])

    # Calculate total time and timesteps
    simTimeTotal = 0
    simStepTotal = 0
    for iFile in range(len(t3dOutputFiles)):
        simTimeTotal += t3dPlots[iFile].time[-1]
        simStepTotal += t3dPlots[iFile].N
        pass
    
    print('simTimeTotal [s]: ', simTimeTotal)
    print('simStepTotal [#]: ', simStepTotal)

    # Plot something on same figure
    #fig, axs = plt.subplots(1, 1)
    #fig, axs = plt.subplots(4, 4, figsize=(, 8))
    fig, axs = plt.subplots(4, 3)
    fig.set_figheight(20)
    fig.set_figwidth(10)
    N_cumulative_last = 0 # The acumulative "N" up to certain log
    time_cumulative_last = 0 # Similar, but for time
    for iFile in range(len(t3dOutputFiles)):
        if iFile == 0:
            t3dPlots[iFile].global_first = True # flag if first sim of a set
        else:
            t3dPlots[iFile].global_first = False
        if iFile == len(t3dOutputFiles)-1:
            t3dPlots[iFile].global_last = True # flag if last sim of a set
        else:
            t3dPlots[iFile].global_last = False

        t3dPlots[iFile].N_cumulative_last = N_cumulative_last
        t3dPlots[iFile].time_cumulative_last = time_cumulative_last
        
        t3dPlots[iFile].warm_map = pylab.cm.autumn(np.linspace(1, 0.25, simStepTotal))
        t3dPlots[iFile].cool_map = pylab.cm.Blues(np.linspace(0.25, 1, simStepTotal))
        t3dPlots[iFile].green_map = pylab.cm.YlGn(np.linspace(0.25, 1, simStepTotal))
        
        #t3dPlots[iFile].plot_gamma(axs)
        t3dPlots[iFile].plot_panel_MY(axs)

        N_cumulative_last += t3dPlots[iFile].N
        time_cumulative_last += t3dPlots[iFile].time[-1]
        
    
    plt.tight_layout(rect=(0,0.02,1,0.99))
    plt.subplots_adjust(hspace=0.6)
    #fig.delaxes(axs[4,2])
    # plt.subplots_adjust(left=0.05,
    #                     bottom=0.1,
    #                     right=0.95,
    #                     top=0.9,
    #                     wspace=0.4,
    #                     hspace=0.4)    
    plt.show()

print('debug')

