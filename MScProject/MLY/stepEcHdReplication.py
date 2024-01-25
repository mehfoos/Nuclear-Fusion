"""
Replication attempt for:
Kennedy 2023 Electromagnetic gyrokinetic instabilities in STEP 
Figure 2
"""

### IMPORTS ###
from pyrokinetics import Pyro, template_dir, normalisation
import matplotlib.pyplot as plt
import math
import pandas as pd
###
    
### IMPORT AND EXTRACT DATA

# Path configurations
inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/stepEcHd_2/input.in"
        
# instantiating the Pyro class as pyro object.
pyro = Pyro(gk_file=inputFilePath, gk_code="GS2")

# Checking if current gyrokinetics context is valid. 
print("GK code check pass?: ", pyro.check_gk_code())

# Load GS2 output data
pyro.load_gk_output()

# Extract signals
ky = pyro.gk_output["ky"]
eigenvalues = pyro.gk_output["eigenvalues"]
growth_rate = pyro.gk_output["growth_rate"]
mode_freq = pyro.gk_output["mode_frequency"]
###

### FILTER /PROCESS DATA

# Get the final growth rate and frequencies
growth_rate_final = growth_rate.isel(kx=0,time=-1)
mode_freq_final = mode_freq.isel(kx=0,time=-1)

# Get the final growth rate adjusted/scaled as in GX (including coordinates)
growth_rate_finalGX = growth_rate_final #*math.sqrt(2)
# growth_rate_finalGX = growth_rate_finalGX.assign_coords(ky=growth_rate_finalGX.ky/math.sqrt(2))
mode_freq_finalGX = mode_freq_final #*math.sqrt(2)

# Check if equilibrium has been reached
growth_rate_4thQaurterMean = growth_rate.isel(kx=0).where(growth_rate.time>growth_rate.time[-1]/4).mean('time')
growth_rate_maxChange4thQaurterMeanToFinalPct = abs(100*(growth_rate_final-growth_rate_4thQaurterMean)/growth_rate_final).max('ky')
if growth_rate_maxChange4thQaurterMeanToFinalPct<10:
    print("Growth rate likely reached equilibrium. Worst case deviation of last quarter of data-series to final value is: ","%.2f" %  growth_rate_maxChange4thQaurterMeanToFinalPct.item(), "%")
else:
    print("Growth rate likely did not get to equilibrium. Worst case deviation of last quarter of data-series to final value is: ","%.2f" %  growth_rate_maxChange4thQaurterMeanToFinalPct.item(), "%")
###

### STEP paper data
DataFig2A = pd.read_csv('MScProject/MLY/Step1Fig2a_psi0.49_growth.csv',skipinitialspace=True)
DataFig2B = pd.read_csv('MScProject/MLY/Step1Fig2b_psi0.49_modefreq.csv',skipinitialspace=True)
###

### PLOT DATA

## Plot all

fig, ax = plt.subplots()
growth_rate.isel(kx=0).plot.line(x='time') 

fig, ax = plt.subplots()
mode_freq.isel(kx=0).plot.line(x='time') 
##

## Plot the last time series, without scaling

fig, ax = plt.subplots()
# Extract growth rate where kx=0 and time is the max of each series and plot against ky.
(1*growth_rate_final).plot.line(x='ky',marker="o",label="GS2 Replication", ax=ax)
DataFig2A.plot(x="x",y="y",label="psi=0.49 (STEP paper)", color="orange",marker="s", markerfacecolor='none', linestyle='None', ax=ax)
plt.suptitle('Growth rate where kx=0 and time is the max of each series')

fig, ax = plt.subplots()
# Extract growth rate where kx=0 and time is the max of each series and plot against ky.
(1*mode_freq_final).plot.line(x='ky', ax=ax)
DataFig2B.plot(x="x",y="y",label="psi=0.49 (STEP paper)", color="orange",marker="s", markerfacecolor='none', linestyle='None', ax=ax)
plt.suptitle('Mode Frequencies where kx=0 and time is the max of each series')
##

plt.show()
###
