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
import matplotlib.scale as mscale
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
DataFig2AIncStable = pd.read_csv('MScProject/MLY/Step1Fig2a_psi0.49_growth.csv',skipinitialspace=True)
DataFig2A = pd.read_csv('MScProject/MLY/Step1Fig2a_psi0.49_growth unstable.csv',skipinitialspace=True)
DataFig2B = pd.read_csv('MScProject/MLY/Step1Fig2b_psi0.49_modefreq.csv',skipinitialspace=True)
###

### PLOT DATA

## Plot all

#fig, ax = plt.subplots()
#growth_rate.isel(kx=0).plot.line(x='time') 

#fig, ax = plt.subplots()
#mode_freq.isel(kx=0).plot.line(x='time') 
##

## Plot the last time series, without scaling

#fig, ax = plt.subplots()
# Extract growth rate where kx=0 and time is the max of each series and plot against ky.
#(1*growth_rate_final).plot.line(x='ky', color="black",marker="o",label="GS2 Replication", ax=ax)
#DataFig2A.plot(x="x",y="y",label="psi=0.49 (STEP paper)", color="orange",marker="s", markerfacecolor='none', linestyle='None', ax=ax)
#plt.suptitle('Growth rate where kx=0 and time is the max of each series')

#fig, ax = plt.subplots()
# Extract growth rate where kx=0 and time is the max of each series and plot against ky.
#(1*mode_freq_final).plot.line(x='ky', ax=ax)
#DataFig2B.plot(x="x",y="y",label="psi=0.49 (STEP paper)", color="orange",marker="s", markerfacecolor='none', linestyle='None', ax=ax)
#plt.suptitle('Mode Frequencies where kx=0 and time is the max of each series')
##

## Plot GX paper figure 2.a.top
fig, axs = plt.subplots(2,1,layout="constrained")
DataFig2AIncStable.plot(x="x",y="y",label="Baseline; stable modes", color="orange",marker="s", markerfacecolor='none', linestyle='None', ax=axs[0])
DataFig2A.plot(x="x",y="y",label="Baseline; unstable modes", color="orange",marker="s", markerfacecolor='orange', linestyle='None', ax=axs[0])
(1*growth_rate_final).plot.line(x='ky', color="black",marker="x",label="Replication", ax=axs[0])
plt.suptitle('Replication attempt: STEP-EC-HD Ψ=0.49')
axs[0].set_title('')
axs[0].set_xlabel('')
axs[0].set_ylabel('Normalised growth rate\n $\\gamma a/c_{s}$')
axs[0].legend()
axs[0].set_xscale('log')
axs[0].grid(linestyle=':')
#plt.ylim((0.0,0.28))
#axs[0].set_yticks((0.00,0.1,0.2))
#axs[0].minorticks_on()

## Plot GX paper figure 2.a.bottom
DataFig2B.plot(x="x",y="y",label="Baseline; unstable modes", color="orange",marker="s", markerfacecolor='orange', linestyle='None', ax=axs[1])
(1*mode_freq_final).plot.line(x='ky',color="black",marker="x",label="Replication", ax=axs[1])
axs[1].sharex(axs[0])
axs[1].set_xscale('log')
axs[1].legend()
plt.xlabel('Normalised binormal wavenumber $k_y\\rho_i$')
axs[1].set_ylabel('Real (Mode) Frequencies\n $\\omega a/c_{s}$')
plt.title('Real (Mode) Frequencies')
axs[1].set_title('')
#plt.xlim((0.0,1.7))
#plt.ylim((-1,1))
#plt.xticks((0.00,0.25,0.50,0.75,1.00,1.25,1.50))
#plt.yticks((-1,-0.5,0,0.5,1))
#plt.minorticks_on()
fig.align_ylabels()
axs[1].grid(linestyle=':')

plt.show()
###
