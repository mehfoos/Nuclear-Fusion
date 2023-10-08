#!/usr/bin/env python3
"""
GX linear CBC results replication using GS2
(Tokamak benchmark: Cyclone base case; kinetic electrons)
"""

__author__ = "Mehfoos"
__version__ = "0.1.0"
__license__ = "MIT"

from pyrokinetics import Pyro, template_dir, normalisation
import matplotlib.pyplot as plt
import math

def main():
    """ Main entry point of the app """ 
    
    ### IMPORT AND EXTRACT DATA

    # Path configurations
    inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/gxLinCbc2/gxLinCbc2.in"
        
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

    ### FILTER /PROCESS DATA
    
    # Get the final growth rate and 
    growth_rate_final = growth_rate.isel(kx=0).max('time')
    
    # Get the final growth rate adjusted/scaled as in GX (including coordinates)
    growth_rate_finalGX = growth_rate_final*math.sqrt(2)
    growth_rate_finalGX = growth_rate_finalGX.assign_coords(ky=growth_rate_finalGX.ky/math.sqrt(2))

    # Check if equilibrium has been reached
    growth_rate_4thQaurterMean = growth_rate.isel(kx=0).where(growth_rate.time>growth_rate.time[-1]/4).mean('time')
    growth_rate_maxChange4thQaurterMeanToFinalPct = abs(100*(growth_rate_final-growth_rate_4thQaurterMean)/growth_rate_final).max('ky')
    if growth_rate_maxChange4thQaurterMeanToFinalPct<10:
        print("Growth rate likely reached equilibrium. Worst case deviation of last quarter of data-series to final value is: ","%.2f" %  growth_rate_maxChange4thQaurterMeanToFinalPct.item(), "%")
    else:
        print("Growth rate likely did not get to equilibrium. Worst case deviation of last quarter of data-series to final value is: ","%.2f" %  growth_rate_maxChange4thQaurterMeanToFinalPct.item(), "%")

    ### PLOT DATA

    ## Plot all

    plt.figure(1)
    growth_rate.isel(kx=0).plot.line(x='time') 

    ## Plot the last time series, without scaling
    plt.figure(2)
    # Extract growth rate where kx=0 and time is the max of each series and plot against ky.
    growth_rate_final.plot.line(x='ky')
    plt.suptitle('Growth rate where kx=0 and time is the max of each series')

    ## Plot the last time series, with scaling
    plt.figure(3)
    (growth_rate_finalGX).plot.line(x='ky',color="orange",marker="s")
    plt.suptitle('Replication attempt: GX Figure2 a) Top')
    plt.xlim((0.0,1.7))
    plt.ylim((0.0,0.28))
    plt.xticks((0.00,0.25,0.50,0.75,1.00,1.25,1.50))
    plt.yticks((0.00,0.1,0.2))
    plt.minorticks_on()

    plt.show()

if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()
