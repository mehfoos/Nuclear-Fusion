#!/usr/bin/env python3
"""
GX results replication using GS2
"""

__author__ = "Mehfoos"
__version__ = "0.1.0"
__license__ = "MIT"

from pyrokinetics import Pyro, template_dir, normalisation
import matplotlib.pyplot as plt
import math

def main():
    """ Main entry point of the app """ 

    # Path configurations
    # inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/gxLinCbc1/gxLinCbc1.in"
    inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/gxLinCbc2/gxLinCbc2.in"
        
    # instantiating the Pyro class as pyro object.
    pyro = Pyro(gk_file=inputFilePath, gk_code="GS2")
    
    # Checking if current gyrokinetics context is valid. https://pyrokinetics.readthedocs.io/en/latest/generated/pyrokinetics.pyro.Pyro.html#pyrokinetics.pyro.Pyro.check_gk_code
    print("GK code check pass?: ", pyro.check_gk_code())

    # Load GS2 output data
    pyro.load_gk_output()

    # Print all the data
    #print(pyro.gk_output.data_vars)

    # Extract signals
    ky = pyro.gk_output["ky"]
    eigenvalues = pyro.gk_output["eigenvalues"]
    growth_rate = pyro.gk_output["growth_rate"]
    mode_freq = pyro.gk_output["mode_frequency"]

    # Plot signals
    plt.figure(1)
    growth_rate.isel(kx=0,time=growth_rate.time.size-1).plot.line(x='ky')
    
    ## Plot after scaling

    # Convert to np array
    growth_rate_np = growth_rate.isel(kx=0,time=growth_rate.time.size-1).to_numpy()
    ky_np = growth_rate.ky.to_numpy()

    # Scale to GX paper levels
    growth_rate_np_GX = growth_rate_np/math.sqrt(2)
    ky_np_GX = ky_np/math.sqrt(2)
    
    # Plot
    plt.figure(2)
    plt.plot(ky_np_GX,growth_rate_np_GX)
    plt.ylabel('growth_rate_np_GX')
    plt.xlabel('ky_np_GX')
    plt.suptitle('Normalised growth rate vs binormal wavenumber')
    plt.title('at kx=0; last time value')


    plt.show()

if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()
