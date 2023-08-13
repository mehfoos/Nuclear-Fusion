#!/usr/bin/env python3
"""
Test/Play with Pyrokinetics + GS2
"""

__author__ = "Mehfoos"
__version__ = "0.1.0"
__license__ = "MIT"

from pyrokinetics import Pyro, template_dir

def main():
    """ Main entry point of the app """ 

    # Path configurations
    # inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/gxLinCbc1/gxLinCbc1.in"
    inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/gxLinCbc2/gxLinCbc2.in"
    # inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/test1/input_MG_MY.in"
    # outputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/gxLinCbc1/gxLinCbc1.out.nc"
    # gs2_template = template_dir / "outputs/GS2_linear/gs2.in"
        
    # instantiating the Pyro class as pyro object.
    pyro = Pyro(gk_file=inputFilePath, gk_code="GS2")
    
    # Checking if current gyrokinetics context is valid. https://pyrokinetics.readthedocs.io/en/latest/generated/pyrokinetics.pyro.Pyro.html#pyrokinetics.pyro.Pyro.check_gk_code
    print("GK code check pass?: ", pyro.check_gk_code())

    # Load GS2 output data
    pyro.load_gk_output()

    # Print all the data
    print(pyro.gk_output.data_vars)

    # Extract signals
    ky = pyro.gk_output["ky"]
    gamma = pyro.gk_output["gamma"]

    # Plot signals
    


    #

if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()