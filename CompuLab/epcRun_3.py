# Code to import and run epc1d functions

from sys import stderr
from epc1d_Opt2 import *
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np
from matplotlib.pyplot import cm
from timeit import default_timer as timer
from scipy.stats import linregress


def repeatRuns()



if __name__ == "__main__":
    # Configuration
    # 2-stream instability
    L = 100
    ncells = 20
    vbeam_ = np.array(range(1,20,1)).astype(float)
    
    # Unique colors
    color = iter(cm.rainbow(np.linspace(0, 1, len(vbeam_))))

    # iterate vbeam
    for vbeam in vbeam_:
        # Generate positions, velocities, run sim
        pos, vel = twostream(10000, L, vbeam)
        s = Summary()
        pos, vel = run(pos, vel, L, ncells, 
                    out=[s], #[p, s],                      # These are called each output
                    output_times=linspace(0.,20,50)) # The times to output    
        
        # Plot firstharmonic damping
        plt.plot(s.t, s.firstharmonic)
        plt.xlabel("Time [Normalised]")
        plt.ylabel("First harmonic amplitude [Normalised]")
        plt.yscale('log')
    

    

    plt.show()

    # To see end status of particles
    plt.figure()
    plt.plot(pos, vel, marker='o' ,linestyle='none')
    plt.xlabel("Position")
    plt.ylabel("Velocity")
    #plt.yscale('log')
    
    plt.show()