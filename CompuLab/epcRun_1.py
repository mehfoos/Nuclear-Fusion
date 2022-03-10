# Code to import and run epc1d functions

#from matplotlib.lines import _LineStyle
from epc1d import *
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np


if __name__ == "__main__":
    # Configuration
    L = 4.*pi
    pos, vel = landau(1000, L)
    s = Summary()
    ncells = 8

    # Run the main sim
    run(pos, vel, L, ncells, [s], linspace(0.,20,50))

    # Find the peaks
    peaks, _ = find_peaks(s.firstharmonic)
    #plt.plot(ArK_Peaks_Offset,ArK_CrosssectionData[ArK_Peaks],marker='o', linestyle='none', label=r"Peaks")
    

    # Plot configurations
    plt.plot(s.t, s.firstharmonic, label=r"First Harmonic Normalised Amplitude")
    #plt.plot(s.t[peaks],s.firstharmonic[peaks],marker='o', lineStyle='none', label=r"Peaks")
    if isinstance(peaks, np.ndarray):
        print("This is a numpy array: peaks")
    if isinstance(s.t, np.ndarray):
        print("This is a numpy array: s.t")
    if isinstance(s.firstharmonic, np.ndarray):
        print("This is a numpy array: s.firstharmonic")
    #print("Peaks indices: ", peaks)
    #print("s.t: ", s.t)
    #print("s.firstharmonic: ", s.firstharmonic)
    plt.xlabel("Time [Normalised]")
    plt.ylabel("First harmonic amplitude [Normalised]")
    plt.yscale('log')
    plt.show()



