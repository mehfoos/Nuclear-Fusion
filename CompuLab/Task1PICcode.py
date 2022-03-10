#!/usr/bin/env python
#
# Electrostatic PIC code in a 1D cyclic domain

import epc1d as originalPic

from numpy import linspace
from numpy import array, pi, std

import numpy as np

from scipy.signal import find_peaks
import matplotlib.pyplot as plt # Matplotlib plotting library

def findNoise(datapeaks, timepeaks, peaks):
    listmine = []
    for pr in range(len(datapeaks[peaks])):
        x = datapeaks[peaks[pr]]
        listmine.append(x)
        
    listmineNew = []
    print("This is listmine", listmine)
    previous=100
    for g in range(len(listmine)):
        if listmine[g]> previous:
            break
        else:
            previous=listmine[g]
            listmineNew.append(previous)
            
    slen= len(listmineNew)
    m=len(listmine)
    p=m-slen
    
    return m, p, slen
    
def NoiseAverage(Noise, peaknumber):
    NoiseAveraged = sum(Noise)/peaknumber
    NoiseSTD = std(Noise)
    return NoiseSTD, NoiseAveraged

def check(datapeaks, timepeaks):  
    
     for p in range(len(datapeaks)):
        x = datapeaks[p]
        print("This is x", x)
      
def PlotPeaks(firstharmonic, time):  
    data = array(firstharmonic)
    timepeak = array(time)
    peaks, _ = find_peaks(data)
    
    return data, timepeak, peaks

def OverallRun():
    if False:
        # 2-stream instability
        L = 100
        ncells = 20
        pos, vel = originalPic.twostream(10000, L, 3.)
    else:
        # Landau damping
        L = 4.*pi
        ncells = 20
        npart = 1000
        pos, vel = originalPic.landau(npart, L)
    
    # Create some output classes
    p = originalPic.Plot(pos, vel, ncells, L) # This displays an animated figure
    s = originalPic.Summary()                 # Calculates, stores and prints summary info
    
    # Run the simulation
    pos, vel = originalPic.run(pos, vel, L, ncells, 
                out=[
                    #p,
                     s],                      # These are called each output
                output_times=linspace(0.,20,500)) # The times to output
    
    # Summary stores an array of the first-harmonic amplitude
    # Make a semilog plot to see exponential damping
    plt.figure()
    plt.plot(s.t, s.firstharmonic)
    y,x,peaktops= PlotPeaks(s.firstharmonic, s.t)
    plt.plot(x[peaktops],y[peaktops],"x")
    #plt.plot(timepeak[peaks],data[peaks])
    check(y,x)
    peakno, backslice, forwardslice = findNoise(y,x, peaktops)
    plt.plot(x[peaktops[:forwardslice]],y[peaktops[:forwardslice]])
    StandardDeviation, NoiseLevel = NoiseAverage(y[peaktops[-backslice:]], peakno)
    plt.suptitle("First harmonic with 1000 particles and 20 cells")
    plt.xlabel("Time [Normalised]")
    plt.ylabel("First harmonic amplitude [Normalised]")
    plt.yscale('log')
    plt.ioff() # This so that the windows stay open
    plt.show()
    
    return StandardDeviation, NoiseLevel



####################################################################

if __name__ == "__main__":
    # Generate initial condition
    #
    StdList =[]
    NoiseList= []
    for x in range(0,3):
        stdfh, Noisefh = OverallRun()
        StdList. append(stdfh)
        NoiseList.append(Noisefh)
    StdListAverage = sum(StdList)/3
    NoiseListAverage =sum(NoiseList)/3

    print("This is the StdList", StdList)
    print("This is the NoiseList", NoiseList)
    print("This is the StdListAverage", StdListAverage)
    print("This is the NoiseListAverage", NoiseListAverage)

    
    