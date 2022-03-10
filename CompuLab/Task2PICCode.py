#!/usr/bin/env python
#
# Electrostatic PIC code in a 1D cyclic domain

import epc1d as originalPic2
import Task1PICcode as FirstTask

from numpy import linspace
from numpy import array, pi, std, exp, diff, sqrt

import numpy as np

from scipy import stats 


from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt # Matplotlib plotting library
from symbol import nonlocal_stmt

def func(x, *params):
    a=params[0]
    b=params[1]
   # c=params[2]
    return a* exp(b*x) 

pguess = [0.1,-0.02]


def OverallRun():
    
    global popt
    if False:
        # 2-stream instability
        L = 100
        ncells = 20
        pos, vel = originalPic2.twostream(10000, L, 3.)
    else:
        # Landau damping
        L = 4.*pi
        ncells = 20
        npart = 1000
        pos, vel = originalPic2.landau(npart, L)
    
    # Create some output classes
    p = originalPic2.Plot(pos, vel, ncells, L) # This displays an animated figure
    s = originalPic2.Summary()                 # Calculates, stores and prints summary info
    
    # Run the simulation
    pos, vel = originalPic2.run(pos, vel, L, ncells, 
                out=[
                    #p,
                     s],                      # These are called each output
                output_times=linspace(0.,20,500)) # The times to output
    
    # Summary stores an array of the first-harmonic amplitude
    # Make a semilog plot to see exponential damping
    Firstharm,FirstHarmonictime,peaktops= FirstTask.PlotPeaks(s.firstharmonic, s.t)
    plt.figure()
    plt.plot(s.t, s.firstharmonic)
    plt.plot(FirstHarmonictime[peaktops],Firstharm[peaktops],"x")
    print("This is the first harmohnic time spacing", FirstHarmonictime[peaktops])
    peakno, backslice, forwardslice = FirstTask.findNoise(Firstharm,FirstHarmonictime, peaktops)
    print("This is peakno, backslice, forwardslice,", peakno, backslice,forwardslice)
   # plt.plot(FirstHarmonictime[peaktops[:forwardslice]],Firstharm[peaktops[:forwardslice]])
    print("This is the first harmohnic time spacing", FirstHarmonictime[peaktops[:forwardslice]])
    diffArray = diff(FirstHarmonictime[peaktops[:forwardslice]])
    print("This is the difference Array", diffArray)
    try:
        popt, pcov = curve_fit(func,FirstHarmonictime[peaktops[:forwardslice]],Firstharm[peaktops[:forwardslice]], pguess, maxfev=2000)
    except (RuntimeError, TypeError):
       print(popt)
    fitcurve= func(FirstHarmonictime[peaktops[:forwardslice]],*popt)
    #print(fitcurve)
    plt.plot(FirstHarmonictime[peaktops[:forwardslice]],fitcurve[:forwardslice])
    #print(fitcurve[:forwardslice])
    plt.suptitle("First harmonic with 1000 particles and 20 cells")
    plt.xlabel("Time [Normalised]")
    plt.ylabel("First harmonic amplitude [Normalised]")
    plt.yscale('log')
    plt.ioff() # This so that the windows stay open
    plt.show()

    return diffArray



####################################################################

if __name__ == "__main__":
    # Generate initial condition
    listofTimeDifferences = []
    #
    for r in range(0,20):
        TimeDifference = OverallRun()
        print("This is TimeSplits", TimeDifference)
        listofTimeDifferences.extend(TimeDifference)
        
print("This is the listofTimeDifferences",  listofTimeDifferences)
p= stats.describe(listofTimeDifferences)
print(p)

numberofTimeDifferences = p[0]
print(numberofTimeDifferences)
meanvalue = p[2]
print("This is the meanvalue", meanvalue)
variance = p[3]
print("This is the variance", variance)
StandardError = sqrt(variance)/sqrt(numberofTimeDifferences)
print("This is the standard Error", StandardError)
AngularFrequency = 2*pi/(2*meanvalue)
print("This is the angular frequency", AngularFrequency)


        
        
    
    
   
    
    
  
