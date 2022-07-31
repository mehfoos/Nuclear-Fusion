# Code to import and run epc1d functions

from sys import stderr
from epc1d_Opt2 import *
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np
from matplotlib.pyplot import cm
from timeit import default_timer as timer
from scipy.stats import linregress


def r_squared(data, fit):
	""" Calculate the R^2 value of a fit, f, to a data set, x """
	
	if len(data) != len(fit):
		print("Data and fit lists must be of equal length")
		return None

	m = np.sum(data)/len(data)
	
	# Sum of squares between data and mean
	ss = np.sum((data-m)**2)
	
	# Sum of residuals between fit and data 
	res = np.sum((data-fit)**2)
		
	# This is the definition of R^2
	return 	1 - res/ss

def plotHelper(textChanging,textPlot,x,y,stdErrs,L,ncells,npart,samples):
    plt.close("all")
    fitCoeffs =  np.polyfit(x, y, 1)
    fitF = np.poly1d(fitCoeffs)
    fitValues = fitF(x)
    rsquared = r_squared(y,fitValues)
    plt.errorbar(x,y,yerr=stdErrs,linestyle='none',uplims=True, lolims=True,label=textPlot + " Error Range")
    plt.plot(x, fitValues, linestyle='dashed',label=r"Polyfit: " + '{}'.format(fitF) + "\n $R^{2}$=" + f"{rsquared:.1f}")
    # plt.title("L={0:.2f}; ncells={1}; simulations={2};".format(L,ncells,samples)) # Whatever's chaning #L={0:.2f};
    plt.title("L={0:.2f}; npart={1}; ncells={2}; simulations={3};".format(L,npart,ncells,samples)) # Whatever's chaning #L={0:.2f};
    plt.suptitle("Two Streams: " + textPlot + " vs " + textChanging)
    plt.xlabel(textChanging) #
    plt.ylabel(textPlot + " Means")
    #plt.yscale('log')
    plt.legend()
    plt.show()
    print(textChanging, ": ", x)
    print(textPlot, " Means: ", y)
    print(textPlot, " Standard Errors: ", stdErrs)

if __name__ == "__main__":
    # Configuration
    # 2-stream instability

    random.seed(1)

    L = 100
    ncells = 20
    vbeam_ = np.array(range(1,20,1)).astype(float)
    # vbeam_ = [17., 3.]
    sims = 5
    npart = 10000
    
    # Unique colors
    # color = iter(cm.Greys(np.linspace(0, 1, len(vbeam_))))
    # color = iter(cm.rainbow(np.linspace(0, 1, sims)))

    # Initialise inefficient lists for quick coding
   

    oscFrequencies = zeros(len(vbeam_))
    oscFrequenciesSe = zeros(len(vbeam_))
    slopes = zeros(len(vbeam_))
    slopesSe = zeros(len(vbeam_))

    # iterate vbeam
    i0=-1
    for vbeam in vbeam_:

        # Per Beam Items
        oscFrequenciesPerBeam = np.array([])
        oscSlopesPerBeam = np.array([])

        # Iterators, initialise, etc.
        i0 += 1
        color = iter(cm.rainbow(np.linspace(0, 1, sims)))

        # iterate sims
        for sim in range(sims):
            # Generate positions, velocities, run sim
            pos, vel = twostream(npart, L, vbeam)
            s = Summary()
            pos, vel = run(pos, vel, L, ncells, 
                        out=[s], #[p, s],                      # These are called each output
                        output_times=linspace(0.,20,50)) # The times to output    
            
            # Iterate color
            c = next(color)

            # Convert to ndarray for convenience
            s.t = array(s.t)
            s.firstharmonic = array(s.firstharmonic)

            # Find peaks
            peaksI, peakProperties = find_peaks(s.firstharmonic) # indexes of peaks
            #print(peaksI)
            peaks = s.firstharmonic[peaksI]
            peakTimes = s.t[peaksI]
            
            # Frequencies
            peakTimeDiffs = np.ediff1d(s.t[peaksI]) # Get time differences
            oscFrequenciesPerBeam = np.append(oscFrequenciesPerBeam, np.pi/peakTimeDiffs)

            # Slopes
            linreg = linregress(peakTimes, np.log(peaks))
            oscSlopesPerBeam = np.append(oscSlopesPerBeam, linreg.slope)
            
            # Debug plot
            plt.plot(s.t, s.firstharmonic, label=r"First harmonic amplitude [Normalised]", c = c)
            plt.yscale('log')
            plt.xlabel("Time [Normalised]")
            plt.ylabel("First harmonic amplitude [Normalised]")
            plt.plot(s.t[peaksI], s.firstharmonic[peaksI], marker='o', linestyle='none',label=r"Peaks", c = c)
            plt.plot(peakTimes, np.exp(linreg.intercept + linreg.slope*peakTimes), label='Fitted line')
            #plt.show()
            plt.close("all")

        oscFrequencies[i0] = np.mean(oscFrequenciesPerBeam)
        oscFrequenciesSe[i0] = np.std(oscFrequenciesPerBeam, ddof=1) / np.sqrt(np.size(oscFrequenciesPerBeam))
        slopes[i0] = np.mean(oscSlopesPerBeam)
        slopesSe[i0] = np.std(oscSlopesPerBeam, ddof=1) / np.sqrt(np.size(oscSlopesPerBeam))

        

        # # Plot firstharmonic damping
        # if vbeam%2==0:
        #     plt.plot(s.t, s.firstharmonic, c=c, label = "vBeam = " + str(vbeam))
        # else:
        #     plt.plot(s.t, s.firstharmonic, c=c, linestyle="dashed", label = "vBeam = " + str(vbeam))
            
        
    # # To plot All vbeam# comment out when not needed
    # plt.yscale('log')
    # plt.xlabel("Time [Normalised]")
    # plt.ylabel("First harmonic amplitude [Normalised]")
    # plt.title("Manual approach to spotting instability")
    # plt.legend()
    # plt.show()
    # plt.close("all")

    # # To see end status of particles
    # plt.figure()
    # plt.plot(pos, vel, marker='o' ,linestyle='none', label = "vBeam = " + str(vbeam))
    # plt.xlabel("Position")
    # plt.ylabel("Velocity")
    # #plt.yscale('log')
    
    # plt.show()

    textChanging = "vBeam"
    textPlot = "Angular Freq"
    x = vbeam_
    y = oscFrequencies
    stdErrs = oscFrequenciesSe
    plotHelper(textChanging,textPlot,x,y,stdErrs,L,ncells,npart,sims)

    textChanging = "vBeam"
    textPlot = "Growth Rate"
    x = vbeam_
    y = slopes
    stdErrs = slopesSe
    plotHelper(textChanging,textPlot,x,y,stdErrs,L,ncells,npart,sims)



    