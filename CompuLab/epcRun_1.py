# Code to import and run epc1d functions

#from matplotlib.lines import _LineStyle
from epc1d import *
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np

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

if __name__ == "__main__":
    random.seed(10)
    PLOTS = True

    # Configuration
    L_ = [4.*pi]
    npart_ = [1000]
    
    ncells_ = [8]

    # Results
    i1 = -1
    i2 = -1
    i3 = -1
    npartR = np.array([])
    ncellsR = np.array([])
    LR = np.array([])
    noiseLevelThresholdR = np.array([])

    for L in L_:
        i1 += 1
        for npart in npart_:
            i2 += 1
            for ncells in ncells_:    
                i3 += 1
            
                # Run the main sim

                pos, vel = landau(npart, L)
                s = Summary()

                run(pos, vel, L, ncells, [s], linspace(0.,20,50))

                # convert to ndarray for convenience
                s.t = array(s.t)
                s.firstharmonic = array(s.firstharmonic)

                plt.plot(s.t, s.firstharmonic, label=r"First harmonic amplitude [Normalised]")
                plt.yscale('log')
                plt.xlabel("Time [Normalised]")
                plt.ylabel("First harmonic amplitude [Normalised]")
                plt.yscale('log')
                plt.title("L={0:.2f}; ncells={1}; npart={2}".format(L,ncells,npart))
                plt.suptitle("Landau Damping First Harmonic")

                # Find peaks
                peaks, peakProperties = find_peaks(s.firstharmonic) # indexes of peaks

                # Plot peaks
                plt.plot(s.t[peaks], s.firstharmonic[peaks], marker='o', linestyle='none',label=r"Peaks")

                # Try fitting line to first three peaks
                fitLine3Pk = np.polyfit(s.t[peaks[:3]], np.log(s.firstharmonic[peaks[:3]]), 1) # taking log of peaks, and doing a line fit
                fitLine3PkF = np.poly1d(fitLine3Pk) # Create a fit function from the polyfit
                fit_Rsq = r_squared(s.firstharmonic[peaks[:3]], np.exp(fitLine3PkF(s.t[peaks[:3]]))) # Find R^2
                

                # 
                #ErrFit3 = abs(s.firstharmonic[peaks[:3]] - np.exp(fitLine3PkF(s.t[peaks[:3]])))
                #errFit = np.polyfit(s.firstharmonic[peaks[:3]], np.log(abs(fitLine3PkF(s.t[peaks[:3]]))), 1) #taking log of first 3 err peaks
                #errFitF = np.poly1d(fitLine3Pk)
                #plt.show()
                #errorExpPeaks = s.firstharmonic[peaks] - np.exp(fitLine3PkF(s.t[peaks])) # max error in fit
                #badPeaks = errorExpPeaks > np.exp(errFitF(s.t[peaks]))
                
                # Find Bad Peaks
                acceptanceThreshold = 1.5
                badPeaks = acceptanceThreshold*np.exp(fitLine3PkF(s.t[peaks])) < s.firstharmonic[peaks] # If the first harmonic is less than fit by 50percent, then it's probably where noise is too high; this is boolean array
                badPeaksI  = np.flatnonzero(badPeaks) # find the indices of values where we think the noise is too high
                if len(badPeaksI) ==0: # If noise isn't high for any peak:
                    lastGoodPeakI = len(peaks)-1 # Set to last peak
                else:
                    lastGoodPeakI = badPeaksI[0]-1
                if fit_Rsq < 0.5:
                    if s.firstharmonic[peaks[1]] > acceptanceThreshold*s.firstharmonic[peaks[0]]: # If the second value is higher than the first, then we can consider first peak to be threshold
                        lastGoodPeakI = 0
                    else: # Consider second value to be threshold
                        lastGoodPeakI = 1
                if s.firstharmonic[peaks[1]] > acceptanceThreshold*s.firstharmonic[peaks[0]]: # If the second value is higher than the first, then we can consider first peak to be threshold
                    lastGoodPeakI = 0

                # Debugs
                print("3pk Fit Rsq: ", fit_Rsq)
                print("badPeaks ", badPeaks)
                print("Fit Amplitudes*coeff", 1.5*np.exp(fitLine3PkF(s.t[peaks])))
                print("s.firstharmonic[peaks]", s.firstharmonic[peaks])
                print("badPeaksI", badPeaksI)
                print("lastGoodPeakI", lastGoodPeakI)
                
                # Results
                noiseLevelThreshold = s.firstharmonic[peaks[lastGoodPeakI]] # 3-d array format
                print("Noise Level Amplitude Threshold for ", L, "L; ", npart, "nparts; ", ncells, "ncells is: ", noiseLevelThreshold)
                
                # Results in 1d arrays
                npartR = np.append(npartR, npart)
                ncellsR = np.append(ncellsR, ncells)
                LR = np.append(LR, L)
                noiseLevelThresholdR = np.append(noiseLevelThresholdR, noiseLevelThreshold)


                #plt.close("all")
                plt.plot(s.t[peaks[:lastGoodPeakI+1]],np.exp(fitLine3PkF(s.t[peaks[:lastGoodPeakI+1]])), linestyle='dashed', label=r"Fit Line for low noise data")


                # Show plot
                plt.legend()
                if PLOTS:
                    plt.show()

        
    print(noiseLevelThreshold)

    # plot 
    #plt.plot(s.t, s.firstharmonic, label=r"First harmonic amplitude [Normalised]")
    plt.xlabel("Number of Particles")
    plt.ylabel("Number of Cells")
    plt.title("Landau Damping Amplitude Noise Threshold vs Particle and Cell Number")
    plt.scatter(npartR, ncellsR, s=noiseLevelThresholdR)




