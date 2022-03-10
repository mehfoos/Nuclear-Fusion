# Code to import and run epc1d functions

from epc1d import *
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np
from matplotlib.pyplot import cm

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

def task1_noiseAmplitude(noiseLevelThresholdR, L, ncells, npart):
    plt.show()
    samples = len(noiseLevelThresholdR)
    plt.boxplot(noiseLevelThresholdR)
    plt.ylabel("Noise Ampitude Threshold")
    plt.title("L={0:.2f}; ncells={1}; npart={2}; samples={3}".format(L,ncells,npart,samples))
    plt.suptitle("Landau Damping Noise Amplitude Threshold: Boxplot")
    print("Mean noiseAmplitudes: ", np.mean(noiseLevelThresholdR))
    print("SD Sample noiseAmplitudes: ", np.std(noiseLevelThresholdR, ddof=1))
    print("SE Sample noiseAmplitudes: ", np.std(noiseLevelThresholdR, ddof=1)/np.sqrt(samples))
    plt.show()

def task1_noiseAmplitudeCells(noiseLevelThresholdR, L, ncellsR, npart, samples):

    # Preprocess, initialise
    ncellsUniqueSorted = np.unique(ncellsR)
    means = np.zeros(len(ncellsUniqueSorted))
    stdErrs = np.zeros(len(ncellsUniqueSorted))
    stdDevsSample = np.zeros(len(ncellsUniqueSorted))

    # Calculate standard errors and means
    for iNcell in range(len(ncellsUniqueSorted)):
        means[iNcell] = np.mean(noiseLevelThresholdR[np.where(ncellsR == ncellsUniqueSorted[iNcell])])
        stdDevsSample[iNcell] = np.std(noiseLevelThresholdR[np.where(ncellsR == ncellsUniqueSorted[iNcell])], ddof=1)
    stdErrs = stdDevsSample/np.sqrt(samples)

    # Polynomial fit
    ncellsNoiseFitCoeffs = np.polyfit(ncellsUniqueSorted, means, 2)
    ncellsNoiseFitF = np.poly1d(ncellsNoiseFitCoeffs)
    fitValues = ncellsNoiseFitF(ncellsUniqueSorted)
    rsquared = r_squared(means,fitValues)

    # Plot
    plt.close("all")
    plt.ylabel("Noise Ampitude Threshold")
    plt.xlabel("Number of Cells")
    plt.errorbar(ncellsUniqueSorted,means,yerr=stdErrs,linestyle='none',uplims=True, lolims=True,label="Noise Amplitude Error Range")
    plt.plot(ncellsUniqueSorted, ncellsNoiseFitF(ncellsUniqueSorted), linestyle='dashed',label=r"Polyfit: " + '{}'.format(ncellsNoiseFitF) + "\n $R^{2}$=" + f"{rsquared:.1f}")
    plt.title("L={0:.2f}; ncells={1}; npart={2}; samples={3}".format(L,ncells,npart,samples))
    plt.suptitle("Landau Damping Noise Amplitude Threshold vs Cells")
    plt.legend()
    plt.show()
    
    print("ncellsUniqueSorted: ", ncellsUniqueSorted)
    print("means: ", means)
    print("stdErrs: ", stdErrs)

def task1_noiseAmplitudeParticles(noiseLevelThresholdR, L, ncells, npartR, samples):

    # Preprocess, initialise
    nParticlesUniqueSorted = np.unique(npartR)
    means = np.zeros(len(nParticlesUniqueSorted))
    stdErrs = np.zeros(len(nParticlesUniqueSorted))
    stdDevsSample = np.zeros(len(nParticlesUniqueSorted))

    # Calculate standard errors and means
    for iNParticle in range(len(nParticlesUniqueSorted)):
        means[iNParticle] = np.mean(noiseLevelThresholdR[np.where(npartR == nParticlesUniqueSorted[iNParticle])])
        stdDevsSample[iNParticle] = np.std(noiseLevelThresholdR[np.where(npartR == nParticlesUniqueSorted[iNParticle])], ddof=1)
    stdErrs = stdDevsSample/np.sqrt(samples)

    # Polynomial fit
    nParticlesNoiseFitCoeffs = np.polyfit(nParticlesUniqueSorted, means, 2)
    nParticlesNoiseFitF = np.poly1d(nParticlesNoiseFitCoeffs)
    fitValues = nParticlesNoiseFitF(nParticlesUniqueSorted)
    rsquared = r_squared(means,fitValues)

    # Plot
    plt.close("all")
    plt.ylabel("Noise Ampitude Threshold")
    plt.xlabel("Number of Particles")
    plt.errorbar(nParticlesUniqueSorted,means,yerr=stdErrs,linestyle='none',uplims=True, lolims=True,label="Noise Amplitude Error Range")
    plt.plot(nParticlesUniqueSorted, nParticlesNoiseFitF(nParticlesUniqueSorted), linestyle='dashed',label=r"Polyfit: " + '{}'.format(nParticlesNoiseFitF) + "\n $R^{2}$=" + f"{rsquared:.1f}")
    plt.title("L={0:.2f}; ncells={1}; npart={2}; samples={3}".format(L,ncells,npart,samples))
    plt.suptitle("Landau Damping Noise Amplitude Threshold vs Particles")
    plt.legend()
    plt.show()
    
    print("ncellsUniqueSorted: ", nParticlesUniqueSorted)
    print("means: ", means)
    print("stdErrs: ", stdErrs)



if __name__ == "__main__":
    
    PLOTS = False

    # Configuration (All configurable to include multiple values to loop through)
    samples = 4
    L_ = [4.*pi]
    npart_ = [500, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000, 50000]
    ncells_ = [8]

    # Initialisations
    i0 = -1
    i1 = -1 # iterators
    i2 = -1
    i3 = -1
    npartR = np.array([]) # results arrays
    rsiR = np.array([])
    ncellsR = np.array([])
    LR = np.array([])
    noiseLevelThresholdR = np.array([])

    # Unique colors
    color = iter(cm.rainbow(np.linspace(0, 1, samples)))

    for rseed in range(samples): # Loop through for repeats; stats
        i0 += 1
        c = next(color) # Plot color iterate
        random.seed(i0) # random seed fix for when needed to compare if optimisation makes a difference from original source

        for L in L_: # Loop through different values of L
            i1 += 1 

            for npart in npart_: # Loop through values of nparts
                i2 += 1
                
                for ncells in ncells_:  # Loop through values of ncells
                    i3 += 1

                    # Run the main sim
                    pos, vel = landau(npart, L)
                    s = Summary()
                    run(pos, vel, L, ncells, [s], linspace(0.,20,50))

                    # Convert to ndarray for convenience
                    s.t = array(s.t)
                    s.firstharmonic = array(s.firstharmonic)

                    plt.plot(s.t, s.firstharmonic, label=r"First harmonic amplitude [Normalised]", c = c)
                    plt.yscale('log')
                    plt.xlabel("Time [Normalised]")
                    plt.ylabel("First harmonic amplitude [Normalised]")
                    plt.yscale('log')
                    plt.title("L={0:.2f}; ncells={1}; npart={2}".format(L,ncells,npart))
                    plt.suptitle("Landau Damping First Harmonic")

                    # Find peaks
                    peaksI, peakProperties = find_peaks(s.firstharmonic) # indexes of peaks
                    peaks = s.firstharmonic[peaksI]

                    # Plot peaks
                    plt.plot(s.t[peaksI], s.firstharmonic[peaksI], marker='o', linestyle='none',label=r"Peaks", c = c)

                    # Find the last good peak / last bad peak
                    peaksDiff = np.ediff1d(s.firstharmonic[peaksI]) # diff between consecutive values
                    risePeaksIndices, = np.where(peaksDiff>0) # offset by one, since it's diff
                    if len(risePeaksIndices) ==0: # if no negative peak difference
                        lastGoodPeakI = len(peaksI) - 1 # Set to last peak
                        firstBadPeak = lastGoodPeakI # Assuming it gets worse after this
                    else:
                        lastGoodPeakI = risePeaksIndices[0] # Set to last non negative peak
                        firstBadPeak = lastGoodPeakI+1
                    #print("lastGoodPeakI", lastGoodPeakI)

                    # Line fit for good peaks
                    fitLineGoodPeaks = np.polyfit(s.t[peaksI[:lastGoodPeakI+1]], np.log(s.firstharmonic[peaksI[:lastGoodPeakI+1]]), 1) # taking log of peaks, and doing a line fit
                    fitLineGoodPeaksF = np.poly1d(fitLineGoodPeaks) # Create a fit function from the polyfit
                    fit_Rsq = r_squared(s.firstharmonic[peaksI[:lastGoodPeakI+1]], np.exp(fitLineGoodPeaksF(s.t[peaksI[:lastGoodPeakI+1]]))) # Find R^2
                    
                    # Results
                    noiseLevelThreshold = s.firstharmonic[peaksI[firstBadPeak]] # 3-d array format
                    #print("Noise Level Amplitude Threshold for ", L, "L; ", npart, "nparts; ", ncells, "ncells is: ", noiseLevelThreshold)
                    
                    # Results in 1d arrays
                    rsiR = np.append(rsiR, i0)
                    npartR = np.append(npartR, npart)
                    ncellsR = np.append(ncellsR, ncells)
                    LR = np.append(LR, L)
                    noiseLevelThresholdR = np.append(noiseLevelThresholdR, noiseLevelThreshold)

                    #plt.close("all")
                    plt.plot(s.t[peaksI[:lastGoodPeakI+1]],np.exp(fitLineGoodPeaksF(s.t[peaksI[:lastGoodPeakI+1]])), linestyle='dashed', label=r"Fit Line for low noise data", c = c)


                    # Show plot
                    if PLOTS:
                        plt.legend()
                        plt.show()

    
    print(noiseLevelThreshold)

    # Function calls for individual tasks (Will be commented out where not applicable)
    # task1_noiseAmplitude(noiseLevelThresholdR, L, ncells, npart) 
    # task1_noiseAmplitudeCells(noiseLevelThresholdR, L, ncellsR, npart, samples)
    task1_noiseAmplitudeParticles(noiseLevelThresholdR, L, ncells, npartR, samples)

    
    
    # plot 
    #plt.plot(s.t, s.firstharmonic, label=r"First harmonic amplitude [Normalised]")
    plt.xlabel("Number of Particles")
    plt.ylabel("Number of Cells")
    plt.title("Landau Damping Amplitude Noise Threshold vs Particle and Cell Number")
    plt.scatter(npartR, npartR, s=noiseLevelThresholdR)




