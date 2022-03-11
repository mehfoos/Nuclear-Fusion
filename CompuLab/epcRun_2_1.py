# Code to import and run epc1d functions

from sys import stderr
from epc1d_Opt2 import *
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np
from matplotlib.pyplot import cm
from timeit import default_timer as timer
from scipy.stats import linregress

FUNCPLOTS = False

def plotHelper(textChanging,textPlot,x,y,stdErrs,L,ncells,npart,samples):
    fitCoeffs =  np.polyfit(x, y, 1)
    fitF = np.poly1d(fitCoeffs)
    fitValues = fitF(x)
    rsquared = r_squared(y,fitValues)
    plt.errorbar(x,y,yerr=stdErrs,linestyle='none',uplims=True, lolims=True,label=textPlot + " Error Range")
    plt.plot(x, fitValues, linestyle='dashed',label=r"Polyfit: " + '{}'.format(fitF) + "\n $R^{2}$=" + f"{rsquared:.1f}")
    # plt.title("L={0:.2f}; ncells={1}; simulations={2};".format(L,ncells,samples)) # Whatever's chaning #L={0:.2f};
    plt.title("L={0:.2f}; npart={1}; simulations={2};".format(L,npart,samples)) # Whatever's chaning #L={0:.2f};
    plt.suptitle("Landau Damping: " + textPlot + " vs " + textChanging)
    plt.xlabel(textChanging) #
    plt.ylabel(textPlot + " Means")
    #plt.yscale('log')
    plt.legend()
    plt.show()
    print(textChanging, ": ", x)
    print(textPlot, " Means: ", y)
    print(textPlot, " Standard Errors: ", stdErrs)

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

def task2_peakFrequency(L, ncells, npart, samples, oscFrequency):
    mean = np.mean(oscFrequency)
    stdDevSample = np.std(oscFrequency, ddof=1)
    stdErrSample = stdDevSample/np.sqrt(len(oscFrequency))

    plt.close("all")
    plt.boxplot(oscFrequency)
    plt.errorbar(1,mean,yerr=stdErrSample,marker='o',linestyle='none',uplims=True, lolims=True,label="Means and Standard Error")
    plt.ylabel("Frequncy (Normalised)")
    plt.title("L={0:.2f}; ncells={1}; npart={2}; simulations={3};".format(L,ncells,npart,samples))
    plt.suptitle("Landau Damping Frequency: Boxplot")

    if FUNCPLOTS:
        print("Mean Frequency: ", np.mean(oscFrequency))
        print("SD Sample Frequency: ", np.std(oscFrequency, ddof=1))
        print("SE Frequency: ", stdErrSample)

    plt.legend()
    if FUNCPLOTS:
        plt.show()
    else:
        plt.close("all")

    oscFreqMean = mean
    oscFreqSe = stdErrSample

    return oscFreqMean, oscFreqSe 

def task2_globalPeakFit(L, ncells, npart, samples, goodPeakAmplitudes, goodPeakTimes, oscFreqMean, oscFreqSe):
    result = linregress(goodPeakTimes,np.log(goodPeakAmplitudes))
    print(result)
    print("On Linear scale, Std Error of slope = ", np.exp(result.stderr))
    print("On Linear scale, Std Error of intercept = ", np.exp(result.intercept_stderr))
    
    plt.close("all")
    plt.plot(goodPeakTimes,goodPeakAmplitudes,'o',label="Low Noise First Harmonic Peaks")
    plt.plot(goodPeakTimes, np.exp(result.intercept + result.slope*goodPeakTimes), label='Fitted line')
    plt.title("L={0:.2f}; npart={1}; ncells={2}; simulations={3}".format(L,npart,ncells,samples))
    plt.suptitle("Landau Damping Dampening rate curve fit")
    plt.xlabel("Time [s]")
    plt.ylabel("First Harmonic Amplitude peaks")
    plt.yscale('log')
    plt.legend()
    if FUNCPLOTS:
        plt.show()
    else:
        plt.close("all")

    # Damping ratio, angular freq
    #dampingRatio = result.slope/oscFreqMean
    dampingRatio = result.slope/sqrt(result.slope**2 + oscFreqMean**2)
    dampingRatioSe = sqrt(oscFreqSe**2 + result.stderr**2)
    decayRate = result.slope
    decayRateSe = result.stderr
    if FUNCPLOTS:
        print("angular freq = ", oscFreqMean, "+-", oscFreqSe)
        print("damping ratio = ", dampingRatio, "+-", dampingRatioSe)
        print("decay rate: ", result.slope, "+-", result.stderr)
    return decayRate, decayRateSe


if __name__ == "__main__":

    # Configuration (All configurable to include multiple values to loop through)
    samples = 10
    L_ = [4.*pi]
    #npart_ = [500, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000, 50000]
    npart_ = [1000]
    # ncells_ = [20]
    ncells_ = [10, 50, 100, 150, 200]

    # Initialisations
    i0 = -1
    i1 = -1 # iterators
    i2 = -1
    i3 = -1
    #random.seed(1)

    # Initialise for all that needs resetting at main level
    oscFreqMeanM = np.array([])
    oscFreqSeM = np.array([]) 
    decayRateMeanM = np.array([])
    decayRateSeM = np.array([]) 

    for L in L_: # Loop through different values of L
        i0 += 1 

        

        for npart in npart_: # Loop through values of nparts
            i2 += 1
            
            for ncells in ncells_:  # Loop through values of ncells
                i3 += 1

                # Emptied initialised every time there's a new set of samples to be run
                oscFrequency = np.array([]) 
                goodPeakAmplitudes = np.array([])
                goodPeakTimes = np.array([])


                for sample in range(samples): # loop through multiple simulations of same condition
                    
                    # Run the main sim
                    timeStart = timer()
                    pos, vel = landau(npart, L)
                    s = Summary()
                    run(pos, vel, L, ncells, [s], linspace(0.,20,50))
                    timeStop = timer()

                    # Convert to ndarray for convenience
                    s.t = array(s.t)
                    s.firstharmonic = array(s.firstharmonic)

                    # Find peaks
                    peaksI, peakProperties = find_peaks(s.firstharmonic) # indexes of peaks
                    peaks = s.firstharmonic[peaksI]

                    # Find the last good peak / last bad peak
                    peaksDiff = np.ediff1d(s.firstharmonic[peaksI]) # diff between consecutive values
                    risePeaksIndices, = np.where(peaksDiff>0) # offset by one, since it's diff
                    if len(risePeaksIndices) ==0: # if no negative peak difference
                        lastGoodPeakI = len(peaksI) - 1 # Set to last peak
                        firstBadPeak = lastGoodPeakI # Assuming it gets worse after this
                    else:
                        lastGoodPeakI = risePeaksIndices[0] # Set to last non negative peak
                        firstBadPeak = lastGoodPeakI+1

                    # Frequency estimate, # can comment out when not needed
                    peakTimeDiffs = np.ediff1d(s.t[peaksI])
                    if lastGoodPeakI > 0: # At least two valid peaks to calc frequency
                        oscFrequency = np.append(oscFrequency, np.pi/peakTimeDiffs[0:lastGoodPeakI])
                    
                    # Peaks collation for fit outside of loop; # can comment out when not needed
                    if lastGoodPeakI > 0: # At least two valid peaks to calc frequency
                        goodPeakAmplitudes = np.append(goodPeakAmplitudes, s.firstharmonic[peaksI[0:lastGoodPeakI+1]])
                        goodPeakTimes = np.append(goodPeakTimes, s.t[peaksI[0:lastGoodPeakI+1]])
                    
                    # end of samples loop
                
                # Get Angular Frequency values and decayRate
                oscFreqMean, oscFreqSe = task2_peakFrequency(L, ncells, npart, samples, oscFrequency)
                decayRateMean, decayRateSe = task2_globalPeakFit(L, ncells, npart, samples, goodPeakAmplitudes, goodPeakTimes, oscFreqMean, oscFreqSe)

                # Store to upper level arrays #Only works if only one values (L_, npart_, ncells_) apart from samples is a non-scalar
                oscFreqMeanM = np.append(oscFreqMeanM, oscFreqMean)
                oscFreqSeM = np.append(oscFreqSeM, oscFreqSe)
                decayRateMeanM = np.append(decayRateMeanM, decayRateMean)
                decayRateSeM = np.append(decayRateSeM, decayRateSe)


                # Add to L level Array



                # end of ncells loop
            # end of npart loop
        
        # Get Angular Frequency values for 

        #end of L loop
    
    # Back to main function

    textChanging = "#Cells"
    textPlot = "Angular Freq"
    x = ncells_
    y = oscFreqMeanM
    stdErrs = oscFreqSeM
    plotHelper(textChanging,textPlot,x,y,stdErrs,L,ncells,npart,samples)

    textChanging = "#Cells"
    textPlot = "Damping Rate"
    x = ncells_
    y = decayRateMeanM
    stdErrs = decayRateSeM
    plotHelper(textChanging,textPlot,x,y,stdErrs,L,ncells,npart,samples)






