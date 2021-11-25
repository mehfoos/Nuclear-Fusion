'''
ICF Lab
'''

import numpy as np
import icf
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.interpolate import interp1d # SciPy 1D interpolant

if __name__ == "__main__":

    # Resonance transition energies
    LyEnergies = np.array([3318.1, 3322.9, 3934.2, 3935.6, 4149.6, 4150.2])
    HeEnergies = np.array([3123.4, 3139.3, 3679.3, 3683.7, 3872.9, 3874.4]) 
    comboEnergies = np.concatenate([LyEnergies, HeEnergies])
    comboIncrementingEnergies = np.sort(comboEnergies)
    comboVisibleEnergies = np.array([3322.9, 3935.6, 3123.4, 3139.3, 3683.7, 3874.4]) 
    comboVisibleEnergies = np.sort(comboVisibleEnergies)
    # Working combos: 0.5/-1100,
    scaleEnergy = 0.39
    offsetEnergy = -1070
    comboScaledEnergies = (comboVisibleEnergies * scaleEnergy) + offsetEnergy
    print("comboIncrementingEnergies = ", comboIncrementingEnergies)

    # streak camera temporal resolution; all units in micrometer
    d = 0.25 * 1000 #um
    v = 1.0/(63.0 / 1000.0) #um/ps
    M = 1.24
    deltaS = 150 #um
    filmRes = 60 #um
    deltaTau = np.sqrt(((d*M)**2) + ((deltaS)**2))/v # ps/pixel
    print("deltaTau (streak camera temporal resolution): ", deltaTau)

    # Conversion to pixel width
    pixelWidth = M*deltaTau*v/filmRes
    print("pixelWidth: ", pixelWidth)

    # Import the cross-section data
    fileName = r"C:\Users\mly509\OneDrive - University of York\Documents\Fusion Lab Practicals\Week2\ArK Cross-section Values.csv"
    ArK_pixelData, ArK_CrosssectionData = icf.load_2col(fileName)
    ArK_CrosssectionData = np.flip(ArK_CrosssectionData)
    xoffset1 = 40.0
    xoffset2 = 579.0
    xmin = icf.find_closest(ArK_pixelData, xoffset1)
    xmax = icf.find_closest(ArK_pixelData, xoffset2)
    ArK_pixelData = ArK_pixelData[xmin:xmax]
    ArK_CrosssectionData = ArK_CrosssectionData[xmin:xmax]

    # Plot the known energy lines
    f1, ax1 = plt.subplots()
    ax1.vlines(x=comboScaledEnergies,ymin=0,ymax=100,colors='gray', label=r"resonance transition energies")

    # Plot the raw data
    ax1.plot(ArK_pixelData, ArK_CrosssectionData, label=r"SSC Data Raw")
    ax1.set_title("SSC Raw Data Plot")
    ax1.set_xlabel("Pixel")
    ax1.set_ylabel("Grey Intensity")

    # Find peaks
    ArK_Peaks, _ = find_peaks(ArK_CrosssectionData,prominence=25)
    # Manual tweak
    ArK_Peaks = np.delete(ArK_Peaks,2)

    # Plot peaks
    ArK_Peaks_Offset = ArK_Peaks+xoffset1
    ax1.plot(ArK_Peaks_Offset,ArK_CrosssectionData[ArK_Peaks],marker='o', linestyle='none', label=r"Peaks")
    print("Peaks x coordinate in pixels: ", ArK_Peaks)
    ax1.legend()

    # Polynomial fit and plot
    energyFitCoeffs = np.polyfit(ArK_Peaks_Offset, comboVisibleEnergies, 2)
    polyEnergyFit = np.poly1d(energyFitCoeffs)
    f2, ax2 = plt.subplots()
    # print(polyEnergyFit(ArK_Peaks_Offset))
    rsquared = icf.r_squared(comboVisibleEnergies, polyEnergyFit(ArK_Peaks_Offset))
    print("R^2 = ", rsquared)
    ax2.plot(ArK_Peaks_Offset, comboVisibleEnergies, 'x', label=r"Visible Resonance transition energies")
    ax2.plot(ArK_pixelData, polyEnergyFit(ArK_pixelData), linestyle='dashed', label=r"Peaks")
    ax2.set_title(r"Polyfit: " + '{}'.format(polyEnergyFit) + "\n $R^{2}$=" + f"{rsquared:.4f}")
    
    # Apply Pixel to Energy transformation
    f3, ax3 = plt.subplots()
    ax3.vlines(x=comboVisibleEnergies,ymin=0,ymax=250,colors='gray', label=r"Visible Resonance transition energies")
    ax3.plot(polyEnergyFit(ArK_pixelData), ArK_CrosssectionData, label=r"SSC Data")
    ax3.plot(polyEnergyFit(ArK_Peaks_Offset),ArK_CrosssectionData[ArK_Peaks],marker='o', linestyle='none', label=r"Peaks")
    ax3.set_title("SSC Raw Data Plotted against fitted Energy Axis")
    ax3.set_xlabel("Energy [eV]")
    ax3.set_ylabel("Grey Intensity")
    ax3.legend()

    # New Arrays
    peakCentre = ArK_Peaks_Offset
    peakEnergy = comboVisibleEnergies
    peakEnergyFit = polyEnergyFit(peakCentre)

    # Step 2.e
    f4, ax4 = plt.subplots()
    ax4.plot(peakCentre, peakEnergy, "x", label=r"Peaks on pixels vs Energy")
    ax4.plot(peakCentre,polyEnergyFit(peakCentre), linestyle='dashed', label=r"Second Order Polynomial Fit")
    ax4.set_xlabel("Peaks in Pixel Location")
    ax4.set_ylabel("Peaks in eV")
    ax4.set_title("Pixels to Energy Transformaiton and fit")
    ax4.legend()

    '''
    Step 3; interpolation
    '''
    # Import all the .dat files for instrument corrections
    fileNameFilterAttenuation = r"C:\Users\mly509\OneDrive - University of York\Documents\Fusion Lab Practicals\Week2\Filter_xray6772.dat"
    filtAttenEnergy, filtAttenCor = icf.load_2col(fileNameFilterAttenuation)
    #print(filtAttenCor)

    fileNameCrystalReflectivity = r"C:\Users\mly509\OneDrive - University of York\Documents\Fusion Lab Practicals\Week2\reflectivity.dat"
    filtCrystalReflEnergy, filtCrystalReflCor = icf.load_2col(fileNameCrystalReflectivity)
    #print(filtCrystalReflCor)

    fileNamePhotocathodeSensitivity = r"C:\Users\mly509\OneDrive - University of York\Documents\Fusion Lab Practicals\Week2\photocathode.dat"
    filtPhotocatSensEnergy, filtPhotocatSensCor = icf.load_2col(fileNamePhotocathodeSensitivity)
    #print(filtPhotocatSensCor)

    # Interpolation
    E_data = polyEnergyFit(ArK_pixelData)

    interpFilterAttenuation = interp1d(filtAttenEnergy, filtAttenCor)
    T_E = interpFilterAttenuation(E_data)

    interpCrystalReflectivity = interp1d(filtCrystalReflEnergy, filtCrystalReflCor)
    T_E = interpCrystalReflectivity(E_data)

    interpPhotocathodeSensitivity = interp1d(filtPhotocatSensEnergy, filtPhotocatSensCor)
    T_E = interpPhotocathodeSensitivity(E_data)

    # 


    # Show all plots
    #plt.show()
