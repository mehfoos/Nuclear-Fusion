'''
ICF Lab
'''

from matplotlib import colors
#from matplotlib.lines import _LineStyle
import numpy as np
import icf
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.interpolate import interp1d # SciPy 1D interpolant
from scipy.optimize import curve_fit

if __name__ == "__main__":
    global gaussianOrder

    # Gaussian generation
    def gaussian(x, *params):
        A = params[0]
        x0 = params[1]
        c = params[2]
        A2 = params[3]
        x02 = params[4]
        c2 = params[5]
        A3 = params[6]
        x03 = params[7]
        c3 = params[8]
        y0 = params[9]
        #return y0 + A*np.exp(-(x-x0)**2/(2*c*c))
        #eturn y0 + A*np.exp(-(x-x0)**2/(2*c*c)) + A2*np.exp(-(x-x02)**2/(2*c2*c2))
        return y0 + A*np.exp(-((x-x0)/(np.sqrt(2)*c))**gaussianOrder) + A2*np.exp(-((x-x02)/(np.sqrt(2)*c2))**gaussianOrder) + A3*np.exp(-((x-x03)/(np.sqrt(2)*c3))**gaussianOrder)

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
    print("comboVisibleEnergies = ", comboVisibleEnergies)

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
    #print(E_data)

    interpFilterAttenuation = interp1d(filtAttenEnergy, filtAttenCor, fill_value="extrapolate")
    T_E = interpFilterAttenuation(E_data)

    interpCrystalReflectivity = interp1d(filtCrystalReflEnergy, filtCrystalReflCor, fill_value="extrapolate")
    R_E = interpCrystalReflectivity(E_data)

    interpPhotocathodeSensitivity = interp1d(filtPhotocatSensEnergy, filtPhotocatSensCor, fill_value="extrapolate")
    C_E = interpPhotocathodeSensitivity(E_data)

    # I Corrected
    I_meas = ArK_CrosssectionData
    I_corr = I_meas/(T_E*R_E*C_E)

    # Plot Correction
    f5, ax5 = plt.subplots()
    ax5.plot(E_data, I_meas, 'b--', label=r"$I_{measured}$")
    ax5.plot([], [], color='red', label=r"$I_{corrected}$ right axis")
    ax5r = ax5.twinx()  # instantiate a second axes that shares the same x-axis
    ax5r.plot(E_data, I_corr, color='red', label=r"$I_{corrected}$")
    ax5.set_ylabel("I Corrected [relative grayscale intesnity]")
    ax5.set_xlabel("Energy [eV]")
    ax5.set_title("Corrected Spectrum")
    ax5.legend()

    ''' 
    Gaussian fit
    '''
    # Gaussian fit
    gaussianOrder = 2
    
    Ly_B = 3935.6
    He_B = 3683.7

    # Trim data
    #E_data_Trim = E_data[310:361] # full: 300:450, First 310:361; second 360:450
    #I_corr_Trim = I_corr[310:361]
    E_data_Trim = E_data[360:450] # full: 300:450, First 310:361; second 360:450
    I_corr_Trim = I_corr[360:450]
    #plt.plot(E_data)

    # Guess parameters: A[0], x0[1], c[2], A2[3], x02[4], c2[5], y0[6]
    #guess = [4620-670, He_B, 1, 2420-670, Ly_B, 2,  500] # Fitted both with 2 gaussians, but not well
    #guess = [2000, 3675, 1, 3000, He_B, 10, 3500, 3686, 5,  450] # fit for HeB alone
    guess = [1500-500, 3870, 13, 1600-500, 3925, 3, 2500-500, Ly_B, 10,  500] # fit for LyB and preceeding peak alone
    #[max(ydata)-min(ydata),(max(xdata)-min(xdata))/2,0.25,percentile(ydata,25)]
    print("Our initial guess is", guess)
    popt, pcov = curve_fit(gaussian, E_data_Trim, I_corr_Trim, p0=guess, maxfev = 10000)
    yfit = gaussian(E_data_Trim, *popt)
    guessFit = gaussian(E_data_Trim, *guess)
    print("Fit parameters : ", popt)
    gaussianCentres = [popt[1], popt[4], popt[7]]
    #plt.close("all")
    rsquared = icf.r_squared(I_corr_Trim, yfit)
    print("R^2 = ", rsquared)
    f6, ax6 = plt.subplots(2, sharex=True)
    #ax6[1].plot(E_data, I_corr, lw=1, label="Data", linestyle='dotted', color="black")	
    ax6[1].plot(E_data_Trim, I_corr_Trim, lw=1, label="Data", linestyle='dotted', color="black")	
    ax6[1].plot(E_data_Trim, yfit, lw=1, label=r"Fit; Order ="+ str(gaussianOrder) + " $R^{2}$=" + f"{rsquared:.2f}")
    ax6[1].vlines(x=gaussianCentres,ymin=0,ymax=4000,colors=['gray','pink','green'], label=r"Gaussian centres gry1, pnk2, grn3")
    #ax6[1].plot(E_data_Trim, guessFit, lw=1, label=r"Guess Fit")
    ax6[0].plot(E_data_Trim, I_corr_Trim-yfit, lw=1, label="Residual; Order ="+str(gaussianOrder)) 
    #ax6[1].vlines(x=comboEnergies,ymin=0,ymax=4000,colors='gray', label=r"Visible Resonance transition energies")
    ax6[1].legend()
    ax6[0].legend()				
    ax6[0].yaxis.tick_right()
    ax6[0].set_ylabel("Residual")
    ax6[1].set_ylabel("Grey Intensity")
    ax6[1].set_xlabel("Energy [eV]")
    ax6[0].set_title("Order: " + str(gaussianOrder))
    f6.suptitle("Gaussian fit")
    f6.subplots_adjust(hspace=0)
    
    # FWHM
    indexC = 2
    fwhm1 = 2*np.sqrt(2)*(np.log(2)**(1.0/gaussianOrder))*popt[indexC]
    indexC = 5
    fwhm2 = 2*np.sqrt(2)*(np.log(2)**(1.0/gaussianOrder))*popt[indexC]
    indexC = 8
    fwhm3 = 2*np.sqrt(2)*(np.log(2)**(1.0/gaussianOrder))*popt[indexC]
    fwhm3_error = np.sqrt(pcov[indexC][indexC])
    #fwhm_error = np.sqrt(pcov[indexC][indexC])
    print("fwhm1 is " + str(fwhm1))
    print("fwhm2 is " + str(fwhm2))
    print("fwhm3 is " + str(fwhm3))
    print("fwhm3_error is " + str(fwhm3_error))

    # Results FWHM (Storing previously calculated values)
    He_B_Fwhm = 24.247888831900895
    Ly_B_Fwhm = 25.53321741504301
    Ly_B_Fwhm_Error = 0.4057924412418137
    He_B_Fwhm_Error = 0.3175458198513755


    # Spatial Resolution Array
    spatial_Res_Array = np.array([400, 150, 150, 400]) # um
    spatial_Res = 150.0

    # R
    R_Fwhm_um = np.sqrt(100.0**2 + spatial_Res**2 + 1.0 + 60.0**2)
    print("R_Fwhm_um " + str(R_Fwhm_um))
    R_Fwhm_pixel = R_Fwhm_um/filmRes
    print("R_Fwhm_pixel " + str(R_Fwhm_pixel))
    Ly_B_Pixel = peakCentre[5]
    He_B_Pixel = peakCentre[3]
    #print(peakCentre[5])
    #print(polyEnergyFit(peakCentre[5]))
    #print(peakCentre[3])
    #print(polyEnergyFit(peakCentre[3]))
    Ly_B_R_Fwhm_ev = polyEnergyFit(Ly_B_Pixel+R_Fwhm_pixel/2)-polyEnergyFit(Ly_B_Pixel-R_Fwhm_pixel/2)
    He_B_R_Fwhm_ev = polyEnergyFit(He_B_Pixel+R_Fwhm_pixel/2)-polyEnergyFit(He_B_Pixel-R_Fwhm_pixel/2)
    print("Ly_B_R_Fwhm_ev " + str(Ly_B_R_Fwhm_ev))
    print("He_B_R_Fwhm_ev " + str(He_B_R_Fwhm_ev))

    # Sfwhm, errors
    Ly_B_S_Fwhm = np.sqrt(Ly_B_Fwhm**2 - Ly_B_R_Fwhm_ev**2)
    He_B_S_Fwhm = np.sqrt(He_B_Fwhm**2 - He_B_R_Fwhm_ev**2)

    print("Ly_B_S_Fwhm " + str(Ly_B_S_Fwhm))
    print("He_B_S_Fwhm " + str(He_B_S_Fwhm))

    # Power in spectral lines
    print(str(popt[6]))
    He_B_Amplitude = 3737.9798524770104
    Ly_B_Amplitude = 1760.7344279131992

    He_B_P = He_B_Amplitude * (He_B_S_Fwhm/2) / 0.3989
    Ly_B_P = Ly_B_Amplitude * (Ly_B_S_Fwhm/2) / 0.3989
    print("He_B_P: " + str(He_B_P))
    print("Ly_B_P: " + str(Ly_B_P))


    # Show all plots
    #plt.close("all")
    plt.show()