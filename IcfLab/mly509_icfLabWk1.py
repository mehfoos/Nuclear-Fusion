'''
ICF Lab
'''

import numpy as np
from numpy import core
from numpy.lib.function_base import percentile
from scipy.optimize import curve_fit
import icf
import matplotlib.pyplot as plt
import statistics as s

global order
# This function will generate a perfect Gaussian 
# You can replace this function with any other you want to fit with...
def gaussian(x, *params):
    A = params[0]
    x0 = params[1]
    c = params[2]
    y0 = params[3]
    #return y0 + A*np.exp(-(x-x0)**2/(2*c*c))
    return y0 + A*np.exp(-((x-x0)/(np.sqrt(2)*c))**order)

if __name__ == "__main__":
    # Constants, configs
    nFiles = 4
    PIXEL_TO_UM = 60.0
    PIXEL_TO_MM = PIXEL_TO_UM/1000.0
    UM_TO_MM = 1.0/1000.0
    nOrders = 4
    rsquared_best = np.zeros(nFiles)
    fwhm = np.zeros(nFiles)
    fwhm_error = np.zeros(nFiles)
    indexC = 2
    order_best = np.zeros(nFiles)
    pinhole_size = 15*UM_TO_MM
    microchannel_size = 27*UM_TO_MM
    microchannel_error = 6*UM_TO_MM
    digitisation_size = 60*UM_TO_MM
    shotSequenceTime = 0.2 #ns
    


    # Import the csv data
    fileNameCommon = r"C:\Users\mly509\OneDrive - University of York\Documents\Fusion Lab Practicals\Week1\CrossSection2Line"
    
    for iFile in range(nFiles):
        fileNo = iFile+1
        fileName = fileNameCommon + str(fileNo) + "Data.csv"
        xdata, ydata = icf.load_2col(fileName)

        # Convert from pixels to mm
        xdata *= PIXEL_TO_MM

        # Trim
        xmin = icf.find_closest(xdata, 0.0)
        xmax = icf.find_closest(xdata, 6.0)
        xdata = xdata[xmin:xmax]
        ydata = ydata[xmin:xmax]

        '''
        Fit section
        '''

        # Guess and other config
        guess = [max(ydata)-min(ydata),(max(xdata)-min(xdata))/2,0.25,percentile(ydata,25)]
        print("Our initial guess is", guess)

        # Plotting based on icf.fit_plot
        f, axarr = plt.subplots(2, sharex=True)
        axarr[1].plot(xdata, ydata, lw=1, label="Data", linestyle='dotted', color="black")	

        # Loop through different orders
        for iOrder in range(nOrders):
            # Gaussian order
            order = 2*(iOrder+1)

            # Try fit
            popt, pcov = curve_fit(gaussian, xdata, ydata, p0=guess)

            # Get fit
            yfit = gaussian(xdata, *popt)

            #Print
            for i in range(len(popt)):
                print ("Order: ",order,"; Parameter",i,":",popt[i],"+/-",np.sqrt(pcov[i][i]))
            print("Fit parameters : ", popt)
            print("Fit standard deviations : ", np.sqrt(np.diag(pcov)))
            rsquared = icf.r_squared(ydata, yfit)
            print("R^2 = ", rsquared)

            # Plot
            axarr[1].plot(xdata, yfit, lw=1, label=r"Fit; Order ="+ str(order) + " $R^{2}$=" + f"{rsquared:.2f}")
            axarr[0].plot(xdata, ydata-yfit, lw=1, label="Residual; Order ="+str(order)) 
            
            # Save the best R^2 FWCM
            if (rsquared > rsquared_best[iFile]):
                rsquared_best[iFile] = rsquared
                fwhm[iFile] = 2*np.sqrt(2)*(np.log(2)**(1.0/order))*popt[indexC]
                fwhm_error[iFile] = np.sqrt(pcov[indexC][indexC])
                order_best[iFile] = order

        # Plotting based on icf.fit_plot
        #f, axarr = plt.subplots(2, sharex=True)
        #axarr[1].plot(xdata, ydata, lw=2, label="Data")
        #axarr[1].plot(xdata, yfit, lw=2, label="Fit ")
        axarr[1].legend()
        axarr[0].legend()
        #axarr[0].plot(xdata, ydata-yfit, lw=2, label="Residual") 				
        axarr[0].yaxis.tick_right()
        axarr[0].set_ylabel("Residual")
        axarr[1].set_ylabel("Grey Intensity in Image line section")
        axarr[1].set_xlabel("Line Distance in mm")
        axarr[0].set_title("Best estimate order: " + str(order_best[iFile]) + "; FWHM = " + f"{fwhm[iFile]:.2f}" + "+- " + f"{fwhm_error[iFile]:.2f}")
        f.suptitle("Gaussian fit for line" + str(fileNo))
        f.subplots_adjust(hspace=0)
    
    # Calc results

    psfFwhm = np.sqrt(pinhole_size**2 + microchannel_size**2 + digitisation_size**2)
    print("psfFwhm: ", psfFwhm)
    coreDiameterFwhm = np.sqrt(s.mean(fwhm)**2 + psfFwhm**2)
    print("coreDiameterFwhm: ", coreDiameterFwhm)
    coreDiameterFwhmMagnificationAdjusted = coreDiameterFwhm/3.5
    print("coreDiameterFwhmMagnificationAdjusted: ", coreDiameterFwhmMagnificationAdjusted)
    coreDiameterError = np.sqrt(microchannel_error**2 + s.stdev(fwhm))/3.5
    print("coreDiameterError: ", coreDiameterError)

    # Print results
    print("fwhm values:", fwhm )
    print("fwhm errors:", fwhm_error )
    print("R^2 Values", rsquared_best)
    print("Orders best", order_best)
    print("Mean:", s.mean(fwhm))
    print("Variance: ", s.variance(fwhm))
    print("Population stdev: ", s.pstdev(fwhm))
    print("Sample stdev: ", s.stdev(fwhm))
    
    plt.show()