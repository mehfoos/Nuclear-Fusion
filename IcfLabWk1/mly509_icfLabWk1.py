'''
ICF Lab
'''

import numpy as np
from scipy.optimize import curve_fit
import icf
import matplotlib.pyplot as plt

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

# Constants
PIXEL_TO_UM = 60.0
PIXEL_TO_MM = PIXEL_TO_UM/1000.0

# Import the csv data
xdata, ydata = icf.load_2col(r"C:\Users\mly509\OneDrive - University of York\Documents\Fusion Lab Practicals\Week1\CrossSectionLine1Data.csv")

# Plot the data to check
#plt.plot(xdata, ydata)
#plt.show()

# Convert from pixels to mm
xdata *= PIXEL_TO_MM

# Fit section

# Gaussian order
order = 4

# This does the fit, and returns the fit parameters and the covariances
guess = [1,1,1,25000]
print("Our initial guess is", guess)
popt, pcov = curve_fit(gaussian, xdata, ydata, p0=guess)

# Print fit parameters
for i in range(len(popt)):
	print ("Parameter",i,":",popt[i],"+/-",np.sqrt(pcov[i][i]))
	
print("Fit parameters : ", popt)
print("Fit standard deviations : ", np.sqrt(np.diag(pcov)))

# This data is therefore the best fit curve 
yfit = gaussian(xdata, *popt)

print("R^2 = ", icf.r_squared(ydata, yfit))

# This will plot the output, both the original data and the best fit, as well as a residual
# Note this is a special plotting routine written for the icf labs, hence the 'icf' prefix
# The source code can be found in icf.py if you want to copy/alter it

icf.fit_plot(xdata, ydata, yfit)
