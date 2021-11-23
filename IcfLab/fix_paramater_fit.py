import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import icf

def gaussian(x, *params):

	#Fix parameters here	
		
	A = params[0]
	x0 = params[1]
	c = params[2]
	
	A2 = params[3]
	x02 = params[4]
	c2 = params[5]
		
	return A*np.exp(-(x-x0)**2/(2*c*c)) + A2*np.exp(-(x-x02)**2/(2*c2*c2))
	

# First we will make a Gaussian
xdata = np.linspace(0,4,100)
ydata = gaussian(xdata, 1, 1.5, 0.4, 3, 2, 0.2 )

# and add some noise
for i in range(len(ydata)):
	ydata[i] +=  0.4*(np.random.random_sample()-0.5)


# This array holds our inital guesses for all the parameters
guesses = np.array([1, 1, 1, 1, 2, 1])

# This 'tuple' defines the minimum and maxmimum allowed values for each parameter
# The first list gives the lowest allowed value for each of the 6 fitting parameters, and the 
# second list lists the maximum allowed values.  Note that the function needs a finite (non-zero)
# range.
# Here we have fixed the fifth parameter (the x0 for the second Gaussian), is fixed to be very close to
# 2, whereas all other parameters are free to vary from 0-10.
# Note that the guesses must lie within the ranges defined here

ranges = [[0,0,0,0,1.9999,0], [10,10,10,10,2.0001,10]]

# This is our fitting command with bounds on the parameters
#popt, pcov = curve_fit(gaussian, xdata, ydata, guesses, bounds=ranges)
popt, pcov = curve_fit(gaussian, xdata, ydata, guesses)

print(popt)

# Grab a copy of the best fit line
yfit = gaussian(xdata, *popt)


# This just plots the result

icf.fit_plot(xdata, ydata, yfit)

