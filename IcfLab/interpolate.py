import numpy as np
import icf
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

# We will make a simple data set of discrete x points
xdata = np.array([1,2,3,4,5,6,7,8,9,10])

# And a y array which is just the square of these points
ydata = xdata**2


# Generate a function, f_square, which interpolates the data
# Here we use a linear interpolation, but take a look at the
# scipy documentation for other options which may give better
# results


f_square  = interpolate.interp1d(xdata, ydata, kind="linear")

# We now define a different set of points, on which we want to 
# know the value of the curve defined by xdata, ydata (in this
# case, we know its just a simple parabola!)
xdata_shifted = np.array([1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5])

# We use the interpolation function to find y points for each 
# x point in the array, these points should lie 'on', or close to, 
# the curve defined by xdata, ydata
ydata_shifted = f_square(xdata_shifted)


# We now just plot the original data and interpolated data
plt.plot(xdata, ydata, 'o', label="Data")
plt.plot(xdata_shifted, ydata_shifted, 'x', label="Interpolation")

plt.legend(bbox_to_anchor=(0.5, 1))

print(xdata_shifted)
print(ydata_shifted)

plt.show()



