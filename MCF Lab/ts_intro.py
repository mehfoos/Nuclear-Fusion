import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian(x, *params):
    A = params[0]
    x0 = params[1]
    c = params[2]
    y0 = params[3]
    #return y0 + A*np.exp(-(x-x0)**2/(2*c*c))
    return y0 + A*np.exp(-((x-x0)/(np.sqrt(2)*c))**order)

def find_closest(data, v):
	""" Find the value in data closest to v """
	return (np.abs(data-v)).argmin()

# read in data
intensity  = np.loadtxt(r'MCF Lab/intensity.dat')	# 2D array of CCD counts data
wavelength = np.loadtxt(r'MCF Lab/lambda.dat')		# 2D array of wavelength data
radius     = np.loadtxt(r'MCF Lab/radius.dat')		# 1D array of major radii
angle      = np.loadtxt(r'MCF Lab/angle.dat')		# 1D array of scattering angles

# Plot intensity data
fig, ax = plt.subplots()
im = ax.imshow(intensity)
plt.colorbar(im,label="Intensity")
ax.set_title("intensity.dat: Imshow Plot")
ax.set_xlabel("X axis")
ax.set_ylabel("Y axis")

fig2, ax2 = plt.subplots()
cont = ax2.contour(intensity)
plt.colorbar(cont,label="Intensity")
ax2.set_title("intensity.dat: Contour Plot")
ax2.set_xlabel("X axis")
ax2.set_ylabel("Y axis")

# Plot wavelength data
fig3, ax3 = plt.subplots()
im3 = ax3.imshow(wavelength)
plt.colorbar(im3,label="Intensity/weightage")
ax3.set_title("lambda.dat: Imshow Plot")
ax3.set_xlabel("X axis")
ax3.set_ylabel("Y axis")

# Shape
size_wavelength = wavelength.shape
midPos = int(size_wavelength[1]/2)
print(size_wavelength)
print(midPos)

# Plot middle of wavelength data
fig4, ax4 = plt.subplots()
ax4.plot(wavelength[midPos,:], intensity[midPos,:], lw=0, marker='.', label=r"Middle Array index = "+ str(midPos))
ax4.set_title("Spectrum plot using actual wavelength from middle of array")
ax4.set_xlabel("Wavelength [assuming $10^{-10}$m]")
ax4.set_ylabel("Intensity [Photo Electrons]")
ax4.legend()

# Part 5
plt.close("all") # Ensure the plots up to now are all closed
global order
xdata = wavelength[midPos,:]
ydata = intensity[midPos,:]
xmin = find_closest(xdata, 6.0)
xmax = find_closest(xdata, 6.0)
xdata = xdata[xmin:xmax]
ydata = ydata[xmin:xmax]

# Guess and other config
guess = [max(ydata)-min(ydata),(max(xdata)-min(xdata))/2,0.25,percentile(ydata,25)]
print("Our initial guess is", guess)



plt.show()