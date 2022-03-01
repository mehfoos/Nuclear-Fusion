import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# read in data
intensity  = np.loadtxt('intensity.dat')	# 2D array of CCD counts data
wavelength = np.loadtxt('lambda.dat')		# 2D array of wavelength data
radius     = np.loadtxt('radius.dat')		# 1D array of major radii
angle      = np.loadtxt('angle.dat')		# 1D array of scattering angles

