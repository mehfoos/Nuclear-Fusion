import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# read in data
intensity  = np.loadtxt('MCF Lab\intensity.dat')	# 2D array of CCD counts data
wavelength = np.loadtxt('MCF Lab\lambda.dat')		# 2D array of wavelength data
radius     = np.loadtxt('MCF Lab\radius.dat')		# 1D array of major radii
angle      = np.loadtxt('MCF Lab\angle.dat')		# 1D array of scattering angles

