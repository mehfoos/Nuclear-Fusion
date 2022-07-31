import wave
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# This function will calculate an R^2 values for your fit	
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

def gaussian(x, *params):
    A = params[0]
    x0 = params[1]
    c = params[2]
    # y0 = params[3]
    return A*np.exp(-(x-x0)**2/(2*c*c))

def gaussian3(x, *params):
    A1 = params[0]
    x01 = params[1]
    c1 = params[2]
    A2 = params[3]
    x02 = params[4]
    c2 = params[5]
    A3 = params[6]
    x03 = params[7]
    c3 = params[8]
    sqrt2 = np.sqrt(2)
    # y0 = params[9]
    return +A1*np.exp(-(x-x01)**2/(2*c1*c1)) - A2*np.exp(-(x-x02)**2/(2*c2*c2)) - A3*np.exp(-(x-x03)**2/(2*c3*c3))
    # return A1*np.exp(-(x-x01)**2/((sqrt2*c1)**2)) + A2*np.exp(-(x-x02)**2/((sqrt2*c2)**2)) + A3*np.exp(-(x-x03)**2/((sqrt2*c3)**2))
    # return y0 + A*np.exp(-(x-x0)**2/(2*c*c))


def find_closest(data, v):
	""" Find the value in data closest to v """
	return (np.abs(data-v)).argmin()

# read in data
intensity  = np.loadtxt(r'MCF Lab/intensity.dat')	# 2D array of CCD counts data
wavelength = np.loadtxt(r'MCF Lab/lambda.dat')		# 2D array of wavelength data
radius     = np.loadtxt(r'MCF Lab/radius.dat')		# 1D array of major radii
angle      = np.loadtxt(r'MCF Lab/angle.dat')		# 1D array of scattering angles

### PART 3: Extracting and Plotting Data

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

# Plot using radius.dat to correct
fig5, ax5 = plt.subplots()
im5 = ax5.pcolor(wavelength, radius, intensity, shading='auto')
plt.colorbar(im5, label="Intensity in CCD")
ax5.set_title("CCD Intensity for Radius vs Wavelength")
ax5.set_xlabel("Wavelength [$10^{-10}$m]")
ax5.set_ylabel("Radius (from Tokomak centre) [m]")

### END OF PART 3 (Extracting and Plotting Data)



### START OF PART 4 (Data analysis: measuring T_e)

plt.close("all") # Ensure the plots up to now are all closed

# Get attributes that can be handy for guess
rangeWavelength = np.max(wavelength[midPos,:]) - np.min(wavelength[midPos,:])
maxI = intensity[midPos,:].argmax()
peakWavelength = wavelength[midPos,maxI]
amplitudeMax = max(intensity[midPos,:])
print("peak Intensity Wavelength: ", peakWavelength)

# Gaussian Fit, attempt1; 1 gaussian
guess1 = [amplitudeMax, peakWavelength, rangeWavelength*0.1]
guessFit1 = gaussian(wavelength[midPos,:], *guess1)
popt1, pcov1 = curve_fit(gaussian, wavelength[midPos,:], intensity[midPos,:], guess1) #, maxfev = 10000)
print("popt1", popt1)
yfit1 = gaussian(wavelength[midPos,:], *popt1)
figP4_1, axP4_1 = plt.subplots()
axP4_1.plot(wavelength[midPos,:], intensity[midPos,:], lw=0, marker='.', label=r"Data; index = "+ str(midPos))
axP4_1.plot(wavelength[midPos,:], yfit1, label=r"Gaussian Fit")
axP4_1.plot(wavelength[midPos,:], guessFit1, lw=1, linestyle = "dashed", label=r"Guess Fit")
axP4_1.set_title("Spectrum Intensity: Gaussian fit attempt1")
axP4_1.set_xlabel("Wavelength [$10^{-10}$m]")
axP4_1.set_ylabel("CCD Intensity")
axP4_1.legend()


# Gaussian Fit, 3 gaussians; attempt2
guess2 = [amplitudeMax*1.2, peakWavelength, rangeWavelength*0.1, \
    amplitudeMax, peakWavelength-150, rangeWavelength*0.008, \
        amplitudeMax, peakWavelength+150, rangeWavelength*0.03]
guessFit2 = gaussian3(wavelength[midPos,:], *guess2)
popt2, pcov2 = curve_fit(gaussian3, wavelength[midPos,:], intensity[midPos,:], guess2, maxfev = 20000)
print("guess2", guess2)
print("popt2", popt2)
yfit2 = gaussian3(wavelength[midPos,:], *popt2)
figP4_2, axP4_2 = plt.subplots()
axP4_2.plot(wavelength[midPos,:], intensity[midPos,:], lw=0, marker='.', label=r"Data; index = "+ str(midPos))
axP4_2.plot(wavelength[midPos,:], yfit2, label=r"Gaussian Fit")
axP4_2.plot(wavelength[midPos,:], guessFit2, lw=1, linestyle = "dashed", label=r"Guess Fit")
axP4_2.set_title("Spectrum Intensity: Gaussian fit attempt2")
axP4_2.set_xlabel("Wavelength [$10^{-10}$m]")
axP4_2.set_ylabel("CCD Intensity")
axP4_2.legend()

# Gaussian Fit, using 1 gaussian but optimal paramaters for 3 gaussians; attempt 3
# plt.close("all")
guess3 = guess2
popt3, pcov3 = curve_fit(gaussian3, wavelength[midPos,:], intensity[midPos,:], guess3)
yfit3_1 = gaussian(wavelength[midPos,:], *popt3)
figP4_3, axP4_3 = plt.subplots()
axP4_3.plot(wavelength[midPos,:], intensity[midPos,:], lw=0, marker='.', label=r"Data; index = "+ str(midPos))
axP4_3.plot(wavelength[midPos,:], yfit3_1, label=r"Gaussian Fit3")
axP4_3.set_title("Spectrum Intensity: Gaussian fit attempt3")
axP4_3.set_xlabel("Wavelength [$10^{-10}$m]")
axP4_3.set_ylabel("CCD Intensity")
axP4_3.legend()

# plt.show()

# Known values
λi = 694.3e-9 # m; 694.3nm
θ = angle[midPos] # Slice angle; assuming radian
kb = 1.381e-23 # J/K
me = 9.109e-31 # kg
c = 3e8 # m/s
ev_to_K = 11606 # K/eV


# Curve-fitted σλ
σλ = popt3[2]*1e-10 # Centre spectrum's "c", for gaussian
σλ_err = 1e-10*np.sqrt(np.diag(pcov3))[2] # / np.sqrt(np.size(wavelength[midPos,:]))

# Calculated βth
βth = σλ / (λi * np.sqrt(2) * np.sin(θ/2))
βth_err = σλ_err / (λi * np.sqrt(2) * np.sin(θ/2))

# Calculated Te
Te_K = ((βth**2)*me*(c**2))/(2*kb)
Te_eV = Te_K/ev_to_K 
Te_eV_err = 2*βth_err*βth*me*(c**2)/(2*kb)/ev_to_K 

# Print everything
print("λi: ", λi)
print("θ: ", θ)
print("kb: ", kb)
print("me: ", me)
print("c: ", c)
print("ev_to_K: ", ev_to_K)
print("σλ: ", σλ)
print("σλ_err: ", σλ_err)
print("βth: ", βth)
print("βth_err: ", βth_err)
print("Te_K: ", Te_K)
print("Te_eV: ", Te_eV)
print("Te_eV_err: ", Te_eV_err)

# Build Temperature profile for full data
plt.close("all")

# Initialise array for results
Te_eV_arr = np.zeros_like(angle)
Te_eV_sdErr_arr = np.zeros_like(angle)
sdErr_arr = np.zeros_like(angle)
rSq_arr = np.zeros_like(angle)

for iPos in range(len(angle)):

    # Get attributes that can be handy for guess
    rangeWavelength = np.percentile(wavelength[iPos,:],99) - np.percentile(wavelength[iPos,:],10)
    maxI = intensity[iPos,:].argmax()
    peakWavelength = wavelength[iPos,maxI]
    amplitudeMax = max(intensity[iPos,:])

    # Gaussian fit
    guess = [amplitudeMax*1.2, peakWavelength, rangeWavelength*0.1, \
        amplitudeMax, peakWavelength-150, rangeWavelength*0.008, \
            amplitudeMax, peakWavelength+150, rangeWavelength*0.03]
    

    try:
        popt, pcov = curve_fit(gaussian3, wavelength[iPos,:], intensity[iPos,:], guess, \
            bounds=([amplitudeMax, np.NINF, np.NINF, np.NINF, np.NINF, np.NINF, np.NINF, np.NINF, np.NINF],\
                [amplitudeMax*1.5, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]), \
                     maxfev = 5000)
    except:
        Te_eV_arr[iPos] = np.nan
        sdErr_arr[iPos] = np.nan
        rSq_arr[iPos] = np.nan
        continue
    sdErr_pos = np.sqrt(np.diag(pcov))[2] # / np.sqrt(np.size(wavelength[iPos,:]))
    rSq_arr[iPos] = r_squared(intensity[iPos,:], gaussian3(wavelength[iPos,:], *popt))

    # if R^2 is too low, then skip
    if rSq_arr[iPos] < 0.9:
        Te_eV_arr[iPos] = np.nan
        sdErr_arr[iPos] = np.nan
        continue


    # Position specific values for calculation
    θ_pos = angle[iPos] # Slice angle; assuming radian
    σλ_pos = popt[2]*1e-10 # Centre spectrum's "c", for gaussian

    # Calculation
    βth_pos = σλ_pos / (λi * np.sqrt(2) * np.sin(θ_pos/2))
    βth_pos_err = sdErr_pos*1e-10 / (λi * np.sqrt(2) * np.sin(θ_pos/2))
    Te_K_pos = ((βth_pos**2)*me*(c**2))/(2*kb)
    Te_eV_pos = Te_K_pos/ev_to_K

    # Add to array
    Te_eV_arr[iPos] = Te_eV_pos
    sdErr_arr[iPos] = sdErr_pos
    Te_eV_sdErr_arr[iPos] = 2*βth_pos_err*βth_pos*me*(c**2)/(2*kb)/ev_to_K 


    # # Debug
    # print("iPos: ", iPos)
    # print("Te_eV_arr[iPos]: ", Te_eV_arr[iPos])
    # print("sdErr_arr[iPos]: ", sdErr_arr[iPos])
    # yfit3_1 = gaussian(wavelength[iPos,:], *popt)
    # guessFit = gaussian3(wavelength[midPos,:], *guess)
    # figP4_3, axP4_3 = plt.subplots()
    # axP4_3.plot(wavelength[iPos,:], intensity[iPos,:], lw=0, marker='.', label=r"Data; index = "+ str(iPos))
    # axP4_3.plot(wavelength[iPos,:], yfit3_1, label=r"Gaussian Fit3" + " $R^{2}$=" + f"{rSq_arr[iPos]:.2f}")
    # axP4_3.plot(wavelength[midPos,:], guessFit, lw=1, linestyle = "dashed", label=r"Guess Fit")
    # axP4_3.set_title("Spectrum Intensity: Gaussian fit attempt3")
    # axP4_3.set_xlabel("Wavelength [$10^{-10}$m]")
    # axP4_3.set_ylabel("CCD Intensity")
    # axP4_3.legend()
    # plt.show()

figP4_4, axP4_4 = plt.subplots()
axP4_4.errorbar(radius,Te_eV_arr,yerr=Te_eV_sdErr_arr,label="Shot 17447", ecolor="red")# ,uplims=True, lolims=True)#,linestyle='none')
axP4_4.set_title("Electron temperature profile for shot 17447")
axP4_4.set_xlabel("Radius [m]")
axP4_4.set_ylabel("$T_{e}$ [eV]")
axP4_4.legend()
axP4_4.set_ylim(0,600)
plt.show()



### END OF PART 4 (Data analysis: measuring T_e)