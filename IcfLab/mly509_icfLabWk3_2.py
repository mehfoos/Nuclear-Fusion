import numpy as np
import icf
import re
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def sigmoid(x, *params):
    A = params[0]
    x0 = params[1]
    k = params[2]
    y0 = params[3]
    #y = (A / (1 + np.exp(-k*(x-x0)))) + y0
    y = (A / (1 + np.exp(x-x0))) + y0
    return (y)

if __name__ == "__main__":

    ### FILE IMPORT
    # Import the csv data
    fileName = r"C:\Users\mly509\OneDrive - University of York\Documents\Fusion Lab Practicals\Week3\capsuledata_Corrected.csv"
    time_delay = []
    diameter_um = []
    err_diameter_um = []
    time_delay_raw = []
    diameter_um_raw = []
    err_diameter_um_raw = []
    
    with open(fileName, 'r') as inputfile:
        lines = inputfile.readlines()
        
        for line in lines:
            data = [_f for _f in re.split('[,]', line) if _f]
            if icf.is_number(data[0]) and icf.is_number(data[3]) and icf.is_number(data[4]):
                time_delay.append(float(data[0]))
                diameter_um.append(float(data[3]))
                err_diameter_um.append(float(data[4]))	

            if icf.is_number(data[0]) and icf.is_number(data[1]) and icf.is_number(data[2]):
                time_delay_raw.append(float(data[0]))
                diameter_um_raw.append(float(data[1]))
                err_diameter_um_raw.append(float(data[2]))	

    # Convert to np.arrays
    diameter_um_raw = np.array(diameter_um_raw)
    err_diameter_um_raw = np.array(err_diameter_um_raw)
    time_delay_raw = np.array(time_delay_raw)
    diameter_um = np.array(diameter_um)
    err_diameter_um = np.array(err_diameter_um)
    time_delay = np.array(time_delay)

    ### Calculations

    # Unit converstion
    u_to_kg = 1.66054e-27
    um_to_m = 1.0/1000000.0

    # Known values
    P_d2 = 50.0 * 101325 # Pa 
    diameter_initial = 440.0 / 1000000.0
    m_d2_u = 2.014 # u
    m_d2 = m_d2_u*u_to_kg # kg
    N_a = 6.02214086e23 # avogradro constant
    R_gasConst = 8.31446261815324 # SI
    Temp_initial = 293 # K

    # Initial Calcs
    V_capsule_initial = (4/3)*np.pi*(diameter_initial/2)**3 # m^3
    n_moles = P_d2*V_capsule_initial/R_gasConst/Temp_initial # mole
    density_initial = m_d2*n_moles*N_a/V_capsule_initial
    number_density_initial = n_moles*N_a/V_capsule_initial
    mass_initial = m_d2*n_moles*N_a # kg

    ## Array calcs and plots
    # Figure setup
    f1, ax1 = plt.subplots(2, sharex=True)

    # Diameter
    #ax1[0].plot(time_delay, diameter_um, marker="x", label=r"Diameter [um]", color="blue", linestyle='none')
    ax1[0].errorbar(time_delay, diameter_um, yerr=err_diameter_um,marker=".", label='Diameter errors [+-um]',linestyle='none', color="red")
    ax1[0].plot(time_delay, diameter_um, marker="x", label='Diameter [um]',linestyle='none', color="blue")
    ax1[0].set_xlabel("Time delay [ns]")
    ax1[0].set_ylabel("Diameters (corrected) [um]")
    ax1[0].set_title("Diameters, Errors")
    ax1[0].legend()

    # Volumes
    volumes_um3 = (4/3)*np.pi*(diameter_um/2)**3 # um^3
    volumes_err_um3 = 3*(err_diameter_um/diameter_um)*volumes_um3
    #volumes_err_um3 = (4/3)*np.pi*(((diameter_um+err_diameter_um)/2)**3 - ((diameter_um)/2)**3) # um^3
    volumes_m3 = (4/3)*np.pi*(diameter_um*um_to_m/2)**3 # um^3
    volumes_err_m3 = 3*(err_diameter_um/diameter_um)*volumes_m3
    #print(str(volumes_um3))
    #print(str(volumes_err_um3))
    ax1[1].errorbar(time_delay, volumes_um3, yerr=volumes_err_um3,marker=".", label='Volume errors [+-$um^{3}$]',linestyle='none', color="red")
    ax1[1].plot(time_delay, volumes_um3, marker="x", label='Volumes [$um^{3}$]',linestyle='none', color="blue")
    ax1[1].set_xlabel("Time delay [ns]")
    ax1[1].set_ylabel("Volumes (corrected) [$um^{3}$]")
    ax1[1].set_title("corrected Volumes")
    ax1[1].legend()


    # Densities
    f2, ax2 = plt.subplots(2, sharex=True)
    densities_kgm3 = mass_initial/volumes_m3
    densities_err_kgm3 = 3*(err_diameter_um/diameter_um)*densities_kgm3
    ax2[0].errorbar(time_delay, densities_kgm3, yerr=densities_err_kgm3,marker=".", label='Densities errors [$kg.m^{-3}$]',linestyle='none', color="red")
    ax2[0].plot(time_delay, densities_kgm3, marker="x", label='Densities [$kg.m^{-3}$]',linestyle='none', color="blue")
    ax2[0].set_xlabel("Time delay [ns]")
    ax2[0].set_ylabel("Densities (corrected) [$kg.m^{-3}$]")
    ax2[0].set_title("corrected Densities, Errors")
    ax2[0].legend()

    # Pressures
    # Assuming Adiabatic
    gamma = 5.0/3.0
    adiabatic_PV_const = P_d2*(V_capsule_initial**gamma)
    pressures_Pa = adiabatic_PV_const/(volumes_m3**gamma)
    pressures_Pa_err = gamma*3*(err_diameter_um/diameter_um)*pressures_Pa
    ax2[1].errorbar(time_delay, pressures_Pa, yerr=pressures_Pa_err, marker=".", label='Pressure errors [$Pa}$]',linestyle='none', color="red")
    ax2[1].plot(time_delay, pressures_Pa, marker="x", label='Pressures [$Pa$]',linestyle='none', color="blue")
    ax2[1].set_xlabel("Time delay [ns]")
    ax2[1].set_ylabel("Pressure (corrected) [$Pa$]")
    ax2[1].set_title("corrected Pressures, Errors")
    ax2[1].legend()

    




    # Simple Results print
    print("V_capsule_initial: " + str(V_capsule_initial))
    print("n_moles: " + str(n_moles))
    print("density_initial: " + str(density_initial))
    print("number_density_initial: " + str(number_density_initial)),
    print("mass intial: ", str(m_d2*n_moles*N_a))


    ## Sigmoid fit: Distance
    #plt.close("all")
    f3, ax3 = plt.subplots(2, sharex=True)
    guess = [max(diameter_um), 1.6, 25, min(diameter_um)]
    popt, pcov = curve_fit(sigmoid, time_delay, diameter_um,guess, maxfev = 10000)
    print(str(popt))
    guessFit = sigmoid(time_delay, *guess)
    yfit = sigmoid(time_delay, *popt)
    ax3[0].plot(time_delay, diameter_um-yfit, lw=1, label="Residual") 
    ax3[1].plot(time_delay, diameter_um, marker="x", label='Diameter [um]',linestyle='none', color="blue")
    #ax3[1].plot(time_delay, guessFit, lw=1, label=r"Guess Fit")
    ax3[1].plot(time_delay, yfit, lw=1, label=r"SciPy Fit",color="orange")
    ax3[0].set_xlabel("Time delay [ns]")
    ax3[1].set_ylabel("Diameters (corrected) [um]")
    ax3[1].set_title("Diameters, Fit")
    ax3[0].set_title("Error in Fit")
    f3.suptitle("Sigmoid fit Diameter")
    ax3[1].legend()
    print("Sigmoid popt Distance: " + str(popt))

    ## Sigmoid fit: Densities
    #plt.close("all")
    f4, ax4 = plt.subplots(2, sharex=True)
    guess = [max(densities_kgm3), 1.6, 25, min(densities_kgm3)]
    popt, pcov = curve_fit(sigmoid, time_delay, densities_kgm3, guess, maxfev = 10000)
    print(str(popt))
    guessFit = sigmoid(time_delay, *guess)
    yfit = sigmoid(time_delay, *popt)
    ax4[0].plot(time_delay, densities_kgm3-yfit, lw=1, label="Residual") 
    ax4[1].plot(time_delay, densities_kgm3, marker="x", label='Density',linestyle='none', color="blue")
    #ax4[1].plot(time_delay, guessFit, lw=1, label=r"Guess Fit")
    ax4[1].plot(time_delay, yfit, lw=1, label=r"SciPy Fit", color="orange")
    ax4[0].set_xlabel("Time delay [ns]")
    ax4[1].set_ylabel("Density [$kg.m^{-3}$]")
    ax4[1].set_title("Densities, Fit")
    ax4[0].set_title("Error in Fit")
    f4.suptitle("Sigmoid fit Densities")
    ax4[1].legend()
    print("Sigmoid popt Densities: " + str(popt))

    ## Sigmoid fit: Pressures
    #plt.close("all")
    f5, ax5 = plt.subplots(2, sharex=True)
    guess = [max(pressures_Pa), 1.6, 25, min(pressures_Pa)]
    popt, pcov = curve_fit(sigmoid, time_delay, pressures_Pa, guess, maxfev = 10000)
    print(str(popt))
    guessFit = sigmoid(time_delay, *guess)
    yfit = sigmoid(time_delay, *popt)
    ax5[0].plot(time_delay, pressures_Pa-yfit, lw=1, label="Residual") 
    ax5[1].plot(time_delay, pressures_Pa, marker="x", label='Pressure',linestyle='none', color="blue")
    #ax5[1].plot(time_delay, guessFit, lw=1, label=r"Guess Fit")
    ax5[1].plot(time_delay, yfit, lw=1, label=r"SciPy Fit", color="orange")
    ax5[0].set_xlabel("Time delay [ns]")
    ax5[1].set_ylabel("Pressure [$Pa}$]")
    ax5[1].set_title("Pressure, Fit")
    ax5[0].set_title("Error in Fit")
    f5.suptitle("Sigmoid fit Pressure")
    ax5[1].legend()
    print("Sigmoid popt Pressures: " + str(popt))

    # Temperatures
    plt.close("all")
    temperatures_K = pressures_Pa*volumes_m3/n_moles/R_gasConst
    temperatures_K_err = temperatures_K*3*(err_diameter_um/diameter_um)
    fT, axT = plt.subplots()
    axT.errorbar(time_delay, temperatures_K, yerr=temperatures_K_err, marker=".", label='Temperature errors K',linestyle='none', color="red")
    axT.plot(time_delay, temperatures_K, marker="x", label='Temperature K',linestyle='none', color="blue")
    axT.set_xlabel("Time delay [ns]")
    axT.set_ylabel("Temperatures K")
    axT.set_title("Calculated Temperatures from Corrected Diameters")
    axT.legend()
    print("Last temp: " + str(temperatures_K[-1]))
    print("Last temp err: " + str(temperatures_K_err[-1]))

    # Doppler broadening
    kb_SI = 1.380649e-23
    T_K1 = 41000000.0 # K
    Ar_M_u = 39.948 #u
    Ar_M_kg = Ar_M_u * u_to_kg
    c = 3.0e8

    deltaFwhmDopplerShift = np.sqrt(8*kb_SI*T_K1*np.log(2)/Ar_M_kg/c/c)

    print("deltaFwhmDopplerShift" + str(deltaFwhmDopplerShift))

    # Show plots
    plt.close("all")
    plt.show()