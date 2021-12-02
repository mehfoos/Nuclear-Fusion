import numpy as np
import icf
import re
import matplotlib.pyplot as plt

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
    #ax1[0].plot(time_delay_raw, diameter_um, marker="x", label=r"Diameter [um]", color="blue", linestyle='none')
    ax1[0].errorbar(time_delay_raw, diameter_um_raw, yerr=err_diameter_um_raw,marker=".", label='Diameter errors [+-um]',linestyle='none', color="red")
    ax1[0].plot(time_delay_raw, diameter_um_raw, marker="x", label='Diameter [um]',linestyle='none', color="blue")
    ax1[0].set_xlabel("Time delay [ns]")
    ax1[0].set_ylabel("Diameters (raw) [um]")
    ax1[0].set_title("Raw Diameters, Errors")
    ax1[0].legend()

    # Volumes
    volumes_um3_raw = (4/3)*np.pi*(diameter_um_raw/2)**3 # um^3
    volumes_err_um3_raw = 3*(4/3)*np.pi*err_diameter_um_raw/diameter_um_raw
    volumes_m3_raw = (4/3)*np.pi*(diameter_um_raw*um_to_m/2)**3 # um^3
    volumes_err_m3_raw = 3*(4/3)*np.pi*err_diameter_um_raw/diameter_um_raw
    ax1[1].errorbar(time_delay_raw, volumes_um3_raw, yerr=volumes_err_um3_raw,marker=".", label='Volume errors [+-$um^{3}$]',linestyle='none', color="red")
    ax1[1].plot(time_delay_raw, volumes_um3_raw, marker="x", label='Volumes [$um^{3}$]',linestyle='none', color="blue")
    ax1[1].set_xlabel("Time delay [ns]")
    ax1[1].set_ylabel("Volumes (raw) [$um^{3}$]")
    ax1[1].set_title("Raw Volumes, Errors")
    ax1[1].legend()


    # Densities
    f2, ax2 = plt.subplots(2, sharex=True)
    densities_kgm3_raw = mass_initial/volumes_m3_raw
    densities_err_kgm3_raw = mass_initial*volumes_err_m3_raw
    ax2[0].errorbar(time_delay_raw, densities_kgm3_raw, yerr=densities_err_kgm3_raw,marker=".", label='Densities errors [$kg.m^{-3}$]',linestyle='none', color="red")
    ax2[0].plot(time_delay_raw, densities_err_kgm3_raw, marker="x", label='Densities [$kg.m^{-3}$]',linestyle='none', color="blue")
    ax2[0].set_xlabel("Time delay [ns]")
    ax2[0].set_ylabel("Densities (raw) [$kg.m^{-3}$]")
    ax2[0].set_title("Raw Densities, Errors")
    ax2[0].legend()

    # Pressures
    # Assuming Adiabatic
    gamma = 5.0/3.0
    adiabatic_PV_const = P_d2*(V_capsule_initial**gamma)
    pressures_Pa_raw = adiabatic_PV_const/(volumes_m3_raw**gamma)
    pressures_Pa_err_raw = gamma*adiabatic_PV_const*volumes_err_m3_raw
    ax2[1].errorbar(time_delay_raw, pressures_Pa_raw, yerr=pressures_Pa_err_raw,marker=".", label='Pressure errors [$Pa}$]',linestyle='none', color="red")
    ax2[1].plot(time_delay_raw, pressures_Pa_raw, marker="x", label='Pressures [$Pa$]',linestyle='none', color="blue")
    ax2[1].set_xlabel("Time delay [ns]")
    ax2[1].set_ylabel("Pressure (raw) [$Pa$]")
    ax2[1].set_title("Raw Pressures, Errors")
    ax2[1].legend()



    # Simple Results print
    print("V_capsule_initial: " + str(V_capsule_initial))
    print("n_moles: " + str(n_moles))
    print("density_initial: " + str(density_initial))
    print("number_density_initial: " + str(number_density_initial)),
    print("mass intial: ", str(m_d2*n_moles*N_a))


    # Show plots
    plt.show()