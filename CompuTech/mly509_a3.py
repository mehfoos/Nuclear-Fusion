'''
Assignment3: Modified Sheath Equations
Rewritten based on feedback from assignment 1 and 2
'''

import matplotlib
from matplotlib import ticker
import matplotlib.pyplot as plt  # Plotting library
from scipy.integrate import odeint # SciPy ODE integration
from scipy.interpolate import interp1d # SciPy 1D interpolant
import numpy as np # numpy functions and constants

def sheathFunc(f, x, vs, L):
    # f is an array of all evolving variables
    phi = f[0] # First element is phi
    E = f[1] # Second element is E
    vi = f[2] # Third element is vi
    
    # Calculate intermediate results
    ne = np.exp(phi)
    ni = vs/vi

    # Calculate the derivatives
    dphi_by_dx = -E
    dE_by_dx = ni - ne
    dvi_by_dx = E/vi - vi/L
    
    # Return the derivatives in the same order as in the input f
    return [dphi_by_dx, dE_by_dx, dvi_by_dx]

def solve(x, vs, phi_Initial, E_Initial, L):
    # Calculate initial vi
    vi_Initial = np.sqrt(vs**2 - 2*phi_Initial)

    # Integrate using Sheath Function
    y = odeint(sheathFunc, [phi_Initial, E_Initial, vi_Initial], x, args = (vs, L))

    # Extract phi
    phi = y[:,0]

    # Extract vi
    vi = y[:,2]

    # Calculate j
    mi_by_me = 1840.0
    j = np.sqrt(mi_by_me/ (2.0*np.pi)) * np.exp(phi) - 1.0

    # Return the value of phi and j
    return phi, j, vi

if __name__ == "__main__":
    # Define all parameter values
    vs = 1.0
    phi_Initial = 0.0
    E_Initial = 0.001
    x = np.linspace(0, 40, 100)
    L_Array = np.logspace(-1, 4, num=6)
    vi_Wall_Array = np.array([]) # Can be pre-allocated, but doing this for clarity

    # Create figure explicitly to obtain references to axes
    fig, ax1 = plt.subplots()

    # Solve for each L
    for L in L_Array:
        # Solve for each L
        phi, j, vi = solve(x, vs, phi_Initial, E_Initial, L) # Solve to obtain j
        
        # Shift distance to set origin to wall distance
        xAsFunctionOfj = interp1d(j, x) # Get distance as an interpolant function of j
        xWall = xAsFunctionOfj(0) # Interpolate to find x where j is 0
        xShifted = x-xWall # Shifted distance coordinates such that j is alway 0 for distance 0

        # Plot vi against x, for each L
        ax1.plot(xShifted, vi, label = r"$\hat{L} =$"+ " {L}".format(L=L))

        # Store value of vs at wall
        viAsFunctionOfj = interp1d(j, vi) # Get vi as an interpolant function of j
        vi_Wall_Array = np.append(vi_Wall_Array, viAsFunctionOfj(0))
        print(vi_Wall_Array)


    # Inset plot and formatting
    left, bottom, width, height = [0.20, 0.5, 0.15, 0.35] # Fractions for inset plot
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(L_Array,vi_Wall_Array,'black')
    ax2.set_xlabel(r'Collision Length $\hat{L}$ [in $λ_D$]')
    ax2.set_ylabel(r'Vi at wall $\hat{L}$ [in $c_s$]')
    ax2.set_ylim([0.5,3])
    ax2.set_xlim([1.e-01,1.e+04])
    ax2.set_xscale('log')
    ax2.tick_params(bottom=True, top=True, right=True, left=True, direction="in", which="both")
    ax2.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1,numticks=10))

    # Main Plot Formatting
    ax1.grid(True)  # Add a background grid
    ax1.set_xlabel(r'Distance from wall [$\hat{x}$, Normalised to Debye Length $λ_D$] ')
    ax1.set_ylabel(r'Ion velocity [Normalised to ion sound speed $c_s$]')
    ax1.set_title(r'Assignment 3: Ion velocity at different collision lengths')
    ax1.legend(loc='upper right')
    ax1.set_ylim(0,5)
    ax1.set_xlim(-20,20)

    plt.show()

