'''
Assignment2: Debye sheath conditions comparison
Rewritten based on model answers for Assignment1
'''

import matplotlib.pyplot as plt  # Plotting library
from scipy.integrate import odeint # SciPy ODE integration
from scipy.interpolate import interp1d # SciPy 1D interpolant
from numpy import linspace, pi, sqrt, exp # numpy functions and constants

def sheathFunc(f, x, vs):
    # f is an array of all evolving variables
    phi = f[0] # First element is phi
    E = f[1] # Second element is E
    
    # Calculate intermediate results
    ne = exp(phi)
    vi = sqrt(vs**2 - 2*phi)
    ni = vs/vi

    # Calculate the time derivatives
    dphi_by_dx = -E
    dE_by_dx = ni - ne
    
    # Return the time derivatives in the same order as in the input f
    return [dphi_by_dx, dE_by_dx]

def solve(x, vs, phi_Initial, E_Initial):
    # Integrate using Sheath Function
    y = odeint(sheathFunc, [phi_Initial, E_Initial], x, args = (vs,))

    # Extract phi
    phi = y[:,0]

    # Calculate j
    j = sqrt(1840.0/ (2.0*pi)) * exp(phi) - 1.0

    # Return the value of phi and j
    return phi, j

if __name__ == "__main__":
    # Define all parameter values
    vsTuple = (1, 1.5, 2) # Allows for adding additional vs
    phi_Initial = 0.0
    E_Initial = 0.001
    x = linspace(0, 40, 100)

    # Loop through each vs
    for vs in vsTuple:
        phi, j = solve(x, vs, phi_Initial, E_Initial) # Solve to obtain j
        xAsFunctionOfj = interp1d(j, x, assume_sorted = False) # Get distance as an interpolant function of j
        xWall = xAsFunctionOfj(0) # Interpolate to find x where j is 0
        xShifted = x-xWall # Shifted distance coordinates such that j is alway 0 for distance 0
        plt.plot(xShifted, j) # Plot

    plt.grid(True)  # Add a background grid
    plt.xlabel(r'Distance [$\hat{x}$, Debye Length $λ_D$] (Offset such that Current Density is 0 at 0 distance)')
    plt.ylabel(r'Current Density [Normalised to ion current density at wall $\hat{j}$]')
    plt.title(r'Assignment2: Comparison of Sheath Conditions')
    plt.legend([r'$\hat{V}_s$ = 1.0',r'$\hat{V}_s$ = 1.5',r'$\hat{V}_s$ = 2.0'])
    plt.ylim(-5,20)
    plt.xlim(-40,30)
    plt.show()